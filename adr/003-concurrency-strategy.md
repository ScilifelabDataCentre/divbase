# 1. ADR 003: Initial DivBase concurrency strategy

## 2. Status
Proposed

## 3. Context
DivBase will be a multi-user service with an async job system and and S3 data store. This means that the system will need to be designed to robustly handle concurrency demands to avoid data corruption and erroneous results. This ADR specifically handles the concurrency challenges that arise due to the asynchronous nature of the job system and it's need to read and write files from S3 buckets. 

As described in [ADR-001](adr/001-initial-system-design.md), the core of the DivBase job mangement system consists of Celery workers, a RabbitMQ message broker, and a persistent results backend. Of note for the discussion below is that [Celery Task-ID are UUIDs](https://docs.celeryq.dev/en/latest/reference/celery.app.task.html) and are generated upon task submission.

## 4. Decision

The decision consists of two parts: a general strategy that can be applied to all Celery tasks that handle files in S3 (Section 4.1), and specific details for individual tasks (Sections 4.2.-4.3). Failure scenarios are covered in Section 4.4.

### 4.1. General strategy that applies to all tasks
To ensure robust operations, every user interaction that goes through the job system need to do the following in sequential order

1. Fetch file version IDs from the bucket of choice for all files in the bucket (or preferably, only the files that are subject to the task, whenever that is trivial knowledge)
2. Submit the task to the job manager queue;  include the file version IDs in addition to other parameters needed for that specific task
3. Idle worker eventually picks up task from queue and checks with S3 that the file version IDs exist in the bucket (e.g. that they have not been deleted or the service somehow has been disrupted).
4. Send message to user: 
	- If the file version IDs were successfully verified, tell the user that the task will use these specific versions of the files (= version ID captured in step 1). For transparency and information secutiry, let the user know whether or not these are the latest files (for the case that files in the bucket have been updated in the time window between the task submission and task start).
	- If the file version IDs failed verification, return an error message to the user informing them of the issue. If the files have been deleted from S3, tell the user to ensure that files are in S3. If the system cannot connect to S3, tell them that.
5. Download the verified versions of files from S3 to worker. 
	- It is possible to control worker concurrency in Celery and in Kubernetes. The decision is to use Celery concurrency=1 together with K8s horizontal scaling. This means that each worker pod only handles one job at a time; if the queue is long and cluster resources are available, new pods can be scaled up to accept more jobs.
	- This has several benefits: simpler file management (no filename collision risk - e.g. if one worker would handle >1 job and they would happend to have the same filenames); resource isolation (job or pod failure doesn't affect other jobs); separated logging and resource monitoring (pod logs become per-job logs); simpler temp file garbage collection (identifying which files belong to which job).
    - File checksums will be verified after all downloads have completed. Which checksum algorithm that will be used will be decided based on the ETAG algoritm used at KTH NetApp.
	- Depending on S3 perforance and bandwith limitations (see [ADR-001: Performance validation](adr/001-initial-system-design.md#performance-validation)), we might want to cache the input VCF files in the workers to enable reuse by subsequent jobs. If so, caching will be done by appending their S3 version ID to the filename. First, the worker checks if the specified file version already exists in its mounted volume; if not, it downloads it to a shared persistent volume claim (ReadWriteMany PVC); caching is also shared between workers in this way. A Time-To-Live based cache eviction will be used: a cron job will be run that delete files older than 24 hours. If we decide to not use input file caching in the workers, input files will be handled by a context manager to ensure that they are deleted after job completion or on job failure.
6. Return results of query and clean up temp files. This is different from task to task, but can involve writing to S3 or to a database table. Redis lock will be used to coordinate access to shared resources to avoid race conditions and data corruption. The `python-redis-lock` library and Redis SETNX-based locks will be the first choice for implementing this.

Task-specific details will be addressed below.

### 4.2. Strategy specific to metadata query tasks

#### 4.2.1. Sample metadata

- Worker needs to read file from S3: Yes
- Worker needs to write file to S3: No

**Context:** 

The sample metadata file is user-provided and is not updated by the system.

**Addition to general sequence in Section 4.1.:**
- 5. Download the verified version of the sample metadata from S3 to worker. Filename is specified by user when submitting the job (or, if omitted by user, a default name is expected. E.g. `sample_metadata.tsv`). Upon download to worker container, add Celery Task-ID to the file name.
- 6. Upon completion of tasks that require the sample metadata file, delete (Task-ID hashed) metadata file from the worker.

#### 4.2.2. Technical metadata

- Worker needs to read from database table: Yes
- Worker needs to write to database table: Yes

**Context:**

A technical metadata file that stores key information about each VCF file in the bucket is also planned. This file will, for each VCF file, store the samples, scaffold names, and number of variants contained in each file, and the file version ID (for verification). By storing this information in a metadata file, the system will only need to download each VCF file once (for each version of the file submitted to S3) and log these parameters, instead of having to download each VCF every time this information is needed to make a logic decision in the system. This will save on resources and be more time efficient for users. For example, there are invalid sample set combinations for `bcftools merge` or `bcftools concat`, and having a technical metadata file allows us to implement _a priori_ checks to see if a bcftools-dependent task is feasible, and ensure an early-exit for the task if needed.    

Unlike the sample metadata file, the technical metadata is intended to be updated by the backend and should never be interacted with by the users, except to perhaps manually schedule a refresh. As VCF files change in the bucket, this can be subject to several concurrent requests. Capturing the technical metadata requires to parse the whole VCF file (not just the headers, since there are row-specific information, such as the scaffold names, that needs to be stored) and thus that should be done by the workers. Updates to the technical metadata will be scheduled by the API each time a new version of a VCF file is uploaded to the bucket. As a failsafe, tasks that deal with VCF files will check if the version IDs of the files indexed in the latest version of the technical metadata are up to date with the files in the bucket; if not, send an error to the users asking them to run the update command.

Since the technical metadata is intended to be used by the backend and not by the users, it will be stored in a PostgreSQL table. This is an ACID-compliant database management system that is compatible SQL SELECT queries, indexing for fast lookups, and Redis locks. Another benefit is that the CloudNativePG on the KTH cluster can be used for automated backups. For continuity with the auth database strategy in [ADR-002: ](adr/002-API-design.md), database models will be defined with SQLalchemy 2. Redis locks will be applied to the given table rows when a update job is queued, and released once the computations are complete and the metadata updated in the row. The planned schema looks like this:

```bash
class VCFMetadata(Base):
	#Inherits Base as defined in ADR-002 (UUID, created_at, updated_at)

    __tablename__ = "vcf_metadata"

    vcf_file_s3_key = Column(String, primary_key=True, index=True) # the unique path or key for the file in the S3 bucket.
    project_id = Column(UUID, ForeignKey('projects.id'), index=True)
    s3_version_id = Column(String, nullable=False, index=True)

    # VCF dimensions
    samples: Mapped[list[str]] = mapped_column(ARRAY(String), index=True)  # List of sample names
    scaffolds: Mapped[list[str]] = mapped_column(ARRAY(String))  # List of Chromosome/scaffold names
    variant_count: Mapped[int] = mapped_column(BigInteger)
    sample_count: Mapped[int] = mapped_column(Integer)

    # Metadata
    file_size_bytes: Mapped[int] = mapped_column(BigInteger)
    indexed_at: Mapped[DateTime] = mapped_column(DateTime, default=func.now()) # Timestamp from last update to row
```

**Addition to general sequence in Section 4.1.:**
- 3. Upon task start, apply a redis lock to the row in the technical metadata table. If errors are raised duing the task, release the redis lock. The redis locks should also have an expiry time to avoid deadlocks, at 2-3 times the worse-case time it takes to perform the operation. The runtime for updating the technical metadata scales proportionally to number of variants in all VCF files in the bucket. In a worst-case time, the deadlock should need no longer than 1-2 h, but this needs to be tested.
- 6. If the task includes making updates to the technical metadata database, write the new version to the table. After writing, release the lock.


### 4.3. Strategy specific to VCF query tasks

- Worker needs to read file from S3: Yes
- Worker needs to write file to S3: Yes, a results file; this is a subset of the input VCFs.

**Context:**

Queries on VCF files in DivBase relies on bcftools to subset VCF files. If a query concerns multiple VCF files, the task logic handles subsetting each file individually, and once that is done, combines all the subsets to a result in the form of a single VCF file. 

**Addition to general sequence in Section 4.1.:**
- Since the task logic for running bcftools-dependent tasks will need to perform checks to see if the query is feasible, the worker will need to interact with the technical metadata. Thus, the additions listed in 4.2.2. apply.  
- The VCF queries can optionally take a sample metadata query as input. In that case, the additions listed in 4.2.1. apply. 
- 6. The output of the query is a subset VCF file that is uploaded to the given S3 bucket. This file should include the Celery Task-ID in the filename for uniqueness and traceability.

### 4.4. Failure scenarios

- Redis lock expires before job completes? 

Lock renewal can be used when a job takes longer time than the initial lock expiry time. The challenge is to not renew locks to jobs that are stuck in futile loops. Tasks will log each operational step, and logic can be implemented to check if the job has progressed since the last lock application or renewal and only then renew the lock. The challenge with this approach is that bcftools is not very verbose during processing; on the other hand, any error raised by bcftools will terminate the task. A hard upper limit for job duration should be implemented so that locks cannot be renewed longer than that. Tests will be needed to understand how an upper limit for heavy jobs will look like, but with lock renewal in place, it might need to be as high as 24h.

- S3 version verification fails mid-job? 

Version verification will run at the start of the job and used to download the specific the input files versions from S3. If verification fails or an error is raised, the job will exit. A checksum verification will be run after completed download to verify file integrity after transfer. The filenames of the downloaded files will be appended with the version ID as an additional safeguard. All processing steps, e.g. with bcftools, will expect the version ID appended filenames.

- File cleanup if worker crashes during file download?

While context managers and input file caching will handle many cleanup cases, they will not _per se_ handle cleanup after worker pod crashes. For input files, cache eviction jobs is considered sufficient to purge input files after a given time, as described earlier in this ADR. For temp file cleanup, a dedicated clean up job will scheduled that identifies dangling temp files based on the Celery Task ID appended to their filenames, their timestamp, and the RabbitMQ queue. Files with Task IDs that belong to tasks that are no longer running can then be safely purged from the worker volume.

- Loss of technical metadata?

 Not a major risk since: a) the database will be backed up using CNPG, and b) the technical metadata can be regenerated by reading the VCF files again.




## 5. Consequences

Positive:

- The concurrency strategy adds robustness to DivBase and is essential since multi-users and async jobs are involved.
- Ensures that the exact file versions in the specified bucket at the time of job submission are the files that will be used. Communicates this to the user for the sake of reproducbility.
- Strategy can be implemented for existing architechture. No need to add any addition services/components other than once already considered (Redis, PostgreSQL).

Negative:

- While crucial for the service to function, the concurrency strategy adds complexity that could complicate maintenance.


## 6. Alternatives Considered:

- Redis lock and not S3 conditional writes:

The benefit of S3 conditional writes is that is already available through the S3 object store. The downside is that it can only lock files in S3, whereas Redis lock can, in addition to files in S3, also act on files in volumes and on database tables. A downside of Redis lock is that it requires that a Redis service is running. However, a Redis instance will be used in the job system results backend, and can thus also be used for locking.

- Store technical metadata in S3 in a YAML or JSON instead of a database table:

For early proof-of-concept, using a YAML file was a good way to scope out the needs of the technical metadata file. However, given that users do not need interact directly with it, the decision is to store the technial metadata in a database table instead. Storing this in a YAML file in S3 has several downsides: S3 transfer needed (production S3 is outside the DivBase k8s namespace network), YAML deserialization and parsing overhead, and it has none of the built-in SQL database benefits (query complexity, column indexing, ACID-compliance, etc.).

- Using Celery Task-ID as hash in temp filenames in worker instead of using file version ID:

Hashing of all files (input files, temp files, results files) with Celery Task-ID removes the possibility to reuse the same input files between tasks. VCF files are potentially large in size, and it would decrese the transfer load from S3 to workers if VCF files could be stored temporarily for a limited time on a worker-mounted PVC instead of transferred once for each task. If the S3 version ID is instead used as a hash, it opens up the possibility to reuse files.

- Input file cache eviction strategy: LRU (Least Recently Used) instead of TTL (Time-To-Live)

The LRU strategy is based on determining which files to remove files when space is needed based on longest time since last accessed. This might be more suitable for DivBase in the long run, since it is predicted that users will use the service in bursts rather than on an everyday basis. The idea is to start with TTL and monitoring the use and investigate the following: how often is the same file needed by the user queries? is it more efficient to store that file for a longer time rather than downloading it again everytime it is needed after TTL eviction?