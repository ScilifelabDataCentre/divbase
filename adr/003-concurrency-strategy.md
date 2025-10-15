# 1. ADR 003: Initial DivBase concurrency strategy

## 2. Status
Proposed

## 3. Context
DivBase will be a multi-user service with an async job system and and S3 data store. This means that the system will need to be designed to robustly handle concurrency demands to avoid data corruption and erroneous results. This ADR specifically handles the concurrency challenges that arise due to the asynchronous nature of the job system and it's need to read and write files from S3 buckets. 

As described in [ADR-001](adr/001-initial-system-design.md), the core of the DivBase job mangement system consists of Celery workers, a RabbitMQ message broker, and a persistent results backend. Of note for the discussion below is that [Celery Task-ID are UUIDs](https://docs.celeryq.dev/en/latest/reference/celery.app.task.html) and are generated upon task submission.

## 4. Decision

The decision consists of two parts: a general strategy that can be applied to all Celery tasks that handle files in S3 (Section 4.1), and specific details for individual tasks (Sections 4.2.-4.3).

### 4.1. General strategy that applies to all tasks
To ensure robust operations, every user interaction that goes through the job system need to do the following in sequential order

1. Fetch file version IDs from the bucket of choice for all files in the bucket (or preferably, only the files that are subject to the task, whenever that is trivial knowledge)
2. Submit the task to the job manager queue;  include the file version IDs in addition to other parameters needed for that specific task
3. Idle worker eventually picks up task from queue and checks with S3 that the file version IDs exist in the bucket (e.g. that they have not been deleted or the service somehow has been disrupted).
4. Send message to user: 
	- If the file version IDs were successfully verified, tell the user that the task will use these specific versions of the files (= version ID captured in step 1). For transparency and information secutiry, let the user know whether or not these are the latest files (for the case that files in the bucket have been updated in the time window between the task submission and task start).
	- If the file version IDs failed verification, return an error message to the user informing them of the issue. If the files have been deleted from S3, tell the user to ensure that files are in S3. If the system cannot connect to S3, tell them that.
5. Download the verified versions of files from S3 to worker. 
    - TO BE DISCUSSED: celery workers can have concurrency (unless specified, default= 4), but it is also possible to set worker concurrency to 1 and handle scaling in kubernetes. If we decide to use celery workers with concurrency >1, unique temp file names becomes crucial. For example, a worker with concurrency=2 picks up two jobs that each act on a different bucket, but happen to have files that share name but not content (e.g. sample_metadata.tsv), errors are likely to occur. Appending the Celery Task-ID to the filename and parsing for that in the tasks is one solution to this issue. 
6. Return results of query and clean up temp files. This is different from task to task, but can involve writing to S3. Task-specific details will be addressed below.

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

- Worker needs to read file from S3: Yes
- Worker needs to write file to S3: Yes

**Context:**

A technical metadata file that stores key information about each VCF file in the bucket is also planned. This file will, for each VCF file, store the samples, scaffold names, and number of variants contained in each file, and the file version ID (for verification). By storing this information in a metadata file, the system will only need to download each VCF file once (for each version of the file submitted to S3) and log these parameters, instead of having to download each VCF every time this information is needed to make a logic decision in the system. This will save on resources and be more time efficient for users. For example, there are invalid sample set combinations for `bcftools merge` or `bcftools concat`, and having a technical metadata file allows us to implement _a priori_ checks to see if a bcftools-dependent task is feasible, and ensure an early-exit for the task if needed.    

Unlike the sample metadata file, this file is intendedn to be updated by the backend and should never be interacted with by the users, except to perhaps manually schedule a refresh. As VCF files change in the bucket, this file can be subject to several concurrent requests. Capturing the metadata requires to parse the whole VCF file (not just the headers, since there are row-specific information, such as the scaffold names, that needs to be stored) and thus that should be done by the workers. Updates to this file will be scheduled by the API each time a new version of a VCF file is uploaded to the bucket. As a failsafe, tasks that deal with VCF files will check if the version IDs of the files indexed in the latest version of the metadata file are up to date with the files in the bucket; if not, send an error to the users asking them to run the update command.

(TO BE DISCUSSED: this technical metadata file can be stored in S3 along with the VCF files, but since the users don't need to directly interact with it, it could as well be stored in a worker-mounted PVC, or in a PostgreSQL table. In either way, redis lock would be convenient since it can handle all three of those use-cases. Loss of this file is not a major risk, since it can be regenerated by reading the VCF files again)


**Addition to general sequence in Section 4.1.:**
- 3. Upon task start, apply a redis lock to the file in the bucket. If errors are raised duing the task, release the redis lock. The redis locks should also have an expiry time to avoid deadlocks, say 24 h.
- 5. Download the verified version of the technical metadata from S3 to worker. Upon download to worker container, add Celery Task-ID to the file name.
- 6. If the task includes making updates to the technical metadata file, write the new version to the S3 bucket. After writing, release the lock. 


### 4.3. Strategy specific to VCF query tasks

- Worker needs to read file from S3: Yes
- Worker needs to write file to S3: Yes, a results file; this is a subset of the input VCFs.

**Context:**

Queries on VCF files in DivBase relies on bcftools to subset VCF files. If a query concerns multiple VCF files, the task logic handles subsetting each file individually, and once that is done, combines all the subsets to a result in the form of a single VCF file. 

**Addition to general sequence in Section 4.1.:**
- Since the task logic for running bcftools-dependent tasks will need to perform checks to see if the query is feasible, the worker will need to interact with the technical metadata. Thus, the additions listed in 4.2.2. apply.  
- The VCF queries can optionally take a sample metadata query as input. In that case, the additions listed in 4.2.1. apply. 
- 6. The output of the query is a subset VCF file that is uploaded to the given S3 bucket. This file should include the Celery Task-ID in the filename for uniqueness and traceability.


## 5. Consequences

Positive:

- The concurrency strategy adds robustness to DivBase and is essential since multi-users and async jobs are involved.
- Ensures that the exact file versions in the specified bucket at the time of job submission are the files that will be used. Communicates this to the user for the sake of reproducbility.
- Strategy can be implemented for existing architechture. No need to add any addition services/components other than once already considered (Redis, PostgreSQL).

Negative:

- While crucial for the service to function, the concurrency strategy adds complexity that could complicate maintenance.


## 6. Alternatives Considered:

- Redis lock and not S3 conditional writes:

The benefit of S3 conditional writes is that is already available through the S3 object store. The downside is that it can only lock files in S3, whereas Redis lock can act on files in volumes and on database tables. If we ever should decide that the technical metadata should not be stored in the S3 bucket (see next bullet below), using Redis lock from the start will simplify the migration to a non-S3 destination. A downside of Redis lock is that it requires that a Redis service is running. If we decide to use Redis for the job system results backend (decision currently pending discussion in ADR001), the same instance can also be used for locking.

- Store technical metadata in a database table instead of a YAML or JSON:

For early proof-of-concept, using a YAML file was a good way to scope out the needs of the technical metadata file. However, given that users do not need interact directly with it, it would make sense to store the data in a database table. PostgreSQL would be a natural choice since is it already used in the architechture. The technical metadata file does not gain anything from being versioned as long as it is kept up to date with the VCF files in the bucket; this could speak towards using a db table instead. Likewise, ACID compliance would be nice to have in the concurrency solution for this file. Potential downsides could be: added burden on the database instance, overkill solution.

- Using Celery Task-ID as hash in temp filenames in worker instead of using file version ID:

Hashing of files with Celery Task-ID removes the possibility to reuse the same input files between tasks. VCF files are potentially large in size, and it would decrese the transfer load from S3 to workers if VCF files could be stored temporarily for a limited time on a worker-mounted PVC instead of transferred once for each task. If the S3 version ID is instead used as a hash, it would open up the possibility to reuse files.