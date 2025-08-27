# Initial Design of DivBase 

## Status

Proposed

## Context

Genetic diversity research, such as population genomics, relies on the analysis and management of large, complex genetic variation datasets, typically stored in VCF ([Variant Call Format](https://samtools.github.io/hts-specs/VCFv4.5.pdf)) files. As studies scale up, effective data management becomes increasingly challenging. This difficulty largely arises from the fragmented nature of research data, which is often scattered across research group members, external collaborators, and government agencies. VCF files and associated sample metadata are rarely stored together, and reliable backups are often lacking, which hinders data analysis and future reuse. Without centralized infrastructure for data storage, versioning, and querying, the risks of meta(data) loss, duplication, and limited accessibility increase significantly. To address these issues, DivBase (Genetic Diversity Database) was proposed by Andreas Wallberg, a researcher from Uppsala University, during the 2024 DSN-EB inventory of needs.

DivBase will be a platform designed for the secure, medium-term management of VCF files and associated sample metadata, offering versioning and basic querying capabilities for non-human genetic variation data. It will support standard file formats that are compatible with the European Variation Archive (EVA), the international repository for genetic variation data. By addressing key data management challenges, DivBase aims to facilitate collaboration, enhance accessibility, and streamline the submission of data to long-term archives.

### Foreseen typical user journey in DivBase:

- A research group studying a species uploads their many VCF files (GB) to the DivBase object store.
- The most common setup will be, per research group: 1 studied species = 1 project = 1 bucket. New species = new bucket for that research group. 
- The research group has a metadata file or -database that contains metadata about these files. This file is uploaded to the bucket in addition to the VCFs.
- For follow-up analysis, researchers want to extract subsets from many files in the store. By using the metadata file and a pre-existing executable (bcftools), they can operate on the files in the object store to create a single (and much much smaller) "results" file.
- DivBase will have a webAPI that can take a query from the researcher, process the query (can take from seconds to hours) and save the output of the query (results file) to the object store. This results file can then be downloaded by the user. These queries should be treated like jobs and users should be able to find out the status of their jobs whenever. 
- The above bullets are iterated over during the research study, until the day the project ends and the study is published. At this point, the data stored in DivBase should be uploaded to a domain-specific end-repository for data sharing. DivBase should preferrably conform to the metadata standard of the key end-repository to facilitate the data publication.  

<p align="center">
  <img src="imgs/001-user-interactions.png" alt="User Interaction with DivBase" width="800">
  <br>
  <em>Fig. 1. Diagram showing how a user would interact with DivBase via a CLI tool, or via a web interface.</em>
</p>

### Core components
- S3 Object store with access control. 
- CLI tool to interact with the DivBase service. This will provide convenient ways to run queries on the files as well as interact with the files directly (e.g. upload/download). 
- WebAPI to accept query “jobs” and provide info about the status of each job. A user can use either the CLI tool or directly send requests to the API.  
- WebAPI talks to some kind of job management system that queues, runs and reports on the status of the jobs. 
- Software for prototype the query jobs. The query logic should be able to take the sample metadata file into account.
- Handling of metadata files. This includes user submitted sample metadata, but also system-generated metadata such as tracking bucket versioning states and indexing VCF dimensions for faster queries. 
- Frontend (not a priority), to manage projects/download API keys. 


## Decision

This initial system design is intended to cover all the architectural components needed to enable users to perform - at a basic level - the workflow outlined in the [Foreseen typical user journey in DivBase](#foreseen-typical-user-journey-in-divbase) section above.

### Core Technology Stack

- S3 object store. For testing and local dev, this will based on MinIO. In prod, this will likely be KTH NetApp S3 store (access have been requested, but is pending; [see sysadmin issue](https://github.com/ScilifelabDataCentre/sysadmin/issues/501)).  
- CLI tool built with [Typer](https://typer.tiangolo.com/).
- WebAPI built with [FastAPI](https://fastapi.tiangolo.com/).  
- WebAPI talks to job management system consisting of: [Celery](https://docs.celeryq.dev) worker(s); [RabbitMQ](https://www.rabbitmq.com/) queue/message broker; [Redis](https://redis.io/open-source/) results backend. Interaction with job management system using [Flower API](https://flower.readthedocs.io/en/latest/). 
- Jobs that consist of logic and syntax to run queries (here: subsets and filters) on sample metadata files (Pandas) and on VCF files using [bcftools](https://samtools.github.io/bcftools/howtos/index.html). 
- Metadata handling for the initial protoype will be done with flat files such as TSV, JSON, or YAML for proof-of-concept. Later on, a more robust solution is probably needed to handle I/O concurrency resulting from DivBase being a multi-user service with an asynchronous job system; most likely an ACID-compliant database will be needed.
- [Docker Compose](https://docs.docker.com/compose/) to orchestrate local dev environments with the above tech stack. The compose stack will also be used to spin up a separate testing environment for running [pytest](https://docs.pytest.org/en/stable/).

<p align="center">
  <img src="imgs/001-initial-design.png" alt="Initial Design of DivBase" width="1200">
  <br>
  <em>Fig. 2. Simplified overview of the data flow in core technology stack. With this approach the user can use the CLI tool to interact with either their projects S3 bucket (upload/download files) or create queries via the webAPI. The webAPI is responsible for auth and communicating with the job management system on the users behalf. Note that Redis and bcftools are not explicitly shown in the image: the former is part of the task history monitoring, the latter is part of the Celery workers.</em>
</p>


### System boundaries

The core functions of the DivBase is not public-facing but locked behind authentication/authorization. DivBase will require storage and compute, so the deployment on SciLifeLab hardware will - in line with scope the DDLS project - be limited to users that are affiliated with Swedish research institutions. Furthermore, DivBase is largely intended for use with unpublished data which also emphasizes the need for access control. Since auth is central to the service, it is likely a big architectural decision that will require its own ADR. 

DivBase is not a data publication platform. At the end of of research project, the users should submit their data and metadata to domain specific end-repositories for publication. This has implications for user roles, the lifespan of a given bucket, support/enforcement of metadata standards for key end-repositories (such as [EVA](https://www.ebi.ac.uk/eva/)). These should be covered in their own ADRs.

In line with the scope of the DDLS research area in Evolution and Biodiversity in which the team operates, DivBase will initially target research projects on non-model organisms. DivBase is not intended for use with sensitive data, such as human genetic variant data.

The investigations that have informed the current tech stack choices are primarily based on experiments in local dev environments using Docker Compose and toy datasets. We are yet to try the stack on the dev and production clusters, and can therefore not fully anticipate performance and implementation challenges on that hardware.

### Out of Scope: design choices that will be covered in separate ADRs

This ADR (001) is for the initial design of a DivBase prototype. More architectural design decisions will come, but will be handled in their own ADRs. Topics that are already foreseen include:

- Auth and user role management

- Metadata schema

- Query syntax

- Handling of metadata files (with a database engine?)

### Why these choices

The guiding principle in choosing the proposed tech stack was to keep the service as lightweight as possible while allowing us to demonstrate proof-of-concept for the basic functionalities of the service. As discussed above, DivBase needs an object store, a query engine, a job management system, user authentication, a webAPI, and a CLI. We believe that the choices fulfill that while also aligning to the DC paved-paths tech whenever possible. 



For readers that are unfamiliar with working with VCF files the importance of bcftools it might be evident from the list of the core tech stack, but its role as the core query engine in DivBase cannot be understaded. Bcftools is an very popular tool with a established syntax. By writing wrapper logic for running bcftools commands, the team will save a substantial amount of work compared to writing a our own custom VCF query engine. The intention is to use bcftools syntax as much as possible for DivBase to benefit from the fact that most researchers working with genetic diversity data analysis are likely to have at least encountered bcftools, if not even being proficient users. 


## Alternatives considered

- **Using Django (/flask) instead of FastAPI:** Not arguing one is better than the other. Due to the low priority of having a frontend and that the frontend is expected to be quite limited, we have decided to try out FastAPI initially. The admin panel that comes with Django would be one benefit of using Django instead (make viewing users/projects easier).   

- **Kueue or SLURM instead of Celery**. Both the devs and the pilot user group have long experience of using SLURM in HPC environmentents and imagine a user experience similar to that of SLURM (queue jobs, check status, position in queue, job logs stored for later inspection). While there are no paved-paths entry for job mangement systems at DC, Celery is already used at DC for Serve. As opposed to the listed alternatives, Celery is written in Python, which is the default language at DC. Should we for some reason change to Django as per the above point, it could be beneficial that Celery has out-of-the-box support for Django. 

- **RabbitMQ also as job system results backend instead of Redis**. While RabbitMQ is primiarly designed to be a job system message broker, it is able to take on the role of job systems result backend as well. Doing so would result in one less component to maintain in the tech stack. However, Redis was chosen as a results backend for its performance, ability to store persistent logs, and its excellent integration with the Flower monitoring system.

- **The proof-of-concept demo for DivBase use flat files for storing sample and system metadata**. The async, concurrent job system and the fact that multiple users might simultaneously work on the same files in the bucket might prove a challenge. Later on, we might need an ACID-compliant database such as PostgreSQL to more reliably store the metadata. However, this should be the topic of its own ADR.

## Consequences

Adopting the proposed tech stack for developing the inital prototype for Divbase is foreseen to have the following consequences:

### Positive

- The proposed tech stack will ensure that the core functionalities of DivBase will be implementable and will, as such, be a major step towards what will eventually be the DivBase MVP.

- The choice of bcftools ensures that all standard VCF subsets and filtering operations are supported by DivBase. This means that the team can focus on developing logic for query functionalities not supported by bcftools: specifically coupling queries to a sample metadata file with the bcftools operations to create more powerful query "pipelines".

- At the time of writing, two team members have research experience working with VCF files and bcftools, and thus have a good understanding of use-cases, specifications, and quirks with the data format and its applications in different research projects.

### Risks & Challenges

- **Team learning curve** The team members have experience of some of the components in tech stack, but most of them are new to the team. Basic implementation does not seem to be an issue, but it would be valuable if other staff at DC or SciLifeLab that could help with advice on best practices and optimization at some point during the development. 

- **Creating a robust system that can handle concurrency demands**. Concurrent operations in a multi-user service will need to be protected against file corruption due to parallel I/O operations. What if a file update is not correctly stored in the bucket? How will that affect the system, the users, their collaborators, and the quality of the research project? This will need to specifically addressed during the development.

- **Performance and implementation challenges on the cluster hardware is currently an unknown** Since we do not yet have been given access to the KTH NetApp S3 store, the proposed tech stack choice has been informed by investigations done using local dev environments. We are aware of the possibility that once we get access to the production hardware and can get an idea of query performance, there is a risk that we need to reconsider some tech choices.