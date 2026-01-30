# Bcftools Celery Task Logic

User submitted queries to check out VCF data from DivBase are run as Celery tasks as defined in `tasks.bcftools_query` in `worker/tasks.py`. This page contains an overview of how the logic is implemented.

The task logic can roughly be divided in DivBase system operations, and bcftools operations. The system operations are used to find which VCF files to use, check if the are compatible with the bcftools logic using the VCF Dimensions data stored for the project in the postgreSQL database, followed by download of the VCF files to the worker that will perform the bcftools operations. The bcftools operations are designed so that each input VCF file is subset seperatelly, intermediate results are saved to temp files, and, in the end, all temp files are combined to a single results file. The input files needs to follow the constrains for bcftools in order to be sucessfully combined (summarised for the DivBase use case in [bcftools Celery Task Constraints](bcftools_task_constraints.md)). This is checked in the systems operations step before downloading the VCF files to the worker. If the check does not pass, the task exits so that no resources are wasted on jobs that are known to not work.

In the subheadings below are flow charts that describe how the bcftools task logic is implemented in DivBase. The are two versions: one describing the full flow, and one describing a simplified version of the flow.

## Diagram of Full bcftools Task Logic Flow

![bcftools Celery Task Logic Full Flow](../assets/diagrams/bcftools_task_logic_flow_full.svg)

Figure 1: Flow chart describing the full signal flow of the bcftools Celery task.

## Diagram of Simplified bcftools Task Logic Flow

![bcftools Celery Task Logic Simplified Flow](../assets/diagrams/bcftools_task_logic_flow_simplified.svg)

Figure 2: Flow chart describing the simplified signal flow of the bcftools Celery task.
