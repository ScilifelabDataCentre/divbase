# Resource Usage Metrics: CPU vs Memory in DivBase

This document explains how resource usage metrics are tracked and reported for Celery tasks and bcftools subprocesses in DivBase, with a focus on the differences between CPU and memory accounting.

DivBase uses `bcftools` for operations that act directly on VCF files. `bcftools` is a compiled binary and is considered to be a very efficient way to process VCF files. When measuring any resource metric for a DivBase task, it’s thus relevant to distinguish between the time spent in the core VCF processing step (handled by `bcftools`) and the time spent in supporting operations (the DivBase overhead). While bcftools sets the lower bound for how fast VCF processing can be, the overall task duration also depends on DivBase-specific steps such as downloading files from S3, checking data compatibility, and managing metadata. Ideally, the DivBase overhead should be as small as possible, giving users a performance similar to that of just running `bcftools` on the files.

This document was written with the Celery tasks that call on `bcftools` in mind. These are the most resource intensive tasks in DivBase and include e.g. downloading of VCF files from the S3 object store to the Celery workers and running of `bcftools`. The text still applies to any Celery task in DivBase: if `bcftools` is not used in a task the element of the calculations will be 0.

## Wall time: total time elapsed from process start to end

Wall time (also known as ["wall clock time", "real elapsed time", and other names](https://en.wikipedia.org/wiki/Elapsed_real_time)) is the total time elapsed from when a process started to when it ended. Wall time includes all time spent on CPU operations, I/O operations, and any periods of idle waiting. From a DivBase user's point-of-view, this is how long time it took from a task to run from start to finish. The wall time value is especially relevant for the user experience, but does not give any detailed hints on what part of the underlying process(es) that can be optimized.

As will be discussed in the section on CPU time below, `bcftools` will be run in one or several subprocesses. Wall time is measured for the Python proccess and for the `bcftools` subprocesses with the `time` python library. The wall time for the Python process is measured from the start to the finish of the task, meaning that:

$$
    \text{Total wall time for a DivBase task} = \text{Python process wall time}
$$

The wall time is also measured seperatelly for the time it takes to run the `bcftools` calls of the task. Since the wall time for total task and the `bcftools` subprocess are measured, the DivBase task overhead can be calculated as:

$$
 \text{DivBase overhead wall time} = \text{Python process wall time} - \sum_{i=1}^N \text{bcftools subprocess}_i \text{ wall time}
$$

To optimize wall time, other metrics such as CPU time and memory usage need to be considered.

## CPU Time: Additive Across Processes

[CPU time](https://en.wikipedia.org/wiki/CPU_time) is the time that the CPU has been actively working on process. Compared to wall time, CPU time is always less or equal to wall time (a process that executes without any breaks will have CPU time equal to wall time; if the process is waiting for I/O or is idle, CPU time will be less than wall time).

CPU time is an additive metric: it is the total amount of time that the CPU spends on executing a process, including distribution across multiple CPU cores (if available). Furthermore, the tasks in DivBase that use `bcftools` need to spawn subprocesses since `bcftools` does not have a Python API and is run as a complied C binary. Every time a `bcftools` command is needed during a DivBase task, the `subprocess` Python library is used to spawn a new process with it unique PID that consumes CPU time independently of the other processes. The Python process waits for each process to finish before continuing processing its own instructions, and as such the Python process can be seen as a master process for all the processes that are spawned and run during the task. The `bcftools` subprocesses use `subprocess.Popen` and `proc.wait()` instead of `subprocess.run` so that `psutlis` can be used to capture CPU and RAM usage during the subprocesses. The monitoring code is setup to measure each task-realted process indidivually: the Python process that runs the Celery task (including all operations that use Python libraries, such as `boto3` for interacting with the S3 object store), and the `bcftools` subprocesses (collected individually but stored as a accumulated total since there is little room for optimization in the DivBase of how `bcftools` runs during the individual subprocesses). This means that:

$$
\text{Total CPU time for a DivBase task} = \text{Python process CPU time} + \sum_{i=1}^N \text{bcftools subprocess}_i \text{ CPU time}
$$

Because all DivBase-specific operations (such as file downloads, metadata checks, and orchestration) are performed in the main Python process, and bcftools subprocesses are measured separately, the DivBase task overhead for CPU time is simply the CPU time of the Python process:

$$
 \text{DivBase overhead CPU time} = \text{Python process CPU time}
$$

Note that this is different from how wall time is measured in DivBase: task start to finish.

### CPU resource considerations for Kubernetes (k8s) deployment

Resource specifications in k8s is enforced per container (not per pod, the pod resouce usage is the sum of all its containers). CPU is measured in cores: 1 = one full core, 500m (millicore) is half a core, etc. Important to know is that when a container exceedes its CPU limitations in k8s, it is throttled but not killed; exceeding memory limitations result in containers being kill, but more on that in the memory section below. Essentially, the CPU requests and limits can be used to control wall time by:

$$
\text{Wall time (seconds)} = \frac{\text{Total CPU time (seconds)}}{\text{Allocated CPU cores}}
$$

**Example:**

$$
\frac{120\ \text{CPU-seconds}}{1\ \text{core}} = 120\ \text{seconds (Wall time)}
$$

To configure this for the k8s deployment use `requests` for the minimum guaranteed CPU resources and `limits` for the maximum allowed CPU resources. Thus, `limits` dictate the fastest possible wall time. Kustomize manifest example:

```
resources:
    requests:
        cpu: "200m"
    limits:
        cpu: "700m"
```

## Memory Usage: Not Additive Across Processes

TODO

## Implementation

TODO
