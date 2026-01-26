# Monitoring: Custom metrics collection for Celery worker

This document explains how resource usage metrics are tracked and reported for the Celery Workers in DivBase. Using a custom `prometheus-client` metrics server, wall time, CPU time, and Memory Usage can be collected on a per-task basis (togglable with the `ENABLE_WORKER_METRICS` environment variable). The metrics are collected with custom code written for DivBase, and this document collects the definition and rationale behind the units and calculations used to capture the metrics. This is specifically for the custom per-task resource monitoring; it is possible to set up system/cluster wide resource monitoring in addition to this (e.g. `cAdvisor`, `node-exporter` etc.)

This document was written with the Celery tasks that call on `bcftools` in mind. These are the most resource intensive tasks in DivBase and include e.g. downloading of VCF files from the S3 object store to the Celery workers and running of `bcftools`. The text still applies to any Celery task in DivBase: if `bcftools` is not used in a task the element of the calculations will be 0.

DivBase uses `bcftools` for operations that act directly on VCF files. `bcftools` is a compiled binary and is considered to be a very efficient way to process VCF files. When measuring any resource metric for a DivBase task, it’s thus relevant to distinguish between the time spent in the core VCF processing step (handled by `bcftools`) and the time spent in supporting operations (the DivBase overhead). While bcftools sets the lower bound for how fast VCF processing can be, the overall task duration also depends on DivBase-specific steps such as downloading files from S3, checking data compatibility, and managing metadata. Ideally, the DivBase overhead should be as small as possible, giving users a performance similar to that of just running `bcftools` on the files.

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

To configure this for the k8s deployment use `requests` for the minimum guaranteed CPU resources and `limits` for the maximum allowed CPU resources. Thus, `limits` dictate the fastest possible wall time.

```
# Kustomize manifest example

resources:
    requests:
        cpu: "200m"
    limits:
        cpu: "700m"
```

## Memory Usage: Not Additive Across Processes

Memory usage monitoring in DivBase is based on measuring RSS (Resident Set Size). How Linux system use memory is a much bigger topic than this document can ambition to describe, but in short, there is [RSS and VSZ (virtual memory)](https://stackoverflow.com/questions/7880784/what-is-rss-and-vsz-in-linux-memory-management). RSS is the memory allocation of a process the in physical memory (RAM). Unlike CPU time, RSS Memory usage is not additive across processes because operating systems allow processes to share memory regions, as for instance discussed in this forum [thread](https://stackoverflow.com/questions/131303/how-can-i-measure-the-actual-memory-usage-of-an-application-or-process) and in this [thread] (<https://unix.stackexchange.com/questions/34795/correctly-determining-memory-usage-in-linux>). See also: [Kerrisk, M. (2010). The Linux programming interface: A Linux and UNIX system programming handbook. No Starch Press](https://www.man7.org/tlpi/). There are other tools and methods to monitor memory, but RSS provides a decent estimate for the needs of DivBase optimisation.

Monitoring RSS memory usage is done per process. This means that for DivBase, the Python process that runs the task has its own memory usage, and each `bcftools` subprocess have their own memory usages. Note that although the Python process waits for the `bcftools` subprocess to finish, the memory used by the `bcftools` subprocesses is not reflected in the memory usage of the Python process during that time.

For DivBase, two memory metrics are calculated: Average memory usage (RSS, bytes) and Peak memory usage (RSS, bytes). The average is used to find the baseline RSS memory usage for the processes running during a task, and the peak is the highest RSS usage that was observed. These are used to set up the kubernetes memory request and limits, as described later below. To collect this data, the RSS memory usage of the Python process and each `bcftools` subprocess is sampled at frequent intervals and stored in a data pool. The average and peak RSS is then calculated for each process based on these data pools. Note that total cumulative memory is not tracked since RSS is not additive.

### Memory resource considerations for Kubernetes (k8s) deployment

As described above in the CPU section, resource specifications for memory in Kubernetes are enforced per container and not per pod. Memory is measured in bytes, and can be specified using units like Mi (mebibytes) or Gi (gibibytes). It is important to know that when a container exceeds its memory limit in Kubernetes, it is terminated (OOMKilled; Out Of Memory Killed). This is different from how CPU resources are handled, as exceeding the CPU limit results in throttling of the container but not termination. Therefore, setting adequate memory allocations is important for kubernetes!

To set appropriate memory constraints, both the average and peak memory usage (RSS) observed during task execution should be considered. As mententioned earlier `requests` is the minimum guaranteed resources and `limits` for the maximum allowed resources. By capturing metrics from a range of real task loads, the `requests` should be set to the encounterd average memory usage (plus a safety margin; 10-30%?) of the benchmarked tasks, and `limits` to the peak memory usage (plus a safety margin; 10-30%?). The `limits` value should be high enough to accommodate rare spikes, but not so high as to waste resources.

**Example:**

If a task's average memory usage is 600 MiB and the peak observed is 900 MiB, applying a 20% safety margin gives:

```
# Kustomize manifest example

resources:
    requests:
        memory: "720Mi"   # average (600 MiB) + 20% safety margin
    limits:
        memory: "1080Mi"  # peak (900 MiB) + 20% safety margin
```

This ensures the container is scheduled on a node with at least 700 MiB available, and will be killed if it ever exceeds 1 GiB. These values should be adjusted after monitoring real use-cases.
