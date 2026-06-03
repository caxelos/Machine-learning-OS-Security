# Preemptive SJF Scheduler Analysis



## Overview



This project explores the behavior of a **preemptive Shortest Job First (SJF)** CPU scheduling algorithm under different workload scenarios.



The objective is to analyze how process characteristics affect scheduling decisions and overall system performance. Various process profiles were tested, including interactive and non-interactive workloads, equal execution times, and increasing/decreasing execution-time distributions.



The scheduler's performance is evaluated using **Average Turnaround Time (ATT)** while also examining fairness and starvation-related issues.



---



## Features



* Preemptive SJF scheduling simulation

* Interactive and non-interactive process workloads

* Support for multiple scheduling profiles

* Comparison of scheduling behavior with and without waiting-time prioritization

* Performance evaluation through turnaround-time measurements

* Analysis of starvation and fairness trade-offs



---

## Experiment Profiles

The following workload configurations were used:

| Profile                 | Description                                                            |

| ----------------------- | ---------------------------------------------------------------------- |

| `onlycalc.conf`         | Different execution times, all interactive processes                   |

| `onlycalcboth.conf`     | Different execution times, mixed interactive/non-interactive processes |

| `allsametime.conf`      | Equal execution times, all non-interactive processes                   |

| `allsametimeinter.conf` | Equal execution times, all interactive processes                       |

| `simple1000.conf`       | Increasing execution times                                             |

| `simple.conf`           | Decreasing execution times                                             |



---



## Evaluation Metric



The primary metric used for evaluation is \*\*Average Turnaround Time (ATT)\*\*:



[

ATT = \\frac{\\sum TurnaroundTime\_i}{N}

]



where:



* `TurnaroundTime\_i` is the completion time of process \*i\*.

* `N` is the total number of processes.



Lower ATT values indicate better scheduling efficiency.



---


## Key Observations



### Fairness vs Performance



The experiments show that introducing a waiting-time component can improve fairness by reducing the likelihood of starvation for long-running processes.



### Workload Sensitivity



Scheduler performance varies significantly depending on the execution-time distribution and process type.



### Starvation in SJF



A known limitation of SJF is that newly arriving short jobs can continuously preempt longer jobs, causing starvation. Incorporating waiting time into scheduling decisions helps alleviate this issue.



---


---

## Performance Results

### Different Execution Times (in scheduler ticks)— All Interactive Processes

| Metric | With WaitingInRQ | Without WaitingInRQ |
|---------|---------|---------|
| Average Turnaround Time | 2898 | 3029 |

### Different Execution Times — Mixed Workload

| Metric | With WaitingInRQ | Without WaitingInRQ |
|---------|---------|---------|
| Average Turnaround Time | 1718 | 2081 |

### Equal Execution Times — All Non-Interactive

| Metric | With WaitingInRQ | Without WaitingInRQ |
|---------|---------|---------|
| Average Turnaround Time | 2308 | 2388 |

### Equal Execution Times — All Interactive

| Metric | With WaitingInRQ | Without WaitingInRQ |
|---------|---------|---------|
| Average Turnaround Time  | 3233 | 3361 |

### Increasing Execution Times

| Metric | With WaitingInRQ | Without WaitingInRQ |
|---------|---------|---------|
| Average Turnaround Time | 10500 | 9558 |

### Decreasing Execution Times

| Metric | With WaitingInRQ | Without WaitingInRQ |
|---------|---------|---------|
| Average Turnaround Time | 1168 | 1278 |

---

## Key Findings

### 1. Waiting-Time Awareness Improves Fairness

In most workload configurations, incorporating the **WaitingInRQ** factor reduced the average turnaround time and improved scheduling fairness.

### 2. Workload Characteristics Matter

The scheduler's effectiveness varied significantly depending on process execution patterns:

- Mixed workloads benefited greatly from WaitingInRQ.
- Equal-length workloads produced relatively stable results.
- Increasing execution-time profiles sometimes favored the standard SJF implementation.

### 3. Starvation Remains a Core Challenge

A known limitation of SJF scheduling is **starvation**.

When shorter jobs continuously arrive, they are repeatedly selected before longer-running processes. As a result, some processes may experience excessive waiting times or indefinite postponement.

### 4. Trade-Off Between Efficiency and Fairness

The experiments demonstrate the classic scheduling trade-off:

- Pure SJF tends to optimize turnaround time.
- Waiting-time-aware scheduling improves fairness.
- The optimal approach depends on workload characteristics and system goals.

---

## Technologies Used

- C
- Linux
- Operating System Scheduling Concepts
- Process Scheduling Simulation
- Performance Evaluation and Benchmarking

---

## Repository Structure

```text
P1_SJF_scheduler/
│
├── configs/
├── src/
├── results/
└── README.md
```

---

## Results Summary



The collected results demonstrate that no single scheduling strategy performs best in every scenario. Depending on the workload characteristics, prioritizing waiting processes may either improve fairness or increase average turnaround time.



This highlights the classic operating-system scheduling trade-off between:



* Maximum throughput

* Low turnaround time

* Fairness among processes



---



## Technologies



* Operating Systems concepts

* CPU Scheduling Algorithms

* Preemptive SJF Scheduling

* Performance Analysis

* Workload Simulation



---



\## Future Improvements



* Comparison with Round Robin scheduling

* Comparison with Priority Scheduling

* Visualization of scheduling timelines (Gantt Charts)

* Additional workload generators

* Statistical analysis of larger process sets



---



## License



This project is provided for educational and research purposes.



