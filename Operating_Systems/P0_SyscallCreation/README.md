# Linux Kernel Development Project

## Overview

This project explores Linux kernel internals through two practical modifications:

1. Implementing a custom system call inside the Linux kernel.
2. Modifying and deploying a custom I/O scheduler kernel module.

---

## Part 1: Custom System Call — `find_roots()`

### Objective

I implemented a new Linux system call named `find_roots()` that traverses the ancestry of the calling process and prints information about each process in the hierarchy.

### Function Prototype

```c
int find_roots(void);
```

### Functionality

When invoked, the system call:

* Retrieves the currently executing process.
* Prints the process identifier (PID) and executable name.
* Iteratively follows the parent process chain.
* Prints the PID and name of each ancestor process.
* Stops when the root process (`init`, PID 1) is reached.

Example output:

```text
find_roots system call called by process 1869

id: 1869, name: find_roots_lib
id: 1830, name: bash
id: 1824, name: konsole
id: 1, name: init
```

### Kernel Concepts Used

* System call registration
* `task_struct`
* Process Control Blocks (PCB)
* Process hierarchy traversal
* Kernel logging (`printk`)
* Kernel recompilation and installation

## Part 2: Custom I/O Scheduler Module

### Objective

To understand Linux block I/O scheduling, I modified the source code of the NOOP scheduler and built it as a custom kernel module.

### Modifications

The scheduler was renamed and isolated from the original implementation by:

* Renaming scheduler functions.
* Updating scheduler registration structures.
* Customizing module metadata.
* Registering the scheduler under a new name.

Additionally, diagnostic logging was added to the request dispatch routine:

```c
printk("In custom_noop_dispatch function\n");
```

This message is emitted whenever the scheduler dispatches an I/O request, allowing verification that the modified scheduler is actively handling disk operations.

### Testing Procedure

1. Build the module using the Linux kernel build system.
2. Load the module with `insmod`.
3. Register the scheduler for a block device.
4. Generate disk I/O activity.
5. Verify scheduler execution through kernel log messages (`dmesg`).
6. Restore the default scheduler and unload the module.

---

## Build Process

### System Call

After modifying the kernel source:

```bash
make
make modules
make install
```

Reboot into the newly compiled kernel.

### Kernel Module

```bash
make
sudo insmod custom-iosched.ko
```

Verify loading:

```bash
lsmod
dmesg | tail
```

Remove the module:

```bash
sudo rmmod custom-iosched
```

---

## Technologies

* C
* Linux Kernel (2.6.x series)
* GNU Make
* GCC
* Ubuntu Linux
* Linux Kernel Build System (kbuild)

---
