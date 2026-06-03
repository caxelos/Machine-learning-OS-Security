# Project 4 – Linux C-LOOK I/O Scheduler

## Overview

This project implements a custom **C-LOOK disk I/O scheduler** for the Linux kernel.

The scheduler manages disk requests by:

* Keeping requests sorted by sector number.
* Servicing requests in one direction.
* Jumping back to the lowest pending request after reaching the highest one.

## Objectives

* Create a new Linux I/O scheduler (`clook-iosched`).
* Integrate it into the Linux kernel build system.
* Implement the C-LOOK scheduling algorithm.
* Compile and load the scheduler as a kernel module.
* Test the scheduler using a dedicated virtual disk.
* Verify correct behavior through logging and experiments.

## Key Features

* Custom request queue management.
* Ordered request dispatch using C-LOOK.
* Kernel logging for request insertion and dispatch:

  * `[CLOOK] add <R|W> <sector>`
  * `[CLOOK] dsp <R|W> <sector>`

## Evaluation

The scheduler is validated through controlled disk I/O experiments that demonstrate requests are serviced in the expected C-LOOK order.

## Deliverables

* Linux kernel patch (`patch_4`)
* Scheduler source code
* Experimental evaluation program(s)
* Sample output logs
* Analysis report
* `README.txt` with build instructions and team information
