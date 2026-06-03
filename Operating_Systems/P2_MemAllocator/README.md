

## Overview

This project modifies the Linux SLOB (Simple List Of Blocks) memory allocator to reduce external memory fragmentation. The default allocation strategy (First-Fit / Next-Fit) was replaced with a Best-Fit approach for both page and block selection.

## Features

* Implementation of the Best-Fit allocation algorithm in the SLOB allocator.
* Runtime logging of allocation decisions and candidate memory blocks.
* Custom system calls for memory statistics:

  * Total allocated memory.
  * Total free memory.
* Comparison framework between the original and modified allocator.

## Objectives

The goal of the project is to evaluate whether the Best-Fit strategy can reduce external fragmentation compared to the original SLOB implementation.

## Project Structure

* Modified Linux kernel source code (SLOB allocator).
* User-space application for collecting allocator statistics.
* Experimental results and observations.
* Patch files containing kernel modifications.

## Build & Usage

1. Apply the provided patch to a Linux 2.6.37.2 kernel source tree.
2. Enable the SLOB allocator in the kernel configuration.
3. Compile and install the modified kernel.
4. Build and run the user-space test application.
5. Use the provided system calls to collect and compare memory allocation statistics.


