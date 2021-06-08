#!/bin/bash

export PGI_ACC_NOTIFY=

jsrun -r1 -a1 -c1 -g1 ../../src/pre_process_code/pre_process
time jsrun -r1 -a1 -c1 -g1 ../../src/simulation_code/simulation
