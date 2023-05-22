#!/bin/bash
options="-pe smp 2 -binding linear:2 -l h_vmem=16G"
source /broad/software/scripts/useuse
use UGER
qsub $options enrichment.sh
