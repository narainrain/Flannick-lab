#!/bin/bash
out_folder="/humgen/diabetes2/users/ellamas/out/"  # make sure the folder is empty before running the script
options="-pe smp 1 -binding linear:1 -l h_vmem=1G"
source /broad/software/scripts/useuse
use Python-3.6
source /home/unix/ellamas/test/bin/activate
python generate_sp.py
use UGER
qsub -t 1-250 $options -o $out_folder sp1.sh
qsub -t 1-250 $options -o $out_folder sp2.sh
qsub -t 1-210 $options -o $out_folder sp3.sh
qsub -t 1-200 $options -o $out_folder sp4.sh
qsub -t 1-190 $options -o $out_folder sp5.sh
qsub -t 1-180 $options -o $out_folder sp6.sh
qsub -t 1-160 $options -o $out_folder sp7.sh
qsub -t 1-150 $options -o $out_folder sp8.sh
qsub -t 1-140 $options -o $out_folder sp9.sh
qsub -t 1-140 $options -o $out_folder sp10.sh
qsub -t 1-140 $options -o $out_folder sp11.sh
qsub -t 1-140 $options -o $out_folder sp12.sh
qsub -t 1-120 $options -o $out_folder sp13.sh
qsub -t 1-110 $options -o $out_folder sp14.sh
qsub -t 1-110 $options -o $out_folder sp15.sh
qsub -t 1-100 $options -o $out_folder sp16.sh
qsub -t 1-90 $options -o $out_folder sp17.sh
qsub -t 1-80 $options -o $out_folder sp18.sh
qsub -t 1-70 $options -o $out_folder sp19.sh
qsub -t 1-70 $options -o $out_folder sp20.sh
qsub -t 1-50 $options -o $out_folder sp21.sh
qsub -t 1-50 $options -o $out_folder sp22.sh
