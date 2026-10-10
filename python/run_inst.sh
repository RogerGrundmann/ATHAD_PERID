#!/bin/bash
# The instruments of ATHAD_COND (ATM_SATADJ_DIAG, ATM_MC_CAP_DIAG, ATM_MC_DIAG) and its ATM_MC_* knobs, ported 2026-10-10.
# Byte checks, 2 iterations, 1 thread:  tiedtke with every knob unset = the binary before the convection switch (cli/perid_pre);
#                                       adjust (default) = the binary before this port (cli/perid_long).
# First readings, 20 iterations from scratch, 8 threads: tiedtke with the three instruments on; adjust with ATM_SATADJ_DIAG.
cd /home/roger/SynologyDrive/Cloudstation_Notebook/ATHAD_PERID/python
for d in run_in_vpre run_in_vt run_in_vlong run_in_vadj; do mkdir -p $d; cp run_ca_vpre/config.xml $d/; done
for d in run_in_tied run_in_adj; do mkdir -p $d; cp run_ca_tied/config.xml $d/; done
( cd run_in_vpre  && OMP_NUM_THREADS=1 ../../cli/perid_pre config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_in_vt    && OMP_NUM_THREADS=1 ATM_CONV_SCHEME=tiedtke ../../cli/perid config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_in_vlong && OMP_NUM_THREADS=1 ../../cli/perid_long config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_in_vadj  && OMP_NUM_THREADS=1 ../../cli/perid config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_in_tied  && env OMP_NUM_THREADS=8 ATM_CONV_SCHEME=tiedtke ATM_SATADJ_DIAG=1 ATM_MC_CAP_DIAG=1 ATM_MC_DIAG=1 ../../cli/perid config.xml > run.log 2>&1; echo $? > exit ) &
( cd run_in_adj   && env OMP_NUM_THREADS=8 ATM_SATADJ_DIAG=1 ../../cli/perid config.xml > run.log 2>&1; echo $? > exit ) &
wait
cmpdir(){ n=0; s=0; for f in $(cd $1 && find . -type f); do n=$((n+1)); cmp -s $1/$f $2/$f && s=$((s+1)); done; echo "ATHAD_PERID: $1 = $2 : $s/$n files identical"; }
cmpdir run_in_vpre/out run_in_vt/out
cmpdir run_in_vlong/out run_in_vadj/out
echo "ATHAD_PERID exits: $(cat run_in_*/exit | tr '\n' ' ')"
