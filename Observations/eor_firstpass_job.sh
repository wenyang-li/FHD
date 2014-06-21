#! /bin/bash
#############################################################################
#Runs one observation at a time in grid engine.  Second level program for 
#running firstpass on the MIT machines. First level program is 
#batch_firstpass.sh
#############################################################################

#$ -V
#$ -N firstpass
#$ -S /bin/bash

# obs_id, outdir, version and nslots expected to be passed from qsub call
# in batch_firstpass.sh

echo JOBID ${JOB_ID}
echo obsid ${obs_id}

/usr/local/bin/idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $nslots -e eor_firstpass_versions -args $obs_id $outdir $version 

if [ $? -eq 0 ]
then
    echo "Finished"
    exit 0
else
    echo "Job Failed"
    exit 1
fi

