#!/bin/bash
# WORKS WELL!
# request an interactive allocation and run OP script once granted
#
WORKDIR="/home/vicentev/projects/def-allen/vicentev/analysis-vicente/OP_nibi"
DRIVER="Model_Driver_nibi_from_yaml"
CONFIG="$WORKDIR/config_files/config_file_OP_P1.yaml"

salloc --time=6:00:00 --mem-per-cpu=8G --ntasks=1 --account=rrg-allen bash -c "
    cd \"$WORKDIR\" && \
    source /home/vicentev/miniforge3/etc/profile.d/conda.sh && \
    conda activate Parcels && \
    python -m $DRIVER $CONFIG
"   