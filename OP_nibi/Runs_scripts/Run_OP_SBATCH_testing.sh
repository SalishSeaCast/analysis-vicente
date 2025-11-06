#!/bin/bash

#SBATCH --job-name=F-M-test-run
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=12000M
#SBATCH --time=50:00:00
#SBATCH --mail-user=vvalenzuela@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=def-allen
# stdout and stderr file paths/names
#SBATCH --output=/home/vicentev/scratch/vicentev/files_out_err/STDOUT
#SBATCH --error=/home/vicentev/scratch/vicentev/files_out_err/STDERR
#
# May need to change something with srun
source /home/vicentev/miniforge3/etc/profile.d/conda.sh
conda activate Parcels
#
WORKDIR="/home/vicentev/projects/def-allen/vicentev/analysis-vicente/OP_nibi"
DRIVER="Model_Driver_nibi_from_yaml"
#CONFIG_P1="$WORKDIR/config_files/config_file_OP_P1.yaml"
CONFIG_P2="$WORKDIR/config_files/config_file_OP_P2.yaml"
CONFIG_P3="$WORKDIR/config_files/config_file_OP_P3.yaml"
out_txt_path="/home/vicentev/scratch/vicentev/out"
#
cd $WORKDIR
#
python -m "$DRIVER" "$CONFIG_P2" &> "$out_txt_path/output_P2.txt" &
python -m "$DRIVER" "$CONFIG_P3" &> "$out_txt_path/output_P3.txt" &
#
wait
echo "All OP runs done"  