#!/bin/bash
#
#SBATCH --job-name=R_1Y_2008
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000M
#SBATCH --time=50:00:00
#SBATCH --mail-user=vvalenzuela@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=def-allen
# stdout and stderr file paths/names
#SBATCH --output=/home/vicentev/scratch/vicentev/files_out_err/STDOUT
#SBATCH --error=/home/vicentev/scratch/vicentev/files_out_err/STDERR
#
source /home/vicentev/miniforge3/etc/profile.d/conda.sh
conda activate Parcels
#
WORKDIR=/home/vicentev/projects/def-allen/vicentev/analysis-vicente/OP_nibi
DRIVER=Model_Driver_nibi_from_yaml
CONFIG_R="$WORKDIR/config_files/config_file_OP_RESTART.yaml"
#
out_txt_path=/home/vicentev/scratch/vicentev/out
#
cd $WORKDIR
#

python -m "$DRIVER" "$CONFIG_R" &> "$out_txt_path/output_R_1Y.txt" &
#
wait
echo "All OP runs done"  