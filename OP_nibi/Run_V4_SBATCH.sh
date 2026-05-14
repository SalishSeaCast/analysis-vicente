#!/bin/bash
#
#SBATCH --job-name=NOV-DEC-V4
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=16000M
#SBATCH --time=100:00:00
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
DRIVER=Model_Driver_V2_final
CONFIG_P1=$WORKDIR/config_files_V4/config_file_V4_P1.yaml
CONFIG_P2=$WORKDIR/config_files_V4/config_file_V4_P2.yaml
CONFIG_P3=$WORKDIR/config_files_V4/config_file_V4_P3.yaml
CONFIG_P4=$WORKDIR/config_files_V4/config_file_V4_P4.yaml
CONFIG_P5=$WORKDIR/config_files_V4/config_file_V4_P5.yaml
CONFIG_P6=$WORKDIR/config_files_V4/config_file_V4_P6.yaml
CONFIG_P7=$WORKDIR/config_files_V4/config_file_V4_P7.yaml
CONFIG_P8=$WORKDIR/config_files_V4/config_file_V4_P8.yaml
CONFIG_P9=$WORKDIR/config_files_V4/config_file_V4_P9.yaml
CONFIG_P10=$WORKDIR/config_files_V4/config_file_V4_P10.yaml
CONFIG_P11=$WORKDIR/config_files_V4/config_file_V4_P11.yaml
CONFIG_P12=$WORKDIR/config_files_V4/config_file_V4_P12.yaml
#
out_txt_path=/home/vicentev/scratch/vicentev/out
#
cd $WORKDIR
#
#python -m "$DRIVER" "$CONFIG_P1" &> "$out_txt_path/output_V4_P1.txt" &
#python -m "$DRIVER" "$CONFIG_P2" &> "$out_txt_path/output_V4_P2.txt" &
#python -m "$DRIVER" "$CONFIG_P3" &> "$out_txt_path/output_V4_P3.txt" &
#python -m "$DRIVER" "$CONFIG_P4" &> "$out_txt_path/output_V4_P4.txt" &
#python -m "$DRIVER" "$CONFIG_P5" &> "$out_txt_path/output_V4_P5.txt" &
#python -m "$DRIVER" "$CONFIG_P6" &> "$out_txt_path/output_V4_P6.txt" &
#python -m "$DRIVER" "$CONFIG_P7" &> "$out_txt_path/output_V4_P7.txt" &
#python -m "$DRIVER" "$CONFIG_P8" &> "$out_txt_path/output_V4_P8.txt" &
#python -m "$DRIVER" "$CONFIG_P9" &> "$out_txt_path/output_V4_P9.txt" &
#python -m "$DRIVER" "$CONFIG_P10" &> "$out_txt_path/output_V4_P10.txt" &
python -m "$DRIVER" "$CONFIG_P11" &> "$out_txt_path/output_V4_P11.txt" &
python -m "$DRIVER" "$CONFIG_P12" &> "$out_txt_path/output_V4_P12.txt" &
#
wait
echo "Version 4 Done! :D"  