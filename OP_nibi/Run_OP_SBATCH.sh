#!/bin/bash
#
#SBATCH --job-name=JUL-2007
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=14000M
#SBATCH --time=60:00:00
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
#CONFIG_P1=$WORKDIR/config_files/config_file_OP_P1.yaml
#CONFIG_P2="$WORKDIR/config_files/config_file_OP_P2.yaml"
#CONFIG_P3="$WORKDIR/config_files/config_file_OP_P3.yaml"
#CONFIG_P4="$WORKDIR/config_files/config_file_OP_P4.yaml"
#CONFIG_P5="$WORKDIR/config_files/config_file_OP_P5.yaml"
#CONFIG_P6="$WORKDIR/config_files/config_file_OP_P6.yaml"
CONFIG_P7="$WORKDIR/config_files/config_file_OP_P7.yaml"
#CONFIG_P8="$WORKDIR/config_files/config_file_OP_P8.yaml"
#CONFIG_P9="$WORKDIR/config_files/config_file_OP_P9.yaml"
#CONFIG_P10="$WORKDIR/config_files/config_file_OP_P10.yaml"
#CONFIG_P11="$WORKDIR/config_files/config_file_OP_P11.yaml"
#CONFIG_P12="$WORKDIR/config_files/config_file_OP_P12.yaml"
#
out_txt_path=/home/vicentev/scratch/vicentev/out
#
cd $WORKDIR
#
#python -m "$DRIVER" "$CONFIG_P1" &> "$out_txt_path/output_P1.txt" &
#python -m "$DRIVER" "$CONFIG_P2" &> "$out_txt_path/output_P2.txt" &
#python -m "$DRIVER" "$CONFIG_P3" &> "$out_txt_path/output_P3.txt" &
#python -m "$DRIVER" "$CONFIG_P4" &> "$out_txt_path/output_P4.txt" &
#python -m "$DRIVER" "$CONFIG_P5" &> "$out_txt_path/output_P5.txt" &
#python -m "$DRIVER" "$CONFIG_P6" &> "$out_txt_path/output_P6.txt" &
python -m "$DRIVER" "$CONFIG_P7" &> "$out_txt_path/output_P7.txt" &
#python -m "$DRIVER" "$CONFIG_P8" &> "$out_txt_path/output_P8.txt" &
#python -m "$DRIVER" "$CONFIG_P9" &> "$out_txt_path/output_P9.txt" &
#python -m "$DRIVER" "$CONFIG_P10" &> "$out_txt_path/output_P10.txt" &
#python -m "$DRIVER" "$CONFIG_P11" &> "$out_txt_path/output_P11.txt" &
#python -m "$DRIVER" "$CONFIG_P12" &> "$out_txt_path/output_P12.txt" &
#
wait
echo "All OP runs done"  