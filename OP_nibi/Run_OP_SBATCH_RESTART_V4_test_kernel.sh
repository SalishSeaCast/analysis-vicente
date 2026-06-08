#!/bin/bash
#
#SBATCH --job-name=test_NOV-DEC-08
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
CONFIG_P1=$WORKDIR/test_kernels/config_file_V4_RP1.yaml
CONFIG_P2=$WORKDIR/test_kernels/config_file_V4_RP2.yaml
CONFIG_P3=$WORKDIR/test_kernels/config_file_V4_RP3.yaml
CONFIG_P4=$WORKDIR/test_kernels/config_file_V4_RP4.yaml
CONFIG_P5=$WORKDIR/test_kernels/config_file_V4_RP5.yaml
CONFIG_P6=$WORKDIR/test_kernels/config_file_V4_RP6.yaml
CONFIG_P7=$WORKDIR/test_kernels/config_file_V4_RP7.yaml
CONFIG_P8=$WORKDIR/test_kernels/config_file_V4_RP8.yaml
CONFIG_P9=$WORKDIR/test_kernels/config_file_V4_RP9.yaml
CONFIG_P10=$WORKDIR/test_kernels/config_file_V4_RP10.yaml
CONFIG_P11=$WORKDIR/test_kernels/config_file_V4_RP11.yaml
CONFIG_P12=$WORKDIR/test_kernels/config_file_V4_RP12.yaml
#
out_txt_path=/home/vicentev/scratch/vicentev/out
#
cd $WORKDIR
#
#python -m "$DRIVER" "$CONFIG_P1" &> "$out_txt_path/output_RV4_P1.txt" &
#python -m "$DRIVER" "$CONFIG_P2" &> "$out_txt_path/output_RV4_P2.txt" &
#python -m "$DRIVER" "$CONFIG_P3" &> "$out_txt_path/output_RV4_P3.txt" &
#python -m "$DRIVER" "$CONFIG_P4" &> "$out_txt_path/output_RV4_P4.txt" &
#python -m "$DRIVER" "$CONFIG_P5" &> "$out_txt_path/output_RV4_P5.txt" &
#python -m "$DRIVER" "$CONFIG_P6" &> "$out_txt_path/output_RV4_P6.txt" &
#python -m "$DRIVER" "$CONFIG_P7" &> "$out_txt_path/output_RV4_P7.txt" &
#python -m "$DRIVER" "$CONFIG_P8" &> "$out_txt_path/output_RV4_P8.txt" &
#python -m "$DRIVER" "$CONFIG_P9" &> "$out_txt_path/output_RV4_P9.txt" &
#python -m "$DRIVER" "$CONFIG_P10" &> "$out_txt_path/output_RV4_P10.txt" &
python -m "$DRIVER" "$CONFIG_P11" &> "$out_txt_path/output_RV4_P11.txt" &
python -m "$DRIVER" "$CONFIG_P12" &> "$out_txt_path/output_RV4_P12.txt" &
#
wait
echo "V4 2008 Kernel Test runs are done :D"  