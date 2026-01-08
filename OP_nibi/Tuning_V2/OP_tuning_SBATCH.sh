#!/bin/bash
#
#SBATCH --job-name=Adsorption-Test_n3-2009
#SBATCH --ntasks=3
#SBATCH --mem-per-cpu=16000M
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
WORKDIR=/home/vicentev/projects/def-allen/vicentev/analysis-vicente/OP_nibi/Tuning_V2
DRIVER=Model_Driver_V2
CONFIG_P1=$WORKDIR/config_tuning/config_adsorption_1.yaml
CONFIG_P2=$WORKDIR/config_tuning/config_adsorption_2.yaml
CONFIG_P3=$WORKDIR/config_tuning/config_adsorption_3.yaml
#
out_txt_path=/home/vicentev/scratch/vicentev/out_tuning
#
cd $WORKDIR
#
python -m "$DRIVER" "$CONFIG_P1" &> "$out_txt_path/output_AD_P1.txt" &
python -m "$DRIVER" "$CONFIG_P2" &> "$out_txt_path/output_AD_P2.txt" &
python -m "$DRIVER" "$CONFIG_P3" &> "$out_txt_path/output_AD_P3.txt" &

#
wait
echo "All OP runs done"  