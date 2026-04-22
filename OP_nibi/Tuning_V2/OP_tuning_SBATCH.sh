#!/bin/bash
#
#SBATCH --job-name=TestS-T1-T2
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=16000M
#SBATCH --time=4:00:00
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

CONFIG_T1=$WORKDIR/yaml_MC_0_2_Tau_0_001_Ads_0_05_Vel_N_TEST_STUCK.yaml
CONFIG_T2=$WORKDIR/config_ADS_0_05_TAU_0_0075_MC_0_2_TEST_STUCK.yaml

#CONFIG_P1=$WORKDIR/config_tuning_V2/yaml_MC_0_2_Tau_0_0025_Ads_0_05_Vel_N.yaml
#CONFIG_P2=$WORKDIR/config_tuning_V2/yaml_MC_0_2_Tau_0_005_Ads_0_05_Vel_N.yaml
#CONFIG_P3=$WORKDIR/config_tuning_V2/yaml_MC_0_2_Tau_0_001_Ads_0_1_Vel_N.yaml
#CONFIG_P4=$WORKDIR/config_tuning_V2/yaml_MC_0_2_Tau_0_001_Ads_0_01_Vel_N.yaml
#CONFIG_P5=$WORKDIR/config_tuning_V2/yaml_MC_0_2_Tau_0_001_Ads_0_01_Vel_Hx1_2.yaml
#CONFIG_P6=$WORKDIR/config_tuning_V2/yaml_MC_0_2_Tau_0_001_Ads_0_1_Vel_Hx1_2.yaml
#CONFIG_P7=$WORKDIR/config_tuning_V2/yaml_MC_0_2_Tau_0_001_Ads_0_05_Vel_Hx1_5.yaml
#CONFIG_P8=$WORKDIR/config_tuning_V2/yaml_MC_0_2_Tau_0_0025_Ads_0_05_Vel_Hx1_5.yaml

#
out_txt_path=/home/vicentev/scratch/vicentev/out_tuning
#
cd $WORKDIR
#

python -m "$DRIVER" "$CONFIG_T1" &> "$out_txt_path/output_Test_1.txt" &
python -m "$DRIVER" "$CONFIG_T2" &> "$out_txt_path/output_Test_2.txt" &
#python -m "$DRIVER" "$CONFIG_P3" &> "$out_txt_path/output_TVR2-P3.txt" &
#python -m "$DRIVER" "$CONFIG_P4" &> "$out_txt_path/output_TVR2-P4.txt" &
#python -m "$DRIVER" "$CONFIG_P5" &> "$out_txt_path/output_TVR2-P5.txt" &
#python -m "$DRIVER" "$CONFIG_P6" &> "$out_txt_path/output_TVR2-P6.txt" &
#python -m "$DRIVER" "$CONFIG_P7" &> "$out_txt_path/output_TVR2-P7.txt" &
#python -m "$DRIVER" "$CONFIG_P8" &> "$out_txt_path/output_TVR2-P8.txt" &

#
wait
echo "Test tuning done! :D"  