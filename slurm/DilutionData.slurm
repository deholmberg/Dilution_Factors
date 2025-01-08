#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=800
#SBATCH --job-name=dilution_factor
#SBATCH --cpus-per-task=1
#SBATCH --output=/farm_out/holmberg/dilution_factor-%j.out
#SBATCH --error=/farm_out/holmberg/dilution_factor-%j.err
${4}/Make_DF_Files.sh ${1} ${2} ${3} ${4}
