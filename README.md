# Dilution Factors

This code is used to calculate quantities called "dilution factors" and "packing fractions" for our experiment in Run Group C of the CLAS12 Collaboration. Here's a quick overview of how to use this code.

NOTES: This code is designed to be run on the Jefferson Lab ifarm server!! If you want to generate the necessary text files for analysis, you'll have to clone this repository to your home or work directory on the farm. However, if you already have the necesasry text files, you should be able to run the "Dilution_Factor_Analysis.C" ROOT macro using ROOT 6.30 or later.

THE DF AND PF CALCULATION FOR ND3 TARGETS DOES NOT WORK YET!! I still need to properly implement the correct derivation of the DF and PF for ND3 targets, so don't try to use it for these targets (yet!).

## To Run the Analysis Code
Log on to the farm and clone this repository to wherever you like:

> git clone https://github.com/deholmberg/Dilution_Factors.git

From there, go into the slurm/ directory and make the following changes:

1) Go into "DilutionData.slurm" and change the slurm "farm_out/holmberg/" to "farm_out/<your_username>/" . This ensures that any errors/output messages go into your "farm_out/" directory.
2) In the "SubmitRuns.sh" script, there are two variables "VERSION" and "PERIOD" which are initialized to "pass1" and "summer22" respectively. If you want to analyze data from other run periods (like from "fall22"), make the adjustments in this file.
3) From the slurm directory, type "./SubmitRuns.sh". This will send all the relevant jobs to the farm. The "Make_DF_Files.sh" script will be called in each of these jobs to run the "Get_DF_By_Sector.C" ROOT macro on each of the relevant runs passed in from "SubmitRuns.sh". This could take up to 20-30 minutes for all the runs to finish. In the "Text_Files/" directory, you should see files of the form "<TARGET\_TYPE>\_<RUN\_NUMBER>\_DF\_Data.txt" and "<TARGET\_TYPE>\_<RUN\_NUMBER>\_DF\_Data\_S\*.txt". These files contain the number of counts of scattered electrons in each of the kinematic bins of Bjorken **x** and **Q^2** as well as the total FC charge for each of the runs. The "S<X>" files contain the same information for each of the individual sectors (1-6).
4) This is something clunky that I still need to iron out, but before running the analysis script, be sure to move all the files of run 16194 from "ET_16194..." to "F_16194..." so the foil-only data can be properly read in.
5) Run the "Dilution_Factor_Analysis.C" script in ROOT ("root Dilution\_Factor\_Analysis.C"). This will generate and plot dilution factors and packing fractions for all **x**, **Q^2** bins for a sample run range in NH3. You can edit this script to make plots for other run "epochs".

I know this code kinda sucks, but it's still a work in progress. I hope this repository will be helpful for anyone trying to figure out how I did things!

## Acknowledgements
I want to thank Sebastian Khun and Darren Upton for their help and feedback in putting this code together. I also want to thank Gregory Matousek, since I used an edited verison of his code in the scripts "SubmitRuns.sh" and "DilutionData.slurm". Without his code as a template, I wouldn't have been able to figure out how to use slurm to send jobs to the farm.

