#!/bin/bash

# Since there isn't QADB for the CD2 targets yet, we need to turn it off for the CD2 targets
#echo ${1} ${2} ${3} ${4} | clas12root -l -q ../Get_DF_By_Sector.C

echo ${1} ${2} ${3} ${4} | clas12root -l -q ../Get_DF_By_Sector.C

