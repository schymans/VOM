#!/bin/sh
# coding: utf-8
#***********************************************************************
#        run_vom.sh
#        Compiles and runs the Vegetation Optimality Model (VOM).
#
#
#-----------------------------------------------------------------------
#        Authors: Remko Nijzink
#        Now at: LIST (Luxembourg Institute of Science and Technology)
#-----------------------------------------------------------------------
#
#  Copyright (C) 2022, Luxembourg Institute of Science and Technology, all rights reserved.
#
#    Code licensed under GPL version 2 or higher
#    SPDX-License-Identifier: GPL-2.0-or-later
#
#***********************************************************************


exe_dir=$1
inputdir=$2
input_weather=$3
input_soil=$4
nml_input=$5
outputdir=$6
restart_dir=$7


date

#compile code
C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/include/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/

#check if paths are empty, use default otherwise
if [ ! -d "$NC_INCLUDE" ]; then
   NC_INCLUDE=/usr/include
   echo $NC_INCLUDE
fi

if [ ! -d "$NC_LIB" ]; then
   NC_LIB=/usr/lib
   echo $NC_LIB
fi



make --directory $exe_dir FC=gfortran NC_INCLUDE=$NC_INCLUDE NC_LIB=$NC_LIB

#check if the outputdir exists and else make one
if [ ! -d "$outputdir" ]; then
   mkdir $outputdir
fi

if [ -f "$restart_dir/sce_lastloop.txt" ]; then
   cp $restart_dir/* $outputdir
fi

if [ ! -d "$inputdir" ]; then
mkdir $inputdir
cp $input_weather $inputdir
cp $input_soil $inputdir

fi

#run the model 
$exe_dir/model.x -i $inputdir -o $outputdir -n $nml_input

#clean again
make clean --directory $exe_dir

date
