#!/bin/bash

# Makes two movie for a given PyECLOUD simulation for two subsequent bunch passages.
# The simulation folder must be a subfolder of this directory.
# This script must be run from its directory.

# Last tested with version 6.3.1 in August 2017

set -e

# PyECLOUD must be in pythonpath, edit this line accordingly
export PYTHONPATH=$PYTHONPATH

# bunch passage at which to make the pngs

if [[ $# == 2 ]] ; then
	folder=$1
	pass_to_save=$2
else
	echo "Incorrect usage of $0.
	Usage: $0 folder bunch_passage
	where folder is a subdirectory of where this script is located.
	and bunch_passage is the bunch passage at which to make the movie"
	exit 1
fi

cd $folder
echo "Using folder $folder"
sleep 2s

echo "

#----created automatically---- by $0
save_simulation_state_time_file = [${pass_to_save}.*25e-9]
stopfile = 'simulation_state_0.pkl'
logfile_path = './logfile'
progress_path = './progress'
" >> simulation_parameters.input

python2 ../000_run_simulation.py || true

echo "stopfile = './stop'" >> simulation_parameters.input

python2 ../002_reload_state_run_and_save_movie.py

cd ..
python2 ./004_VIDEO_high_res.py --folder_sim $folder 

