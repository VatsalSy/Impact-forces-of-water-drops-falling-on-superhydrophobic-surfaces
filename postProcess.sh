#!/bin/bash

ci="1000"
We="20.0"
Ohd="0.00333"
Rhor="0.001"  
Ohs="1e-5"
Bo="0.25"
Lo="4.0"
nR="512"
tsnap="0.01"

Video="0"
a="1"

if [ $Video = $a ]
then

ffmpeg -framerate 30 -pattern_type glob -i 'VideoD2_v1/*.png' -vf scale=850:880 -c:v libx264 -r 30 -pix_fmt yuv420p $ci.mp4 -y &
python getEpsNForce.py $ci $Ohd $We &
python getEnergyScript.py $ci $Rhor $Ohd $Ohs $Bo $We &
wait

else
python VideoFullDomain.py $Lo $nR $We $Ohd $tsnap
fi

