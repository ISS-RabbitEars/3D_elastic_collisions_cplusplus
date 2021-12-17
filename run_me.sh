#!/bin/bash

mkdir pos
mkdir pos/cam
g++ collisions.cpp -o col
./col
rm col
g++ cam_path.cpp -o cam
./cam
rm cam
cp render.sh pos/
cp pos_template.pov pos/
cd pos
./render.sh
#ffmpeg -framerate 60 -i frame_%d.png 2dcol_pbc.gif
ffmpeg -framerate 60 -i frame_%d.png 2dcol_pbc.mp4
rmdir cam
