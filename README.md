# Plates
Matlab Plate Recognition

This project aims to read any car plate from any image. Some Peter Kovesi's functions
are provided but you can find them here: http://www.peterkovesi.com/matlabfns/#spatial.

The next step is to shape the ROI rectification (basically a homography) based on 4 good Harris
corners and good parallel lines detected by a Hough transfom.
