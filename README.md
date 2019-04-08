## Ellipse Detection by Hilbert-Edge Detection and Ranging (HEDAR)
This repository contains a MATLAB implementation of our HEDAR method as used in the paper "Combining Mathematical Morphology and the Hilbert Transform for Fully Automatic Nuclei Detection in Fluorescence Microscopy" (submitted to LNCS)

- by **Carl J. Nelson, Philip T. G. Jackson and Boguslaw Obara**

## Overview
This repository provides a MATLAB implementation for detecting ellipses and ellipse-like objects in 2D greyscale images.

This project contains three MATLAB functions:
1. Main HEDAR code, see Algorithm 2 in the paper (hedar.m)
3. A helper function that creates binary images of ellipses (ellipse2.m)

For general usage:
- Open MATLAB, right click project directory and add to path.
- From the command prompt, call `data=hedar(im);` to run on image `im` with default options.
- Alternatively, define the parameters, as discussed in the paper (submitted) and the code preamble.

## License
This code is licensed under the GNU General Public License Version 3.
- For alternative licenses, please contact *carl.nelson@durham.ac.uk*.
