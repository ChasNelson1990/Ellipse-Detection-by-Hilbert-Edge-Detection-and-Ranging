## Ellipse Detection by Hilbert-Edge Detection and Ranging (HEDAR)
This repository contains our implementation for the paper "Ellipse Detection by Hilbert-Edge Detection and Ranging (HEDAR)"

- by **Carl J. Nelson, Philip T. G. Jackson and Boguslaw Obara**

## Overview
This repository provides a MATLAB implementation for detecting ellipses and ellipse-like objectsin 2D greyscale images.

This project contains three MATLAB functions:
1. Main HEDAR code, see Algorithm 2 in the paper (hedar.m)
2. A modified version of MATLAB's imerode (moderode.m)
3. A helper function that creates binary images of ellipses (ellipse2.m)

For general usage:
- Open MATLAB, right click project directory and add to path.
- From the command prompt, call `data=hedar(im);` to run on image `im` with default options.
- Alternatively, define the parameters, as discussed in the paper and the code preamble.

## License
This code is licensed under the GNU General Public License Version 3.
- For alternative licenses, please contact *carl.nelson@durham.ac.uk*.
