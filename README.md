# SNIC superpixels

This repository offers the code for the Simple Non-Iterative Clustering (SNIC) algorithm with a Python interface. One can produce superpixels for grey, color, as well as images with any other number of channels. A demo file is provided, which should be easy to use. When using for the first time, the C files need to be compiled using:
```
python compile_snic_lib.py
```

The demo can then be run using:
```
python SNICDemo.py
```
The output of the demo should be as follows:

<p float="center">
  <img src="https://github.com/achanta/SNIC/blob/master/snic_python_interface/bee.png" width="400" />
  <img src="https://github.com/achanta/SNIC/blob/master/snic_python_interface/bee_snic.png" width="400" /> 
</p>

If you use the code, please cite the following publication:

"Superpixels and Polygons Using Simple Non-Iterative Clustering", R. Achanta and S. Süsstrunk, CVPR 2017.

```
@InProceedings{Achanta_2017_CVPR,
author = {Achanta, Radhakrishna and Susstrunk, Sabine},
title = {Superpixels and Polygons Using Simple Non-Iterative Clustering},
booktitle = {The IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
month = {July},
year = {2017}
}
```
