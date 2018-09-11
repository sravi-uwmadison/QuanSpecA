# QuanSpecA
This github repository contains the code for replicating experiments of our CVPR 2018 paper: A Biresolution Spectral framework for Product Quantization. The paper can be accessed freely from: http://openaccess.thecvf.com/content_cvpr_2018/CameraReady/1103.pdf

In particular, our algorithm can be found in optimize_subspaces_blockcord.m.

The main file to call is mmf_test. This is setup to run on random data.
Your can do the following calls:
mmf_test
or
mmf_test(8,100,10)
or
mmf_test(8,100,10,16)


Need to install the yael-master codebase, since some functions are used from there
https://github.com/jackculpepper/yael
may need to recomile some mex files (used here) as well

Feel free to use your favorite package for Eigenvector and k-Means calculation.

