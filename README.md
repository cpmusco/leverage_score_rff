This directory contains all of the code and data necessary to reproduce all experiments and plots from the NeurIPS 2020 paper [Fourier Sparse Leverage Scores and Approximate Kernel Learning](https://arxiv.org/abs/2006.07340). All code was run in MATLAB version 2020a -- it should be compatible with older versions of MATLAB, except for some of the plotting functions used (like exportgraphics).

* empirical_upper_bounds.m will reproduce Figure 1. 
* experiments.m will reproduce all Figures in Section 4. 
* gaussianKernelMRFF.m contains our algorithm for the Gaussian kernel.
* gaussianKernelRFF.m contains the original random Fourier features algorithm.
* cauchyKernelMRFF.m contains our algorithm for the Cauchy kernel.
* cauchyKernelRFF.m contains the original random Fourier features algorithm.

