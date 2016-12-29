# Demo for the Rotated Compressive Sensing Method
## Xiu Yang

This is an introduction how to use the MATLAB demo code of the rotated compressive sensing method for the Hermite polynomial expansions. The algorithm is described in Yang X et al, *J Comput Phys*, 2016, doi:[10.1016/j.jcp.2015.11.038](http://dx.doi.org/10.1016/j.jcp.2015.11.038).

In this demo we use the first example in our paper, i.e.,

```
f(\xi_1,\xi_2,\cdots,\xi_d) = \sum_{i=1}^d \xi_i + 0.25 \left(\sum_{i=1}^d\xi_i\right)^2 + 0.025 \left(\sum_{i=1}^d\xi_i\right)^3,

d = 12
```
where $\xi_i$ are i.i.d. standard normal random variables.
We generate 160 samples from and then use $120$ of them to construct the Hermite polynomial expansion (surrogate model) and use the remaining 40 to examine the accuracy.
In our paper, we use the exact solution to examine the error.
In most practical problems, statistical tests (e.g., cross-validation, AIC, BIC, etc.) are used to evaluate the accuracy of the surrogate model.

The demo should be run with the following steps:

1. At the command line:
  1. Run ```git submodule init``` and ```git submodule update``` to make sure you have the latest version of SPGL1 from <https://github.com/mpf/spgl1> in the repository.
2. In MatLab:
  1. Run ```gen\_data.m``` to generate the inputs and outputs.  Here the inputs are samples of random vector $(\xi_1,\xi_2,\cdots,\xi_d)$ and outputs are the corresponding values of $f$.  The data is stored in ```data.mat```.
  2. Run ```gen\_kernel.m``` to generate the kernel matrices $K_{ij}$.  You need specify the dimension (`dim`) and polynomial order (```poly_order```).  The kernel matrices are stored in ```kernel_*.mat```.
4. Run ```demo.m```.  The relative root mean square error is computed to compare the accuracy of the standard $\ell_1$ result and the rotated $\ell_1$ results.  Here the inputs and outputs are loaded from ```data.mat```. The users should specify the dimension of the problem (```dim```), polynomial order (```poly_order```), the total number of samples (```num_sample```; i.e., size of the training set), and the number of rotation iterations (```num_iteration```).

Xiu reports the the results from MATLAB 2016a (Mac OS X?) are:

    >> demo
    finish reading data
    RMSE of standard l1

    ans =

    0.5676

    Start rotation
    RMSE after rotations

    ans =

    0.0331

However, Nathan gets a smaller value (~0.02) for the rotated RMSE with 2016a Matlab on Windows 7.
