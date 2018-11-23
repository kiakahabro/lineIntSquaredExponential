# lineIntSquaredExponential
For evaluating double line integrals of the squared exponential covariance function

For details see 

Two example matlab scripts are included
1. comparison_set_2.m: Runs set two from the paper and plots the histogram
2. comparison_all_sets.m: Runs all sets from the paper and plots the output. Note, the log10 of the error is being plotted which for set 4 is -inf and so displays as a blank point on the plot

Source code for the mex function "intTwoK.c" is included and has been mexxed and tested on OSX 10.12.6 and UNIX using gsl 2.5. These files are
1. intTwoK.mexmaci64
2. intTwoK.mexa64

On OSX 10.12.6 this can be mexxed using
mex -O intTwoK.c -I<path_to_gsl> -lmwblas -lgsl

On Unix <version number> this can be mexxed using
mex -O intTwoK.c -lmwblas -lgsl -lgslcblas -lm

The mex function should be called as

OUT = intTwoK(u_ij,w_i,w_j,V)
    OUT is a vector containing 4 elemen ts
    OUT(1): results
    OUT(2): absolute error estimated by the gsl integration routine
    OUT(3): number of function evaluations required by the gsl integration routine, max = 81, if OUT(3) == -1, then L1 = L2 = 0, and the function evaluated a constant.
    OUT(4): gsl integration routine error code, if code == 0 then no issues. Otherwise see gsl
    see referenced paper for description of the inputs
