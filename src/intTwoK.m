function out = intTwoK(uij, wi, wj, V)%#ok<STOUT,INUSD>
% intTwoK solves the double line integral of the kernel squared exponential function
%
%                             /1 /1   1 
% K_{ij} = ||w_{i}|| ||w_{j}|||  |  - - exp((r_{i}(t) - r_{j}(s)).' * V * (r_{i}(t) - r_{j}(s)) ds dt,
%                             /0 /0   2
%
% where r_{i}(t) = p_{i} + t * w_{i} and r_{j}(s) = p_{j} + s * w_{j}. The
% origins are combined such that u_{ij} = p_{i} - p_{j}.
%
% out = intTwoK(uij, wi, wj, V)
% out is a vector containing 4 elements
% out(1): results
% out(2): absolute error estimated by the gsl integration routine
% out(3): number of function evaluations required by the gsl integration 
%         routine, max = 81, if out(3) == -1, then L1 = L2 = 0, and the 
%         function evaluated a constant.
% out(4): gsl integration routine error code, if code == 0 then no issues. 
%         Otherwise see gsl see referenced paper for description of the inputs
%
% Associated publication:
% Hendriks, J. N., Jidling, C., Wills, A., & Schön, T. B. (2018). 
% Evaluating the squared-exponential covariance function in Gaussian processes 
% with integral observations: Technical report, arXiv:1812.07319.
%
% Copyright (c) 2018-${BUILD_YEAR} Johannes Hendriks and others
%
% Permission is hereby granted, free of charge, to any person obtaining
% a copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to
% permit persons to whom the Software is furnished to do so, subject to
% the following conditions:
%
% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
% LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
% OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
% WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% DO NOT EDIT! Automatically generated:
% template file: ${MATLAB_TEMPLATE_FILE}
%    build date: ${BUILD_DATE}
%    build time: ${BUILD_TIME}

