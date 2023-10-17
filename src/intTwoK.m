function out = intTwoK(u, w1, w2, V)
% intTwoK solves the double line integral of the kernel function K(x, xp)
%
% Double Integral of exp(-0.5*(u+t*w1-s*w2).'*(V\(u+t*w1-s*w2)))
% between 0 and 1 on both integrals
%
%
%
% Copyright (c) 2018-${BUILD_YEAR} Johannes Hendricks and others
% 
% Permission is hereby granted, free of charge, to any person obtaining
% a copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to
% permit persons to whom the Software is furnished to do so, subject to
% the following conditions:

% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
% LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
% OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
% WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%
% DO NOT EDIT! Automatically generated: 
% template file: ${MATLAB_TEMPLATE_FILE} 
%    build date: ${BUILD_DATE} 
%    build time: ${BUILD_TIME} 

