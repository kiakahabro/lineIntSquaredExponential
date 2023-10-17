function out = intTwoSimps(u, w1, w2, V)
% intTwoSimps solves the double line integral of the kernel function K(x, xp)
%
% Double Integral of exp(-0.5*(u+t*w1-s*w2).'*(V\(u+t*w1-s*w2)))
% between 0 and 1 on both integrals
