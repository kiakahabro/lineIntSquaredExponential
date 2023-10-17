// mex -O intTwoSimps.c

#include "mex.h"
#include "math.h"
#include <stdlib.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n, Nx, Ny, i, j;
    double a, b, c, d, e, f;
    double hx, hy, vv, x, y, weight;
    double *u, *si, *sj, *V, *out;
    double L1 = 0, L2 = 0;

    n = mxGetM(prhs[0]);
    u = mxGetPr(prhs[0]);
    si = mxGetPr(prhs[1]);
    sj = mxGetPr(prhs[2]);
    V = mxGetPr(prhs[3]);
    Nx = 2 * ceil((mxGetPr(prhs[4])[0]) / 2.0);
    Ny = 2 * ceil((mxGetPr(prhs[5])[0]) / 2.0);

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    out = mxGetPr(plhs[0]);

    hx = 1.0 / Nx;
    hy = 1.0 / Ny;

    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;
    e = 0.0;
    f = 0.0;
    for (i = 0; i < n; i++)
    {
        vv = V[i + n * i];
        f += u[i] * vv * u[i];
        d += u[i] * vv * si[i];
        e += u[i] * vv * sj[i];
        b += si[i] * vv * sj[i];
        a += si[i] * vv * si[i];
        c += sj[i] * vv * sj[i];
        L1 += si[i] * si[i];
        L2 += sj[i] * sj[i];
    }
    L1 = sqrt(L1);
    L2 = sqrt(L2);

    out[0] = 0.0;
    for (i = 0; i <= Nx; i++)
    {
        x = i * hx;
        for (j = 0; j <= Ny; j++)
        {
            if ((i == 0) || (i == Nx))
            {
                weight = 1.0;
            }
            else if (i % 2 == 0)
            {
                weight = 2.0;
            }
            else
            {
                weight = 4.0;
            }
            if ((j == 0) || (j == Ny))
            {
                weight = 1.0 * weight;
            }
            else if (j % 2 == 0)
            {
                weight = 2.0 * weight;
            }
            else
            {
                weight = 4.0 * weight;
            }

            y = j * hy;
            out[0] += weight * exp(-0.5 * (a * x * x - b * x * y - b * y * x + c * y * y + 2 * d * x - 2 * e * y + f));
        }
    }
    out[0] = out[0] * hx * hy / 9.0 * L1 * L2;
}
