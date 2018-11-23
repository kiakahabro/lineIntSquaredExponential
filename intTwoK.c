// for OSX use 
// mex -O intTwoK.c -I/usr/local/include/ -lmwblas -lgsl
// for UNIX use 
// mex -O intTwoK.c -lmwblas -lgsl -lgslcblas -lm

// called from matlab as intTwoK(u,w,x,V)

#include "mex.h"
#include "blas.h"

#ifdef __APPLE__
    #include "libc.h"
#elif __unix__
    #include <unistd.h>
    #include <string.h>
    #include <pthread.h>
#else
    // other versions to be added
#endif

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>


#define DBL_EPS 2.220446049250313080847e-16

struct PARAMS {
    double a;
    double b;
    double c;
    double d;
    double e;
    double f;
};

double intfunc(double s,double a,double b, double c, double d,double e, double f);
double wrap_intfunc(double s,void * gsl_params);
void err_handler(const char * reason, const char * file, int line, int gsl_errno);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *u,*w,*x,*V;                /* pointers to input matrices */
    double *A,*R, *tmp;      /* working things*/
    double *OUT;
    ptrdiff_t m;                        /* matrix dimensions and lda = m*/
    int errcode = 0;            // for gsl functions
    ptrdiff_t n = 3;
    char *TRANS = "N";
    double one = 1.0; double zero = 0.0;  // for use in dgemm
    double c0,Ival,abserr;             // results
    size_t neval;               // number of evaluations used by solver
    ptrdiff_t ione = 1;


    /* Check for proper number of arguments. */
    if ( nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:intTwoK_mex:rhs",
            "This function requires 4 inputs.");
    }

    u = mxGetPr(prhs[0]); /* pointer to vector u */
    m = mxGetM(prhs[0]);  // dimension of problem
    w = mxGetPr(prhs[1]); // pointer to vector w
    x = mxGetPr(prhs[2]); // pointer to vector w
    V = mxGetPr(prhs[3]);// pointer to scaling matrix

    // check inputs
    if (mxGetM(prhs[1]) != m || mxGetM(prhs[2]) != m || mxGetN(prhs[0]) != 1 || mxGetN(prhs[1]) != 1 || mxGetN(prhs[2]) != 1){
        mexErrMsgIdAndTxt("MATLAB:intTwoK_mex:rhs",
            "u,w,x must be [m,1].");
    }
    // check sV and put into vector sv
    if (mxGetM(prhs[3]) != m || mxGetN(prhs[3]) != m){
        mexErrMsgIdAndTxt("MATLAB:intTwoK_mex:rhs",
            "sV must be [m,m].");
        }

    // check number of outputs
    if (nlhs != 1)
    {
        mexErrMsgIdAndTxt("MATLAB:qlessqr_lapack:lhs",
            "This function requires one output matrices");
    }

    // // create output matrices
    plhs[0] = mxCreateDoubleMatrix(4,1,mxREAL);
    OUT = mxGetPr(plhs[0]);

    // put matrices together
    A = malloc(m*n*sizeof(double));
    memcpy(A,w,m*sizeof(double));
    memcpy(A+m,x,m*sizeof(double));
    memcpy(A+2*m,u,m*sizeof(double));

    // calculte lengths
    double L1 = dnrm2(&m,w,&ione);
    double L2 = dnrm2(&m,x,&ione);

    // perform sV*[w,x,u]
    R = malloc(m*n*sizeof(double));
    dgemm(TRANS,TRANS,&m,&n,&m,&one,V,&m,A,&m,&zero,R,&m);

    double ww = ddot(&m,w,&ione,R,&ione);
    double xx = ddot(&m,x,&ione,R+m,&ione);
    double uu = ddot(&m,u,&ione,R+2*m,&ione);
    double b,c,d;
    if (ww < xx && xx > DBL_EPS){
    tmp = malloc(sizeof(double));
    memcpy(tmp,&ww,sizeof(double));
    memcpy(&ww,&xx,sizeof(double));
    memcpy(&xx,tmp,sizeof(double));
    b = -2*ddot(&m,u,&ione,R,&ione);
    c = -2*ddot(&m,u,&ione,R+m,&ione);
    d =2*ddot(&m,w,&ione,R+m,&ione);

    }else if (ww < DBL_EPS && xx < DBL_EPS){
        OUT[0] = L1*L2*exp(-uu/2);
        OUT[1] = 0;
        OUT[2] = -1;
        return;
    }else{
        b = 2*ddot(&m,u,&ione,R+m,&ione);
        c = 2*ddot(&m,u,&ione,R,&ione);
        d =2*ddot(&m,x,&ione,R,&ione);
    }
 
    // Calculate the constant value
    c0 = L1*L2*M_SQRTPI*M_SQRT1_2/sqrt(ww);


    // prepare for gsl
    gsl_function F;
    struct PARAMS params;
    params.a = uu;
    params.b = b;
    params.c = c;
    params.d = d;
    params.e = ww;
    params.f = xx;

    F.function = &wrap_intfunc;
    F.params = (void *) &params;

    // NON ADAPTIVE GAUSS CONROD
    gsl_set_error_handler_off();
    errcode= gsl_integration_qng(&F,0.0,1.0,sqrt(DBL_EPS),sqrt(DBL_EPS),&Ival,&abserr,&neval);
    // mexPrintf("nevals = %d \n",neval);

    OUT[0] = c0*Ival;
    OUT[1] = abserr;
    OUT[2] = (double) neval;
    OUT[3] = (double) errcode;

    free(R);
    free(A);
}


double intfunc(double s,double a,double b, double c, double d, double e, double f){
    double p1,p2;

    p1 = erf((c-d*s+2*e)/M_SQRT2/2/sqrt(e))-erf((c-d*s)/M_SQRT2/2/sqrt(e));
    p2 = exp(0.5*(-f*s*s+b*s-a+(c-d*s)*(c-d*s)/4/e));

    return p1*p2;
}

double wrap_intfunc(double s,void * gsl_params){
    struct PARAMS * params;
    params = (struct PARAMS *) gsl_params;


    return intfunc(s,params->a,params->b,params->c,params->d,params->e,params->f);
}

