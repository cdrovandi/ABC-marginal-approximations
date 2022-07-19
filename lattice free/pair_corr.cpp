#include <math.h>
#include <mex.h>
#include <cstring>
#include <matrix.h>

// Universal functions and constants
#define pi          3.14159265358979323846
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#define max(a,b)    (((a) < (b)) ? (b) : (a))

// Simulation parameters
double L    = 1440;

// Distance2
double distance2(double x1, double x2, double y1, double y2, double L) {

    double dx = fabs(x1 - x2);
    double dy = fabs(y1 - y2);

    dx = min(dx,L - dx);
    dy = min(dy,L - dy);

    double d = pow(dx,2) + pow(dy,2);
    return d;

}

// Pair correlation function
double pair_corr(double *X, double *Y, size_t N, double PC_dr) {

    double pair_corr = 0.0;


    // Loop through all agents
    for (unsigned int i = 0; i < N; i++) {

        // Get big loop agent location
        double x0 = *(X + i);
        double y0 = *(Y + i);

        // Store distances of other agents
        for (unsigned int j = 0; j < N; j++) { if (i != j) {
            double x1 = *(X + j);
            double y1 = *(Y + j);
            if (sqrt(distance2(x0,x1,y0,y1,L)) < PC_dr) {
                pair_corr += 1;
            }
        }}


    } // End loop through all agents

    double density = N / pow(L,2);

    // Calculate pair correlation
    pair_corr /= N * pi * density * PC_dr * PC_dr;

    return pair_corr;

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
    // input
    double *data = (double *) mxGetData(prhs[0]);
    double PC_dr = (double) mxGetScalar(prhs[1]);
    size_t N = mxGetM(prhs[0]);
    double *X = (double *) mxMalloc(N * sizeof(double));
    double *Y = (double *) mxMalloc(N * sizeof(double));
    memcpy(X, data, N * sizeof(double));
    memcpy(Y, data + N, N * sizeof(double));
    
    // output
	double out = pair_corr(X, Y, N, PC_dr);
    plhs[0] = mxCreateDoubleScalar(out);
    
	mxFree(X);
	mxFree(Y);
	
    return;
}