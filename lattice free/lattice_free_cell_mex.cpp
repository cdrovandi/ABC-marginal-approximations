#include <iostream>
#include <math.h>
#include <vector>
#include <mex.h>
#include <stdlib.h>
#include <cstring>

using namespace std;

// Stop after population
int N_catch = 4000;

// Times to output SS
const int ntimes = 4;
double outTimes[ntimes] = { 0,12,24,36 };

// Universal functions and constants
#define pi          3.14159265358979323846
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#define max(a,b)    (((a) < (b)) ? (b) : (a))

// Simulation parameters
double L = 1440;

// Constant physical parameters
double sigma2_m = pow(24, 2);
double sigma2_p = pow(24, 2);
double sigma2_b = pow(24, 2);

double mu_s = 24;
double sigma_d = 24;

// Pair correlation parameter
double PC_dr = 50;


// Uniform sample on (0,1]
double rand_scale = 1 / (double(RAND_MAX) + 1);
double sampleU() {
    double random = (rand() + 1) * rand_scale;
    return random;
}

// Distance2
double distance2(double x1, double x2, double y1, double y2, double L) {

	double dx = fabs(x1 - x2);
	double dy = fabs(y1 - y2);

	dx = min(dx, L - dx);
	dy = min(dy, L - dy);

	double d = pow(dx, 2) + pow(dy, 2);
	return d;

}

// Kernel
double kernel(double r2, double sigma2, double gamma) {

	double y = 0;

	if (r2 < 9 * sigma2) {

		y = gamma * exp(-r2 / (2 * sigma2));

	}

	return y;

}

// VM unnormalised density function
double VM(double x, double mu, double kappa) {

	double out = exp(kappa * cos(x - mu));
	return out;

}

// Von mises sample
double sampleVM(double mu, double kappa) {

	// Sample from non-normalised density with rejection sampling
	double fmax = VM(mu, mu, kappa);

	double x = sampleU() * 2 * pi - pi;
	double y = VM(x, mu, kappa);
	double u = sampleU()*fmax;

	while (u > y) {
		x = sampleU() * 2 * pi - pi;
		y = VM(x, mu, kappa);
		u = sampleU()*fmax;
	}

	return x;

}

// Choose agent
int choose_agent(vector<double> M, double M_tot) {

	int i = 0;
	double Mc = max(0, M[i]);
	double alpha = sampleU() * M_tot;
	while (alpha > Mc) {
		i += 1;
		Mc += max(0, M[i]);
	}
	return i;

}

// Periodic modulus
double mod(double x, double L) {

	if (x <= 0) {
		x += L;
	}
	else if (x >= L) {
		x -= L;
	}

	return x;

}

// 1D periodic displacement (x1 -> x2)
double disp(double x1, double x2, double L) {

	double s = 0;

	double dx = x2 - x1;
	double adx = fabs(dx);
	double dxL = L - adx;

	if (adx < dxL) {
		s = dx;
	}
	else {
		if (x1 < L / 2) {
			s = -dxL;
		}
		else {
			s = dxL;
		}
	}

	return s;

}


// Simulate function
int sim_obs(double params[], double *timeOffset, double timeOut, double *X, double *Y, int N) {
    int Nout = N;
	if (Nout == 0) {
		return(Nout);
	}

	double m = params[0];
	double p = params[1];

	double gamma_m = -m / 6 * exp(1 / 2);
	double gamma_p = -p / 6 * exp(1 / 2);
	double gamma_b = params[2];

	double t = *timeOffset;

	//vector<double> Xout(Nout);   // X locations
    //vector<double> Yout(Nout);   // Y locations
	vector<double> Xout(X, X + Nout);   // X locations
	vector<double> Yout(Y, Y + Nout);   // Y locations
	vector<double> M(Nout);  // Movement rates
	vector<double> P(Nout);  // Proliferation rates
	vector<double> Bx(Nout);  // Bias - x-component
	vector<double> By(Nout);  // Bias - y-component

	// Initialise rates
	for (int i = 0; i < Nout; i++) {

		double MrS = m;
		double PrS = p;
		double BxS = 0;
		double ByS = 0;

		double x1 = Xout[i];
		double y1 = Xout[i];

		for (int j = 0; j < Nout; j++) {
			if (i != j) {

				double x2 = Xout[j];
				double y2 = Yout[j];
				double r2 = distance2(x1, x2, y1, y2, L);

				double b_ = kernel(r2, sigma2_b, gamma_b);
				MrS += kernel(r2, sigma2_m, gamma_m);
				PrS += kernel(r2, sigma2_p, gamma_p);

				if (b_ != 0) {
					BxS += b_ * disp(x1, x2, L);
					ByS += b_ * disp(y1, y2, L);
				}

			}
		}

		M[i] = MrS;
		P[i] = PrS;
		Bx[i] = -BxS / sigma2_b;
		By[i] = -ByS / sigma2_b;

	}

	// SIMULATE THROUGH TIME
	while (t < timeOut && Nout != 0 && Nout < N_catch) {

		double M_tot = 0;
		double P_tot = 0;

		// Calculate total event rate
		for (int i = 0; i < Nout; i++) {

			M_tot += max(0, M[i]);
			P_tot += max(0, P[i]);

		}

		// Decide what should happen
		double alpha = sampleU() * (M_tot + P_tot);
		if (alpha < M_tot) {
			// Movement

			// Choose agent
			int i = choose_agent(M, M_tot);

			// Location of agent
			double xc = Xout[i];
			double yc = Yout[i];

			// Reduce rates of surrounding agents
			for (int j = 0; j < Nout; j++) {
				if (i != j) {

					double x2 = Xout[j];
					double y2 = Yout[j];

					double r2 = distance2(xc, x2, yc, y2, L);
					double m_ = kernel(r2, sigma2_m, gamma_m);
					double p_ = kernel(r2, sigma2_p, gamma_p);
					double b_ = kernel(r2, sigma2_b, gamma_b);
					if (m_ != 0) {
						M[j] -= m_;
					}
					if (p_ != 0) {
						P[j] -= p_;
					}
					if (b_ != 0) {
						b_ = -1 / sigma2_b * b_;
						Bx[j] -= b_ * disp(x2, xc, L);
						By[j] -= b_ * disp(y2, yc, L);
					}

				}
			}

			// Move somewhere
			double md = mu_s;

			// Include bias to get mu, kappa or VM
			double Bx_i = Bx[i];
			double By_i = By[i];
			double vm_mu = atan2(By_i, Bx_i);
			double vm_kappa = sqrt(pow(Bx_i, 2) + pow(By_i, 2));

			// Angle to move
			double theta = sampleVM(vm_mu, vm_kappa);

			// New position
			double xp = mod(xc + md * cos(theta), L);
			double yp = mod(yc + md * sin(theta), L);

			// Update rates
			double MrS = m;
			double PrS = p;
			double BxS = 0;
			double ByS = 0;

			for (int j = 0; j < Nout; j++) {
				if (i != j) {

					double x2 = Xout[j];
					double y2 = Yout[j];

					double r2 = distance2(xp, x2, yp, y2, L);
					double m_ = kernel(r2, sigma2_m, gamma_m);
					double p_ = kernel(r2, sigma2_p, gamma_p);
					double b_ = kernel(r2, sigma2_b, gamma_b);

					if (m_ != 0) {
						M[j] += m_;
						MrS += m_;
					}
					if (p_ != 0) {
						P[j] += p_;
						PrS += p_;
					}
					if (b_ != 0) {
						// Displacements go from xp -> x2
						b_ = -1 * b_ / sigma2_b;
						double sx = disp(xp, x2, L);
						double sy = disp(yp, y2, L);
						BxS += b_ * sx;
						ByS += b_ * sy;
						Bx[j] -= b_ * sx;
						By[j] -= b_ * sy;
					}

				}
			}

			// "Move" agent
			Xout[i] = xp;
			Yout[i] = yp;
			M[i] = MrS;
			P[i] = PrS;
			Bx[i] = BxS;
			By[i] = ByS;

		} else if (alpha < (M_tot + P_tot)) {
			// Proliferation

			// Choose agent
			int i = choose_agent(P, P_tot);

			// Location of agent
			double xc = Xout[i];
			double yc = Yout[i];

			// New location
			double u1 = sampleU();
			double u2 = sampleU();

			double xp = mod(xc + sigma_d * sqrt(-2 * log(u1)) * cos(2 * pi*u2), L);
			double yp = mod(yc + sigma_d * sqrt(-2 * log(u1)) * sin(2 * pi*u2), L);

			// Update rates
			double MrS = m;
			double PrS = p;
			double BxS = 0;
			double ByS = 0;

			for (int j = 0; j < Nout; j++) {

				double x2 = Xout[j];
				double y2 = Yout[j];

				double r2 = distance2(xp, x2, yp, y2, L);
				double m_ = kernel(r2, sigma2_m, gamma_m);
				double p_ = kernel(r2, sigma2_p, gamma_p);
				double b_ = kernel(r2, sigma2_b, gamma_b);

				if (m_ != 0) {
					M[j] += m_;
					MrS += m_;
				}
				if (p_ != 0) {
					P[j] += p_;
					PrS += p_;
				}
				if (b_ != 0) {

					b_ = -1 / sigma2_b * b_;

					// Displacements go from xp -> x2
					double sx = disp(xp, x2, L);
					double sy = disp(yp, y2, L);
					BxS += b_ * sx;
					ByS += b_ * sy;
					Bx[j] -= b_ * sx;
					By[j] -= b_ * sy;

				}

			}

			// "Create" new agent
			Xout.push_back(xp);
            Yout.push_back(yp);
            M.push_back(MrS);
            P.push_back(PrS);
            Bx.push_back(BxS);
            By.push_back(ByS);
            Nout  += 1;

		}

		// Sample timestep
		double tau = -log(sampleU()) / (M_tot + P_tot);
		t += tau;

	} // end through time loop

	for (int i=0; i<Nout; i++) {
		*(X + i) = Xout[i];
		*(Y + i) = Yout[i];
	}
	*timeOffset = t;
	
	return(Nout);

} // end sim_obs function



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int it, i, j;
	double t;

	// input
	double *params = (double *)mxGetData(prhs[0]);
	double *X = (double *)mxMalloc(N_catch * sizeof(double));
	double *Y = (double *)mxMalloc(N_catch * sizeof(double));
	double *X0 = (double *)mxGetData(prhs[1]);
	double *Y0 = (double *)mxGetData(prhs[2]);
	int N = (int)mxGetScalar(prhs[3]);
	int seed = (int)mxGetScalar(prhs[4]);
	memcpy(X, X0, N * sizeof(double));
	memcpy(Y, Y0, N * sizeof(double));

	// output
	mxArray *C = mxCreateCellMatrix(ntimes, 1);

	// set the seed
	srand(seed);

	double offset = 0;
	double *pOffset = &offset;
	
	for (it = 0;it < ntimes;it++) {
		t = outTimes[it];
		N = sim_obs(params, pOffset, t, X, Y, N);
		// assign to output
		mxArray *Mat;
		Mat = mxCreateDoubleMatrix(N, 2, mxREAL);
		memcpy(mxGetPr(Mat), X, sizeof(double)*N);
		memcpy(mxGetPr(Mat)+N, Y, sizeof(double)*N);
		mxSetCell(C, it, Mat);
	}
	plhs[0] = C;
	
	mxFree(X);
	mxFree(Y);
	
	return;
}