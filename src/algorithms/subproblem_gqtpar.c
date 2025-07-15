#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#define N 20
#define ZERO 0.0
#define ONE 1.0
#define HALF 0.5
#define SMALL 0.001

// Helper functions to compute norms and dot products
double nrm2(double vec[N], int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += vec[i] * vec[i];
    }
    return sqrt(sum);
}

double dot(double vec1[N], double vec2[N], int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += vec1[i] * vec2[i];
    }
    return sum;
}

// This is rpobably wrong function
void trsv(char uplo, char trans, char diag, double A[N][N], double x[N], int n) {
    // A simple triangular solver assuming A is triangular
    for (int i = 0; i < n; i++) {
        if (uplo == 'U') {
            for (int j = i + 1; j < n; j++) {
                x[i] -= A[i][j] * x[j];
            }
            if (diag == 'N') {
                x[i] /= A[i][i];
            }
        } else if (uplo == 'L') {
            for (int j = 0; j < i; j++) {
                x[i] -= A[i][j] * x[j];
            }
            if (diag == 'N') {
                x[i] /= A[i][i];
            }
        }
    }
}

void copy_matrix(double dest[N][N], double src[N][N], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

void zero_vector(double vec[N]) {
    for (int i = 0; i < N; i++) {
        vec[i] = 0.0;
    }
}

// Main function implementing the trust region algorithm
void gqtpar(double A[N][N], char uplo, double b[N], double delta, double rtol, double atol, int itmax, double par,
			double x[N], double z[N], double d[N], double w[N]) {

	double alpha	= 0.0;
	double anorm	= 0.0;
	double bnorm 	= 0.0;
	double prod 	= 0.0;
	double rxnorm 	= 0.0;
	double rznorm 	= 0.0;
	double xnorm 	= 0.0;
	double parc 	= 0.0;
	double parf 	= 0.0;
	double pars 	= 0.0;
	double parl 	= 0.0;
	double paru 	= 0.0;
	double temp 	= 0.0;
	double temp1 	= 0.0;
	double temp2 	= 0.0;
	double temp3 	= 0.0;

	int n = N;

	// Initialize arrays
	zero_vector(x);
	zero_vector(z);

	if (uplo == 'U') {
		// Copy upper triangular part of A to lower triangle and save diagonal
		for (int j = 0; j < n; j++) {
			d[j] = A[j][j];
			for (int i = j + 1; i < n; i++) {
				A[i][j] = A[j][i];
			}
		}
	}
	else if (uplo == 'L') {
		// Copy lower triangular part of A to upper triangle and save diagonal
		for (int j = 0; j < n; j++) {
			d[j] = A[j][j];
			for (int i = j + 1; i < n; i++) {
				A[j][i] = A[i][j];
			}
		}
	}
	else {
		printf("Invalid value for uplo='%c'\n", uplo);
		return;
	}

	// Calculate the norms of A and b
	for (int j = 0; j < n; j++) {
		temp = 0.0;
		for (int i = 0; i < n; i++) {
			temp += fabs(A[i][j]);
		}
		anorm = fmax(anorm, temp);
		w[j] = temp - fabs(d[j]);
	}

	// Calculate a lower bound, PARS, for the domain of the problem.
	// Also calculate an upper bound, PARU, and a lower bound, PARL,
	// for the Lagrange multiplier.
	pars = parl = paru = -anorm;
	for (int j = 0; j < n; j++) {
		pars = fmax(pars, -d[j]);
		parl = fmax(parl, w[j] + d[j]);
		paru = fmax(paru, w[j] - d[j]);
	}

	temp = bnorm / delta;
	parl = fmax(fmax(0.0, temp - parl), pars);
	paru = fmax(0.0, temp + paru);

	par = fmin(fmax(par, parl), paru);
	paru = fmax(paru, (1.0 + rtol) * parl);

	int info = 0;
	int iter = 0;

	// Begin iteration
	while (1) {
		// Safeguard par
		if (par <= pars && paru > 0.0) {
			par = fmax(sqrt(parl / paru), SMALL) * paru;
		}

		// Compute A + par * I and attempt Cholesky factorization
		for (int j = 0; j < n; j++) {
			A[j][j] = d[j] + par;
			for (int i = j + 1; i < n; i++) {
				A[j][i] = A[i][j];
			}
		}


		int		indef = is_positive_definite(A);
		bool	rednc = false;

		if (indef == 0) {
			parf = par;
			for (int j = 0; j < n; j++) {
				w[j] = b[j];
			}

			// trsv function is probably wrong
			trsv('U', 'T', 'N', A, w, n);
			rxnorm = nrm2(w, n);
			trsv('U', 'N', 'N', A, w, n);
			for (int j = 0; j < n; j++) {
				x[j] = -w[j];
			}
			xnorm = nrm2(x, n);

			if (fabs(xnorm - delta) <= rtol * delta
				|| (par == 0.0 && xnorm <= (1.0 + rtol) * delta)) {
				info = 1;
			}

			// tu moze byc problem bo A must be pointer
			rznorm = estsv(A, z, N) ;  // Negative curvature direction
			pars = fmax(pars, par - rznorm * rznorm);

			// Negative curvature solution: x + alpha * z, where norm(x + alpha * z) = delta
			if (xnorm < delta) {
				prod = dot(z, x, n) / delta;
				temp = ((delta + xnorm) / delta) * (delta - xnorm);
				alpha = temp / (fabs(prod) + sqrt(prod * prod + temp / delta));
				alpha = copysign(alpha, prod);

				rznorm *= fabs(alpha);
				temp1 = (rznorm / delta) * (rznorm / delta);
				rednc = (temp1 + par * (xnorm / delta) * (xnorm / delta) <= par);

				temp2 = par + (rxnorm / delta) * (rxnorm / delta);
				if (HALF * temp1 <= rtol * (1.0 - HALF * rtol) * temp2) {
					info = 1;
				}
				if (info == 0 && HALF * temp2 <= (atol / delta) / delta) {
					info = 2;
				}
				if (xnorm == 0.0) {
					info = 1;
				}
			}

			// Compute the Newton correction for par
			if (xnorm == 0.0) {
				parc = -par;
			} else {
				temp = 1.0 / xnorm;
				// be careful here also
				for (int j = 0; j < n; j++) {
					w[j] = x[j] * temp;
				}
				trsv('U', 'T', 'N', A, w, n);
				temp = nrm2(w, n);
				parc = (((xnorm - delta) / delta) / temp) / temp;
			}

			// Update PARL or PARU.
			if (xnorm > delta && parl < par){
				parl = par;
			}
			if (xnorm < delta && paru > par){
				paru = par;
			}

		} else {
			// Case where A + par * I is not positive definite, handle accordingly
			// (this part can involve updating the matrix A and solving again)
			// ###################################################
			//  ##  Case 2: A + par⋅I is not positive definite.  ##
			// ###################################################

			// Use the rank information from the Cholesky
			// decomposition to update PAR.

			if (indef > 1) {
				// Restore column indef to A + par * I
				for (int j = 0; j < indef - 1; j++) {
					A[j][indef] = A[indef][j];
				}
				A[indef][indef] = d[indef] + par;

				// Compute PARC
				for (int j = 0; j < indef - 1; j++) {
					w[j] = A[j][indef];
				}
				trsv('U', 'T', 'N', indef - 1, A, w);
				// Update A[indef][indef] based on w
				A[indef][indef] -= nrm2(indef - 1, w, 1) * nrm2(indef - 1, w, 1);
				trsv('U', 'N', 'N', indef - 1, A, w);
			}

			// Now handle the rest of the operations for parc and paru
			w[indef] = -ONE;
			double temp = nrm2(indef, w, 1);
			parc = -(A[indef][indef] / temp) / temp;
			pars = fmax(pars, fmax(par, par + parc));

			// Adjust PARU if necessary
			paru = fmax(paru, (ONE + rtol) * pars);
		}
	}
