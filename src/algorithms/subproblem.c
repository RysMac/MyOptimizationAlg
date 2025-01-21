
#include "utilities.h"

// subproblem function
void	subproblem(const T hess[N][N], const T grad[N], T sol[N], T delta, T tol, int max_iters, int verbose) {
	T regularized_hess[N][N];
	T chol[N][N];
	T temp[N];
	T lambda = 0.0;

	// Initialize solution to zero
	for (int i = 0; i < N; i++) {
		sol[i] = 0.0;
	}

	// Check if the Hessian is positive definite
	if (!is_positive_definite(hess)) {
		if (verbose > 1) {
			printf("Hessian is not positive definite.\n");
		}

		// Estimate minimum eigenvalue (simplified for performance)
		T eigen_min = 0.0;  // Replace this with actual eigenvalue computation if needed
		lambda = fabs(eigen_min) + 0.0001;

		if (verbose > 1) {
			printf("Regularizing with lambda = %.6f\n", lambda);
		}
	}

	// Main loop
	for (int iter = 0; iter < max_iters; iter++) {
		// Regularize the Hessian: regularized_hess = 0.5 * (hess + hess^T) + lambda * I
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				regularized_hess[i][j] = 0.5 * (hess[i][j] + hess[j][i]);
			}
			regularized_hess[i][i] += lambda;
		}

		// Perform Cholesky decomposition: chol * chol^T = regularized_hess
		cholesky_decomposition(regularized_hess, chol);

		// Solve the Newton step: sol = regularized_hess \ -grad
		for (int i = 0; i < N; i++) {
			temp[i] = -grad[i];
		}

		// Forward substitution: solve chol * y = -grad
		T y[N];
		forward_substitution(chol, temp, y);

		// Backward substitution: solve chol^T * x = y
		backward_substitution(chol, y, sol);

		// Check norm of the solution
		T sol_norm = vector_norm(sol);
		if (sol_norm <= delta) {
			if (lambda == 0 || fabs(sol_norm - delta) < tol * delta) {
				if (verbose > 0) {
					printf("Newton step converged in iteration %d.\n", iter + 1);
				}
				return;
			} else {
				if (verbose > 1) {
					printf("Hard case solution (not fully implemented).\n");
				}
				// Hard case handling can be added here if needed
				return;
			}
		} else {
			// Update lambda
			T ql_norm = vector_norm(temp);  // Use forward_substitution result as an approximation
			lambda += pow(sol_norm / ql_norm, 2) * ((sol_norm - delta) / delta);
		}
	}

	if (verbose > 0) {
		printf("Maximum iterations reached without convergence.\n");
	}
}
