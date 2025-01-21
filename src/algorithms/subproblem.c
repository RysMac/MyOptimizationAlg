
#include "utilities.h"


// subproblem function
void	subproblem(const T hess[N][N], const T grad[N], T sol[N], T *delta) {
	T regularized_hess[N][N];
	T chol[N][N];
	T temp[N];
	T lambda = 0.0;
	T smallest_eigenvector[N];


	//Initialize solution to zero
	for (int i = 0; i < N; i++) {
		sol[i] = 0.0;
	}

	// Check if the Hessian is positive definite
	if (!is_positive_definite(hess)) {
		// Estimate minimum eigenvalue (simplified for performance)
		T eigen_min = find_smallest_eigenvalue(hess, smallest_eigenvector);
		lambda = fabs(eigen_min) + 0.0001;
	}

	// Main loop
	for (int iter = 0; iter < subproblem_max_iters; iter++) {
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
		if (sol_norm <= *delta) {
			if (lambda == 0 || fabs(sol_norm - *delta) < subproblem_tol * *delta) {
				return;
			} else {
				lambda = find_root(sol, smallest_eigenvector, *delta, 0.0);
				// Hard case handling can be added here if needed
				return;
			}
		} else {
			// Update lambda
			T ql_norm = vector_norm(temp);  // Use forward_substitution result as an approximation
			lambda += pow(sol_norm / ql_norm, 2) * ((sol_norm - *delta) / *delta);
		}
	}
}
