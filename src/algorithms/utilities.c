#include "utilities.h"

// Matrix operations
void	matrix_add_identity(T matrix[N][N], T lambda) {
	for (int i = 0; i < N; i++) {
		matrix[i][i] += lambda;
	}
}

// Matrix-vector multiplication A_ij vec_j = result_i
void	mat_vec_mult(const double A[N][N], const double vec[N], double result[N]) {
	for (int i = 0; i < N; i++) {
		result[i] = 0.0;
		for (int j = 0; j < N; j++) {
			result[i] += A[i][j] * vec[j];
		}
	}
}

// vector 2-norm returns sqrt(v.v)
T	vector_norm(const T vector[N]) {
	T norm = 0.0;
	for (int i = 0; i < N; i++) {
		norm += vector[i] * vector[i];
	}
	return sqrt(norm);
}

// Normalize a vector
void normalize(double vec[N]) {
	double norm = vector_norm(vec);
	for (int i = 0; i < N; i++) {
		vec[i] /= norm;
	}
}

// basically this is Cholesky decomposition but returns 0 if failed
int	is_positive_definite(const T matrix[N][N]) {
	T chol[N][N] = {0};
	for (int i = 0; i < N; i++) {
		for (int j = 0; j <= i; j++) {
			T sum = matrix[i][j];
			for (int k = 0; k < j; k++) {
				sum -= chol[i][k] * chol[j][k];
			}
			if (i == j) {
				if (sum <= 0.0) {
					return 0;
				}
				chol[i][j] = sqrt(sum);
			} else {
				chol[i][j] = sum / chol[j][j];
			}
		}
	}
	return 1;
}

// Cholesky decomposition
void	cholesky_decomposition(const T matrix[N][N], T chol[N][N]) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j <= i; j++) {
			T sum = matrix[i][j];
			for (int k = 0; k < j; k++) {
				sum -= chol[i][k] * chol[j][k];
			}
			if (i == j) {
				chol[i][j] = sqrt(sum);
			} else {
				chol[i][j] = sum / chol[j][j];
			}
		}
	}
}

// Solves L.y = b where L is a lower triangular matrix
void	forward_substitution(const T chol[N][N], const T b[N], T y[N]) {
	for (int i = 0; i < N; i++) {
		y[i] = b[i];
		for (int j = 0; j < i; j++) {
			y[i] -= chol[i][j] * y[j];
		}
		y[i] /= chol[i][i];
	}
}

// Solves U.x = y where U is an upper triangular matrix
void	backward_substitution(const T chol[N][N], const T y[N], T x[N]) {
	for (int i = N - 1; i >= 0; i--) {
		x[i] = y[i];
		for (int j = i + 1; j < N; j++) {
			x[i] -= chol[j][i] * x[j];
		}
		x[i] /= chol[i][i];
	}
}

// Solve a linear system A * x = b for a symmetric positive definite matrix
void solve_linear_system(const double A[N][N], const double b[N], double x[N]) {
	double L[N][N] = {0}; // Cholesky factorization of A

	// Perform Cholesky decomposition: A = L * L^T
	for (int i = 0; i < N; i++) {
		for (int j = 0; j <= i; j++) {
			double sum = A[i][j];
			for (int k = 0; k < j; k++) {
				sum -= L[i][k] * L[j][k];
			}
			if (i == j) {
				L[i][j] = sqrt(sum);
			} else {
				L[i][j] = sum / L[j][j];
			}
		}
	}

	// Forward substitution: Solve L * y = b
	double y[N];
	for (int i = 0; i < N; i++) {
		y[i] = b[i];
		for (int j = 0; j < i; j++) {
			y[i] -= L[i][j] * y[j];
		}
		y[i] /= L[i][i];
	}

	// Backward substitution: Solve L^T * x = y
	for (int i = N - 1; i >= 0; i--) {
		x[i] = y[i];
		for (int j = i + 1; j < N; j++) {
			x[i] -= L[j][i] * x[j];
		}
		x[i] /= L[i][i];
	}
}

// Program for finding eigenvalue
double inverse_power_method(const double A[N][N], double eigenvector[N], int max_iters, double tol) {
	double v[N], Av[N], lambda = 0.0;
	double prev_lambda = 0.0;
	double diff = 0.0;

	// Step 1: Initialize a random vector v (you can make this deterministic)
	for (int i = 0; i < N; i++) {
		v[i] = 1.0; // Start with all ones for simplicity
	}
	normalize(v);

	// Step 2: Iterative process
	for (int iter = 0; iter < max_iters; iter++) {
		// Solve (A - σI)x = v (σ = 0, so just solve Ax = v)
		solve_linear_system(A, v, Av);

		// Normalize the result to get the next eigenvector approximation
		normalize(Av);

		// Approximate the eigenvalue λ = (v^T A v) / (v^T v)
		lambda = 0.0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				lambda += v[i] * A[i][j] * v[j];
			}
		}

		// Check for convergence
		diff = fabs(lambda - prev_lambda);
		if (diff < tol) {
			break;
		}

		// Update v and eigenvalue
		for (int i = 0; i < N; i++) {
			v[i] = Av[i];
		}
		prev_lambda = lambda;
	}

	// Store the resulting eigenvector
	for (int i = 0; i < N; i++) {
		eigenvector[i] = v[i];
	}

	return lambda;
}
