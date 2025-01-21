#ifndef UTILITIES_H
# define UTILITIES_H

#include <stdio.h>
#include <math.h>

#define T double
#define N 24 // You can adjust this value as needed
#define MAX_ITERS 1000
#define TOL 1e-6

// Helper functions
void	matrix_add_identity(T matrix[N][N], T lambda);
T		vector_norm(const T vector[N]);

int		is_positive_definite(const T matrix[N][N]);
void	cholesky_decomposition(const T matrix[N][N], T chol[N][N]);
void	forward_substitution(const T chol[N][N], const T b[N], T y[N]);
void	backward_substitution(const T chol[N][N], const T y[N], T x[N]);
void	solve_linear_system(const T A[N][N], const T b[N], T x[N]);

T		inverse_power_method(const T A[N][N], T eigenvector[N], int max_iters, T tol);


#endif
