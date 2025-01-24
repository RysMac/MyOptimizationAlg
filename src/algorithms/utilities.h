#ifndef UTILITIES_H
# define UTILITIES_H

#include <stdio.h>
#include <math.h>

#define T double
#define N 2 // You can adjust this value as needed
#define MAX_ITERS 1000
#define ROOT_MAX_ITERS 1000
#define subproblem_max_iters 100

#define TOL 1e-6
#define ROOT_TOL 1e-6
#define subproblem_tol 1e-4



// Helper functions
void	matrix_add_identity(T matrix[N][N], T lambda);
void	mat_vec_mult(const T A[N][N], const T vec[N], T result[N]);
T		v_dot_v(const T v[N]);
T		v_dot_w(const T v[N], const T w[N]);
T		vector_norm(const T vector[N]);
void 	normalize(double vec[N]);
int		is_positive_definite(const T matrix[N][N]);
void	cholesky_decomposition(const T matrix[N][N], T chol[N][N]);
void	forward_substitution(const T chol[N][N], const T b[N], T y[N]);
void	backward_substitution(const T chol[N][N], const T y[N], T x[N]);
void	solve_linear_system(const T A[N][N], const T b[N], T x[N]);

T	inverse_power_method(const T A[N][N], T eigenvector[N], int max_iters, T tol);
T	inverse_power_method_shift(const T A[N][N], T eigenvector[N], int max_iters, T tol, T sigma);
T	power_iteration(const T a[N][N], T eigenvector[N]);
T	find_smallest_eigenvalue(const T A[N][N], T smallest_eigenvector[N]);
T	find_root(const T pl[N], const T eigvector[N], T delta, T initial_guess);

void	subproblem(const T hess[N][N], const T grad[N], T sol[N], T *delta);
void	SubproblemAceGenTest(T v[359],T solution[2], T lb[2]);

#endif
