/*************************************************************
* AceGen    8.202 Linux (7 Oct 24)                           *
*           Co. J. Korelc  2020           22 Jan 25 00:32:35 *
**************************************************************
User     : Limited evaluation version
Notebook : Subproblem
Evaluation time                 : 2 s     Mode  : Optimal
Number of formulae              : 3       Method: Automatic
Subroutine                      : SubproblemAceGenTest size: 93
Total size of Mathematica  code : 93 subexpressions
Total size of C code            : 364 bytes */
#include "utilities.h"
#include "sms.h"

#include "stdlib.h"
#include "stdio.h"
#include "mathlink.h"
double workingvector[219];
void SubproblemAceGenTest(double v[119],double hess[6][6],double grad[6],double solution[6]);



void SubproblemAceGenTestMathLink(){
long ll;
int i1000,i1001,i1002,i1003,i1004,i1j1,i1j2,i2j1,i3j1,i1s1,i1s2,i2s1,i3s1;
const char *b1; double *b2;int *b3;
double hess[6][6];
double grad[6];
double solution[6];
++MathLinkCallCount[0];

/* read from link */
MLGetFunction(stdlink,&b1,&i1j1);
MLReleaseSymbol(stdlink,b1);
for(i1j1=0;i1j1<6;i1j1++){
MLGetRealList(stdlink,&b2,&ll);
for(i1j2=0;i1j2<6;i1j2++){
  hess[i1j1][i1j2]=b2[i1j2];
}
MLReleaseReal64List(stdlink,b2,ll);
};
MLGetRealList(stdlink,&b2,&ll);
for(i2j1=0;i2j1<6;i2j1++){
  grad[i2j1]=b2[i2j1];
}
MLReleaseReal64List(stdlink,b2,ll);
MLGetRealList(stdlink,&b2,&ll);
for(i3j1=0;i3j1<6;i3j1++){
  solution[i3j1]=b2[i3j1];
}
MLReleaseReal64List(stdlink,b2,ll);

/* allocate output parameters */
i1s1=6;
i1s2=6;
i2s1=6;
i3s1=6;

/* call module */
SubproblemAceGenTest(workingvector,hess,grad,solution);

/* write to link */
MLPutFunction(stdlink,"List",3);
PutRealArray(&b2,hess,i1s1,i1s2,MathLinkOptions[CO_SparseArray],0);
PutRealList(grad,i2s1);
PutRealList(solution,i3s1);
};

void MathLinkInitialize()
{
  MathLinkOptions[CO_NoSubroutines]=1;
  printf("MathLink module: %s\n","SubproblemAceGenTest");
};



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


#include "utilities.h"

// Matrix operations
void	matrix_add_identity(T matrix[N][N], T lambda) {
	for (int i = 0; i < N; i++) {
		matrix[i][i] += lambda;
	}
}

// Matrix-vector multiplication A_ij vec_j = result_i
void	mat_vec_mult(const T A[N][N], const T vec[N], T result[N]) {
	for (int i = 0; i < N; i++) {
		result[i] = 0.0;
		for (int j = 0; j < N; j++) {
			result[i] += A[i][j] * vec[j];
		}
	}
}

T v_dot_v(const T v[N]) {
	double  result = 0.0;
	for (int i = 0; i < N; i++) {
		result += v[i] * v[i]; // Sum of squares
	}
	return (result);
}

T v_dot_w(const T v[N], const T w[N]) {
	double  result = 0.0;
	for (int i = 0; i < N; i++) {
		result += v[i] * w[i]; // Sum of squares
	}
	return (result);
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
void solve_linear_system(const T A[N][N], const T b[N], T x[N]) {
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
double	inverse_power_method(const T A[N][N], T eigenvector[N], int max_iters, T tol) {
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
		// Solve (A - I"\203I)x = v (I"\203 = 0, so just solve Ax = v)
		solve_linear_system(A, v, Av);

		// Normalize the result to get the next eigenvector approximation
		normalize(Av);

		// Approximate the eigenvalue I^>> = (v^T A v) / (v^T v)
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

double	inverse_power_method_shift(const T A[N][N], T eigenvector[N], int max_iters, T tol, T sigma) {
	double v[N], Av[N], lambda = 0.0;
	double prev_lambda = 0.0;
	double diff = 0.0;
	double shifted_A[N][N];

	// Step 1: Create the shifted matrix (A - sigma * I)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			shifted_A[i][j] = A[i][j];
			if (i == j) {
				shifted_A[i][j] -= sigma; // Subtract sigma from diagonal elements
			}
		}
	}

	// Step 2: Initialize a random vector v (you can make this deterministic)
	for (int i = 0; i < N; i++) {
		v[i] = 1.0; // Start with all ones for simplicity
	}
	normalize(v);

	// Step 3: Iterative process
	for (int iter = 0; iter < max_iters; iter++) {
		// Solve (A - sigma * I)x = v
		solve_linear_system(shifted_A, v, Av);

		// Normalize the result to get the next eigenvector approximation
		normalize(Av);

		// Approximate the eigenvalue of the original matrix A
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

// find the biggest abs of eigenvalue and corresponding eigenvector
double	power_iteration(const T a[N][N], T eigenvector[N]) {
	double w[N];
	double v[N];
	double lambda = 0.0, lambda_old = 0.0;
	double epsilon = 1e-9;

	// Step 1: Initialize vector v with random values
	for (int i = 0; i < N; i++) {
		v[i] = 0.01 + 0.01 * i;  // Initialization, adjust as needed
	}

	// Step 2: Normalize the initial vector v
	normalize(v);

	// Step 3: Iterative process
	for (int j = 0; j < 1000; j++) {
		// Compute w = A * v
		mat_vec_mult(a, v, w);

		// Normalize w to get the next eigenvector approximation
		double wnorm = vector_norm(w);
		for (int i = 0; i < N; i++) {
			v[i] = w[i] / wnorm;  // Normalize w to become the next v
		}

		// Calculate the Rayleigh quotient: I^>> = (v^T A v) / (v^T v)
		mat_vec_mult(a, v, w);
		lambda_old = lambda;
		lambda = v_dot_w(v, w) / v_dot_v(v);

		// Check for convergence
		if ((lambda - lambda_old) * (lambda - lambda_old) < epsilon * epsilon) {
			break;
		}
	}

	// Store the resulting eigenvalue and eigenvector
	for (int i = 0; i < N; i++) {
		eigenvector[i] = v[i];
	}
	return lambda;
}

double	find_smallest_eigenvalue(const T A[N][N], T smallest_eigenvector[N]) {
	double	largest_eigenvector[N];
	double	shifted_A[N][N];

	// Step 1: Find the largest eigenvalue and eigenvector
	double largest_eigenvalue = power_iteration(A, largest_eigenvector);

	// Step 2: Create the shifted matrix A' = A - I^>>_max * I
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			shifted_A[i][j] = A[i][j];
			if (i == j) {
				shifted_A[i][j] -= largest_eigenvalue;
			}
		}
	}

	// Step 3: Use power_iteration on the shifted matrix to find the smallest eigenvalue of A'
	double	shifted_eigenvector[N];
	double	shifted_eigenvalue = power_iteration(shifted_A, shifted_eigenvector);

	// Step 4: Calculate the smallest eigenvalue of the original matrix
	double smallest_eigenvalue = shifted_eigenvalue + largest_eigenvalue;

	// Store the corresponding eigenvector
	for (int i = 0; i < N; i++) {
		smallest_eigenvector[i] = shifted_eigenvector[i];
	}
	return smallest_eigenvalue;
}

// find root with the use of Newton-Raphson method
double find_root(const T pl[N], const T eigvector[N], T delta, T initial_guess) {
	double x = initial_guess;
	double eigvector_norm_squared = 0.0;

	// Precompute ||eigvector||^2
	for (int i = 0; i < N; i++) {
		eigvector_norm_squared += eigvector[i] * eigvector[i];
	}

	for (int iter = 0; iter < ROOT_MAX_ITERS; iter++) {
		// Compute f(x) and f'(x)
		double norm_squared = 0.0, dot_product = 0.0;
		for (int i = 0; i < N; i++) {
			double temp = pl[i] + x * eigvector[i];
			norm_squared += temp * temp;
			dot_product += temp * eigvector[i];
		}

		double fx = sqrt(norm_squared) - delta;
		double fx_prime = dot_product / sqrt(norm_squared);

		// Avoid division by zero
		if (fabs(fx_prime) < ROOT_TOL/1000) {
			break;
		}

		// Newton-Raphson update
		double x_new = x - fx / fx_prime;

		// Check for convergence
		if (fabs(x_new - x) < ROOT_TOL) {
			return x_new; // Converged
		}

		x = x_new;
	}
	return x;
}


/******************* S U B R O U T I N E *********************/
void SubproblemAceGenTest(double v[119],double hess[6][6]
     ,double grad[6],double solution[6])
{
int i01;double v02;double sol[6];
v02=1e0;
subproblem(hess,grad,sol,&v02);
solution[0]=sol[0];
solution[1]=sol[1];
solution[2]=sol[2];
solution[3]=sol[3];
solution[4]=sol[4];
solution[5]=sol[5];
};
