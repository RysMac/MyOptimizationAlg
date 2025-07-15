/*************************************************************
* AceGen    8.202 Linux (7 Oct 24)                           *
*           Co. J. Korelc  2020           31 May 25 00:09:21 *
**************************************************************
User     : Limited evaluation version
Notebook : CrystalPlasticity24SlipsPoint
Evaluation time                 : 3 s     Mode  : Optimal
Number of formulae              : 128     Method: Automatic
Subroutine                      : SubproblemAceGenTest size: 867
Total size of Mathematica  code : 867 subexpressions
Total size of C code            : 4615 bytes */
#include "utilities.h"
#include "sms.h"


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
			forward_substitution(chol, sol, y);
			T ql_norm = vector_norm(y);  // Use forward_substitution result as an approximation
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

// Function to scale a vector by a given factor
void vscale(double *z, int n, double factor) {
	for (int i = 0; i < n; i++) {
		z[i] *= factor;
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

// Function to estimate the smallest singular value and corresponding vector
double estsv(double **R, double *z, int n) {
	double e, s, sm, temp, w, wm, ynorm, znorm;

	// Initialize vector z to zero
	for (int i = 0; i < n; i++) {
		z[i] = 0.0;
	}

	// Set initial value of e based on R[0][0]
	e = fabs(R[0][0]);
	if (e == 0.0) {
		z[0] = 1.0;
		return 0.0;
	}

	// Solve R' * y = e
	for (int i = 0; i < n; i++) {
		e = copysign(e, -z[i]);

		if (fabs(e - z[i]) > fabs(R[i][i])) {
			temp = fmin(fabs(R[i][i] / (e - z[i])), SMALL);
			vscale(z, n, temp);
			e *= temp;
		}

		if (R[i][i] == 0.0) {
			w = 1.0;
			wm = 1.0;
		} else {
			w = (e - z[i]) / R[i][i];
			wm = -(e + z[i]) / R[i][i];
		}

		s = fabs(e - z[i]);
		sm = fabs(e + z[i]);

		for (int j = i + 1; j < n; j++) {
			sm += fabs(z[j] + wm * R[i][j]);
		}

		for (int j = i + 1; j < n; j++) {
			temp = z[j] + w * R[i][j];
			z[j] = temp;
			s += fabs(temp);
		}

		if (s < sm) {
			temp = wm - w;
			w = wm;
			if (temp != 0.0) {
				for (int j = i + 1; j < n; j++) {
					z[j] += temp * R[i][j];
				}
			}
		}
		z[i] = w;
	}
	ynorm = vector_norm(z);

	// Solve R * z = y
	for (int j = n - 1; j >= 0; j--) {
		if (fabs(z[j]) > fabs(R[j][j])) {
			temp = fmin(fabs(R[j][j] / z[j]), SMALL);
			vscale(z, n, temp);
			ynorm *= temp;
		}

		if (R[j][j] == 0.0) {
			z[j] = 1.0;
		} else {
			z[j] /= R[j][j];
		}

		if ((temp = z[j]) != 0.0) {
			for (int i = 0; i < j; i++) {
				z[i] -= temp * R[i][j];
			}
		}
	}

	// Normalize z and return svmin
	znorm = 1.0 / vector_norm(z);
	vscale(z, n, znorm);
	return ynorm * znorm;
}


/******************* S U B R O U T I N E *********************/
void SubproblemAceGenTest(double v[384],double solution[2]
     ,double lb[2])
{
int i01;double v04;double sol[2];double v02[2][2];double v03[2];
int i24,i52,b11,b14,b33,b34,b36,b37,b43,b45,b65,b66,b68,b69,b76,b79,b80,b82,b87
     ,b89,b96,b99,b101,b107,b110,b115,b119,b122,b130,b131,b133,b134,b139,b140,b142;
v[1]=solution[0];
v[157]=(v[1]*v[1]);
v[2]=solution[1];
v[158]=(v[2]*v[2]);
v[3]=lb[0];
v[4]=lb[1];
v[5]=100e0;
v[6]=v[1]-v[3];
v[7]=v[2]-v[4];
v[8]=v[5]*v[6];
v[9]=v[5]*v[7];
if(v[8]<0e0){
 v[12]=v[8];
} else {
 v[12]=0e0;
};
if(v[9]<0e0){
 v[15]=v[9];
} else {
 v[15]=0e0;
};
v[18]=v[12];
v[19]=v[15];
v[20]=Power(v[5],-0.1e0);
v[22]=v[1];
v[23]=v[2];
for(i24=1;i24<=100;i24++){
 v[150]=v[5]*v[7];
 v[149]=v[5]*v[6];
 v[25]=v[149]+v[18];
 b33=v[25]<0e0;
 v[26]=v[150]+v[19];
 b36=v[26]<0e0;
 v[29]=2e0+(b33?v[5]:0e0);
 v[30]=0e0;
 v[31]=0e0;
 v[32]=2e0+(b36?v[5]:0e0);
 if(b33){
  v[35]=v[25];
 } else {
  v[35]=0e0;
 };
 if(b36){
  v[38]=v[26];
 } else {
  v[38]=0e0;
 };
 v[41]=2e0*v[1]+v[35];
 v[42]=2e0*v[2]+v[38];
 if(b33){
  v[44]=(0.5e0*v[149]+v[18])*v[6];
 } else {
  v[44]=0e0;
 };
 if(b36){
  v[46]=(0.5e0*v[150]+v[19])*v[7];
 } else {
  v[46]=0e0;
 };
 v[48]=v[157]+v[158]+v[44]+v[46];
 v[49]=v[1];
 v[50]=v[2];
 v[51]=0.5e0;
 for(i52=1;i52<=10000;i52++){
  v02[0][0]=v[29];
  v02[0][1]=v[30];
  v02[1][0]=v[31];
  v02[1][1]=v[32];
  v03[0]=v[41];
  v03[1]=v[42];
  v04=v[51];
  subproblem(v02,v03,sol,&v04);
  v[57]=sol[0];
  v[153]=0.5e0*v[57];
  v[58]=sol[1];
  v[59]=v[49]+v[57];
  v[60]=v[50]+v[58];
  v[154]=(v[59]*v[59])+(v[60]*v[60]);
  v[61]=-v[3]+v[59];
  v[151]=v[5]*v[61];
  v[62]=-v[4]+v[60];
  v[152]=v[5]*v[62];
  v[63]=v[151]+v[18];
  v[64]=v[152]+v[19];
  b65=v[63]<0e0;
  if(b65){
   v[67]=(0.5e0*v[151]+v[18])*v[61];
  } else {
   v[67]=0e0;
  };
  b68=v[64]<0e0;
  if(b68){
   v[70]=(0.5e0*v[152]+v[19])*v[62];
  } else {
   v[70]=0e0;
  };
  v[73]=v[154]-v[48]+v[67]+v[70];
  v[74]=v[41]*v[57]+0.5e0*v[29]*(v[57]*v[57])+v[58]*(v[153]*v[30]+v[153]*v[31]+v[42]
   +0.5e0*v[32]*v[58]);
  if(fabs(v[73])<0.1e-15 && fabs(v[73]-v[74])<0.1e-15){
   v[77]=1e0;
  } else {
   v[77]=v[73]/v[74];
  };
  if(v[77]>0.2e0){
   v[49]=v[59];
   v[50]=v[60];
   if(b65){
    v[81]=v[61]*(v[18]+0.5e0*v[5]*v[61]);
   } else {
    v[81]=0e0;
   };
   if(b68){
    v[83]=v[62]*(v[19]+0.5e0*v[5]*v[62]);
   } else {
    v[83]=0e0;
   };
   v[48]=v[154]+v[81]+v[83];
   if(b65){
    v[88]=v[63];
   } else {
    v[88]=0e0;
   };
   if(b68){
    v[90]=v[64];
   } else {
    v[90]=0e0;
   };
   v[41]=2e0*v[49]+v[88];
   v[42]=2e0*v[50]+v[90];
   v[29]=2e0+(b65?v[5]:0e0);
   v[30]=0e0;
   v[31]=0e0;
   v[32]=2e0+(b68?v[5]:0e0);
  } else {
  };
  if(v[77]<0.25e0){
   v[51]=v[51]/4e0;
  } else {
  };
  printf("\n%s %g ","Delta = ",(double)(v[51]));
  if(v[77]>0.8e0 && fabs(v[51]-sqrt(Power(fabs(v[49]),2)+Power(fabs(v[50]),2)))>=(v[51]/100000000e0
   ) && v[51]<10e0){
   v[51]=2e0*v[51];
  } else {
  };
  if(sqrt(Power(fabs(v[41]),2)+Power(fabs(v[42]),2))<0.1e-11 || i52==10000){
   printf("\n%s %g %s %g ","Gradient = ",(double)(sqrt(Power(fabs(v[41]),2)+Power(fabs(v[42]),2)))
    ,"  Outer iterations = ",(double)(i52));
   break;
  } else {
  };
 };/* end for */
 v[22]=v[49];
 v[23]=v[50];
 v[104]=v[22]-v[3];
 v[105]=v[23]-v[4];
 if(v[22]<v[3]){
  v[108]=v[104];
 } else {
  v[108]=0e0;
 };
 if(v[23]<v[4]){
  v[111]=v[105];
 } else {
  v[111]=0e0;
 };
 if(sqrt(Power(fabs(v[108]),2)+Power(fabs(v[111]),2))<=v[20]){
  v[116]=v[18]+v[104]*v[5];
  v[117]=v[19]+v[105]*v[5];
  if(v[116]<0e0){
   v[120]=v[116];
  } else {
   v[120]=0e0;
  };
  if(v[117]<0e0){
   v[123]=v[117];
  } else {
   v[123]=0e0;
  };
  v[18]=v[120];
  v[19]=v[123];
  printf("\n%s %g %g ","Lgr update = ",(double)(v[18]),(double)(v[19]));
  v[20]=v[20]/Power(v[5],0.8e0);
 } else {
  v[5]=2e0*v[5];
  printf("\n%s %g ","Penalty update = ",(double)(v[5]));
  v[20]=Power(v[5],-0.1e0);
 };
 v[128]=v[18]+v[104]*v[5];
 v[129]=v[19]+v[105]*v[5];
 b130=v[128]<0e0;
 if(b130){
  v[132]=v[128];
 } else {
  v[132]=0e0;
 };
 b133=v[129]<0e0;
 if(b133){
  v[135]=v[129];
 } else {
  v[135]=0e0;
 };
 if(sqrt(Power(fabs(v[132]+2e0*v[22]),2)+Power(fabs(v[135]+2e0*v[23]),2))<=0.1e-11){
  if(b130){
   v[141]=v[128];
  } else {
   v[141]=0e0;
  };
  if(b133){
   v[143]=v[129];
  } else {
   v[143]=0e0;
  };
  printf("\n%s %g %s %g ","Gradient final= ",(double)(sqrt(Power(fabs(v[141]+2e0*v[22]),2)+Power
   (fabs(v[143]+2e0*v[23]),2))),"  Outer iterations = ",(double)(i24));
  break;
 } else {
 };
};/* end for */
solution[0]=v[22];
solution[1]=v[23];
};
