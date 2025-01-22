(*************************************************************
* AceGen    8.202 Linux (7 Oct 24)                           *
*           Co. J. Korelc  2020           22 Jan 25 19:28:46 *
**************************************************************
User     : Limited evaluation version
Notebook : Subproblem
Evaluation time                 : 2 s     Mode  : Debug
Number of formulae              : 48      Method: Automatic
Module                          : SubproblemAceGenTest size: 549
Total size of Mathematica  code : 549 subexpressions       *)


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

(*********************** M O D U L E **************************)
SetAttributes[SubproblemAceGenTest,HoldAll];
SubproblemAceGenTest[solution$$_]:=
Module[{i01,v04$$
     ,sol$$=Table[Null,{2}],v02$$=Table[Null,{2},{2}]
     ,v03$$=Table[Null,{2}]},
SMSExecuteBreakPoint["1","SubproblemAceGenTest",1,1];
$VV[1]=0;(*debug*)
(*2= hessVal_1|1 *)
$VV[2]=0.101400000000000022737367544323205947875976562500000000000000*10^4;
(*3= hessVal_1|2 *)
$VV[3]=-0.440000000000000056843418860808014869689941406250000000000000*10^3;
(*4= hessVal_2|1 *)
$VV[4]=-0.440000000000000056843418860808014869689941406250000000000000*10^3;
(*5= hessVal_2|2 *)
$VV[5]=0.200000000000000000000000000000000000000000000000000000000000*10^3;
(*6= gradVal_1 *)
$VV[6]=0.486000000000000511590769747272133827209472656250000000000000*10^2;
(*7= gradVal_2 *)
$VV[7]=-0.220000000000000213162820728030055761337280273437500000000000*10^2;
(*8= funcVal *)
$VV[8]=0.122000000000000219380069665930932387709617614746093750000000*10^1;
(*9= sol_1 *)
$VV[9]=0.110000000000000008881784197001252323389053344726562500000000*10^1;
(*10= sol_2 *)
$VV[10]=0.110000000000000008881784197001252323389053344726562500000000*10^1;
SMSExecuteBreakPoint["2","SubproblemAceGenTest",1,2];
$VV[11]=0;(*debug*)
Do[
 v02$$[[1,1]]=$VV[2];
 v02$$[[1,2]]=$VV[3];
 v02$$[[2,1]]=$VV[4];
 v02$$[[2,2]]=$VV[5];
 $VV[13]=0;(*debug*)
 v03$$[[1]]=$VV[6];
 v03$$[[2]]=$VV[7];
 $VV[14]=0;(*debug*)
 v04$$=0.100000000000000000000000000000000000000000000000000000000000*10^1;
 $VV[15]=0;(*debug*)
 subproblem[v02$$,v03$$,sol$$,v04$$];
 $VV[16]=0;(*debug*)
 (*17= solInner_1 *)
 $VV[17]=sol$$[[1]];
 (*18= solInner_2 *)
 $VV[18]=sol$$[[2]];
 (*19= solTry_1 *)
 $VV[19]=$VV[9]+$VV[17];
 (*20= solTry_2 *)
 $VV[20]=$VV[10]+$VV[18];
 SMSExecuteBreakPoint["3","SubproblemAceGenTest",1,3];
 $VV[21]=0;(*debug*)
 $VV[22]=Abs[-$VV[8]+(1-$VV[19])^2+100*(-$VV[19]^2+$VV[20])^2] < 1/10000000000000000 && Abs[
  -$VV[8]-$VV[6]*$VV[17]-$VV[7]*$VV[18]+-0.5*10^0*($VV[17]*($VV[2]*$VV[17]+$VV[4]*$VV[18])
  +$VV[18]*($VV[3]*$VV[17]+$VV[5]*$VV[18]))+(1-$VV[19])^2+100*(-$VV[19]^2+$VV[20])^2] < 1
  /10000000000000000;
 $VV[23]=$VV[22];
 If[$VV[23],
  (*24= rho *)
  $VV[24]=1;
  SMSExecuteBreakPoint["4","SubproblemAceGenTest",1,4];
  $VV[25]=0;(*debug*)
 , (* else *)
  (*24= rho *)
  $VV[24]=(-$VV[8]+(1-$VV[19])^2+100*(-$VV[19]^2+$VV[20])^2)/($VV[6]*$VV[17]+$VV[7]*$VV[18]
   +0.500000000000000000000000000000000000000000000000000000000000*10^0*($VV[17]*($VV[2]*$VV[17]
   +$VV[4]*$VV[18])+$VV[18]*($VV[3]*$VV[17]+$VV[5]*$VV[18])));
  SMSExecuteBreakPoint["5","SubproblemAceGenTest",1,5];
  $VV[26]=0;(*debug*)
 ]; (* endif *)
 SMSExecuteBreakPoint["6","SubproblemAceGenTest",1,6];
 $VV[27]=0;(*debug*)
 $VV[28]=$VV[24] > 0.1*10^-1;
 $VV[29]=$VV[28];
 If[$VV[29],
  (*9= sol_1 *)
  $VV[9]=$VV[19];
  $VV[32]=$VV[9]^2;
  $VV[30]=1-$VV[9];
  (*10= sol_2 *)
  $VV[10]=$VV[20];
  $VV[31]=$VV[10]-$VV[32];
  (*8= funcVal *)
  $VV[8]=$VV[30]^2+100*$VV[31]^2;
  (*6= gradVal_1 *)
  $VV[6]=-2*$VV[30]-400*$VV[9]*$VV[31];
  (*7= gradVal_2 *)
  $VV[7]=200*$VV[31];
  (*2= hessVal_1|1 *)
  $VV[2]=0.200000000000000000000000000000000000000000000000000000000000*10^1-400*$VV[10]
   +1200*$VV[32];
  (*3= hessVal_1|2 *)
  $VV[3]=-400*$VV[9];
  (*4= hessVal_2|1 *)
  $VV[4]=$VV[3];
  (*5= hessVal_2|2 *)
  $VV[5]=0.200000000000000000000000000000000000000000000000000000000000*10^3;
  SMSExecuteBreakPoint["7","SubproblemAceGenTest",1,7];
  $VV[33]=0;(*debug*)
 ]; (* endif *)
 SMSExecuteBreakPoint["8","SubproblemAceGenTest",1,8];
 $VV[34]=0;(*debug*)
 $VV[35]=Sqrt[Abs[$VV[6]]^2+Abs[$VV[7]]^2] < 1/1000000000000;
 $VV[36]=$VV[35];
 If[$VV[36],
  Print["Gradient = "," ",Sqrt[Abs[$VV[6]]^2+Abs[$VV[7]]^2]," ","  Outer iterations = "," "
   ,2];
  $VV[37]=0;(*debug*)
  Break[];
  $VV[38]=0;(*debug*)
  SMSExecuteBreakPoint["9","SubproblemAceGenTest",1,9];
  $VV[39]=0;(*debug*)
 ]; (* endif *)
 SMSExecuteBreakPoint["10","SubproblemAceGenTest",1,10];
 $VV[40]=0;(*debug*)
,{$VV[12],1,100,1}]; (*EndDo*)
solution$$[[1]]=$VV[9];
solution$$[[2]]=$VV[10];
$VV[41]=0;(*debug*)
SMSExecuteBreakPoint["11","SubproblemAceGenTest",1,11];
$VV[42]=0;(*debug*)
];
