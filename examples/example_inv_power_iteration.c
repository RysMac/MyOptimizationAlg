#include "utilities.h"

// Example usage
int main() {
	//double A[N][N] = {0};
	double eigenvector[N];

	// Initialize a test matrix (symmetric positive definite)
	//for (int i = 0; i < N; i++) {
	//	for (int j = 0; j < N; j++) {
	//		A[i][j] = (i == j) ? 2.0 : 1.0; // Example: Diagonal dominance
	//	}
	//}

	// Example 5x5 matrix
	// double A[5][5] = {
	// 	{2.0, 1.0, 0.0, 0.0, 0.0},
	// 	{1.0, 2.0, 1.0, 0.0, 0.0},
	// 	{0.0, 1.0, 2.0, 1.0, 0.0},
	// 	{0.0, 0.0, 1.0, 2.0, 1.0},
	// 	{0.0, 0.0, 0.0, 1.0, 2.0}
	// };
	double A[5][5] =	{{4.0, 1.0, 0.0, 0.0, 0.0},
						{1.0, -2.0, 1.0, 0.0, 0.0},
						{0.0, 1.0, 3.0, 1.0, 0.0},
						{0.0, 0.0, 1.0, 0.0, 1.0},
						{0.0, 0.0, 0.0, 1.0, -1.0}};

	double tol = 1e-6;
	int max_iters = 1000;
	double min_eigenvalue = inverse_power_method(A, eigenvector, max_iters, tol);

	printf("inverse_power_method\n");
	printf("Minimum Eigenvalue: %.6f\n", min_eigenvalue);
	printf("Eigenvector:\n");
	for (int i = 0; i < N; i++) {
		printf("%.6f ", eigenvector[i]);
	}
	printf("\n");

	min_eigenvalue = inverse_power_method_shift(A, eigenvector, max_iters, tol, -4.19);

	printf("inverse_power_method_shift\n");
	printf("Minimum Eigenvalue: %.6f\n", min_eigenvalue);
	printf("Eigenvector:\n");
	for (int i = 0; i < N; i++) {
		printf("%.6f ", eigenvector[i]);
	}
	printf("\n");

	min_eigenvalue = power_iteration(A, eigenvector);

	printf("power_iteration\n");
	printf("Maximum Eigenvalue: %.6f\n", min_eigenvalue);
	printf("Eigenvector:\n");
	for (int i = 0; i < N; i++) {
		printf("%.6f ", eigenvector[i]);
	}
	printf("\n");


	min_eigenvalue = find_smallest_eigenvalue(A, eigenvector);

	printf("find_smallest_eigenvalue\n");
	printf("Minimum Eigenvalue: %.6f\n", min_eigenvalue);
	printf("Eigenvector:\n");
	for (int i = 0; i < N; i++) {
		printf("%.6f ", eigenvector[i]);
	}
	printf("\n");

	return 0;
}
