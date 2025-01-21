#include "utilities.h"

// Example usage
int main() {
	double A[N][N] = {0};
	double eigenvector[N];

	// Initialize a test matrix (symmetric positive definite)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			A[i][j] = (i == j) ? 2.0 : 1.0; // Example: Diagonal dominance
		}
	}

	double tol = 1e-6;
	int max_iters = 1000;
	double min_eigenvalue = inverse_power_method(A, eigenvector, max_iters, tol);

	printf("Minimum Eigenvalue: %.6f\n", min_eigenvalue);
	printf("Eigenvector:\n");
	for (int i = 0; i < N; i++) {
		printf("%.6f ", eigenvector[i]);
	}
	printf("\n");

	return 0;
}
