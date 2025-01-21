#include "utilities.h"

// Example usage of subproblem subroutine
int main() {
	T hess[N][N] = {0};
	T grad[N] = {0};
	T sol[N] = {0};

	// Fill hess and grad with test values
	for (int i = 0; i < N; i++) {
		grad[i] = i + 1;
		for (int j = 0; j < N; j++) {
			hess[i][j] = (i == j) ? 4.0 : 1.0;
		}
	}

	T delta = 10.0;
	T tol = 1e-4;
	int max_iters = 100;
	int verbose = 1;

	subproblem(hess, grad, sol, delta, tol, max_iters, verbose);

	printf("Solution:\n");
	for (int i = 0; i < N; i++) {
		printf("%.6f ", sol[i]);
	}
	printf("\n");

	return 0;
}
