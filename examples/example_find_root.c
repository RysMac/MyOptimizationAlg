// gcc -Wall -Wextra -Werror -I../src/algorithms example_find_root.c ../src/algorithms/utilities.c -o example_find_root -lm

#include "utilities.h"

int main() {
	double pl[N] = {1.0, 2.0, 3.0};       // Example vector
	double eigvector[N] = {0.5, 0.5, 0.5}; // Example eigenvector
	double delta = 2.5;                    // Example delta
	double initial_guess = 0.0;            // Initial guess

	double root = find_root(pl, eigvector, delta, initial_guess);
	printf("Root: %f\n", root);
	return 0;
}
