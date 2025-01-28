#include "utilities.h"

void printVector(const double *vector, int size) {
    printf("[");
    for (int i = 0; i < size; i++) {
        printf("%f", vector[i]); // Use %f for floating-point numbers
        if (i < size - 1) {
            printf(", "); // Add a comma and space between elements
        }
    }
    printf("]\n"); // Close the vector and go to a new line
}


int	main()
{
	double	initial[2];
	double	lb[] = {5., 2.};
	double v[359];

	for (int i = 0; i < 1; i++) {
		initial[0] = (double) 20.0;
		initial[1] = (double) 2.0;
		SubproblemAceGenTest(v, initial, lb);
		printf("\n");
		printVector(initial, 2);
	}
}
