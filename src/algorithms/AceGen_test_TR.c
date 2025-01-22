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
	double	initial[] = {-3, -4};
	double v[197];

	SubproblemAceGenTest(v, initial);
	printf("\n");
	printVector(initial, 2);

}
