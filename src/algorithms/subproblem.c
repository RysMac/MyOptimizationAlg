#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define T double
#define N 24  // Fixed size for Hessian and gradient

// Helper functions for matrix operations
void matrix_add_identity(T matrix[N][N], T lambda);
bool is_positive_definite(const T matrix[N][N]);
void cholesky_decomposition(const T matrix[N][N], T chol[N][N]);
void forward_substitution(const T chol[N][N], const T b[N], T y[N]);
void backward_substitution(const T chol[N][N], const T y[N], T x[N]);
T vector_norm(const T vector[N]);
void matrix_vector_product(const T matrix[N][N], const T vector[N], T result[N]);

// Core subproblem function
void subproblem(const T hess[N][N], const T grad[N], T sol[N], T delta, T tol, int max_iters, int verbose) {
    T regularized_hess[N][N];
    T chol[N][N];
    T temp[N];
    T lambda = 0.0;

    // Initialize solution to zero
    for (int i = 0; i < N; i++) {
        sol[i] = 0.0;
    }

    // Check if the Hessian is positive definite
    if (!is_positive_definite(hess)) {
        if (verbose > 1) {
            printf("Hessian is not positive definite.\n");
        }

        // Estimate minimum eigenvalue (simplified for performance)
        T eigen_min = 0.0;  // Replace this with actual eigenvalue computation if needed
        lambda = fabs(eigen_min) + 0.0001;

        if (verbose > 1) {
            printf("Regularizing with lambda = %.6f\n", lambda);
        }
    }

    // Main loop
    for (int iter = 0; iter < max_iters; iter++) {
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
        if (sol_norm <= delta) {
            if (lambda == 0 || fabs(sol_norm - delta) < tol * delta) {
                if (verbose > 0) {
                    printf("Newton step converged in iteration %d.\n", iter + 1);
                }
                return;
            } else {
                if (verbose > 1) {
                    printf("Hard case solution (not fully implemented).\n");
                }
                // Hard case handling can be added here if needed
                return;
            }
        } else {
            // Update lambda
            T ql_norm = vector_norm(temp);  // Use forward_substitution result as an approximation
            lambda += pow(sol_norm / ql_norm, 2) * ((sol_norm - delta) / delta);
        }
    }

    if (verbose > 0) {
        printf("Maximum iterations reached without convergence.\n");
    }
}

// Matrix operations
void matrix_add_identity(T matrix[N][N], T lambda) {
    for (int i = 0; i < N; i++) {
        matrix[i][i] += lambda;
    }
}

bool is_positive_definite(const T matrix[N][N]) {
    T chol[N][N] = {0};
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            T sum = matrix[i][j];
            for (int k = 0; k < j; k++) {
                sum -= chol[i][k] * chol[j][k];
            }
            if (i == j) {
                if (sum <= 0.0) {
                    return false;
                }
                chol[i][j] = sqrt(sum);
            } else {
                chol[i][j] = sum / chol[j][j];
            }
        }
    }
    return true;
}

void cholesky_decomposition(const T matrix[N][N], T chol[N][N]) {
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

void forward_substitution(const T chol[N][N], const T b[N], T y[N]) {
    for (int i = 0; i < N; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= chol[i][j] * y[j];
        }
        y[i] /= chol[i][i];
    }
}

void backward_substitution(const T chol[N][N], const T y[N], T x[N]) {
    for (int i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < N; j++) {
            x[i] -= chol[j][i] * x[j];
        }
        x[i] /= chol[i][i];
    }
}

T vector_norm(const T vector[N]) {
    T norm = 0.0;
    for (int i = 0; i < N; i++) {
        norm += vector[i] * vector[i];
    }
    return sqrt(norm);
}

// Example usage
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



#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#define N 24
#define MAX_ITERS 1000
#define TOL 1e-6

void mat_vec_mult(const double A[N][N], const double vec[N], double result[N]);
void normalize(double vec[N]);
double vector_norm(const double vec[N]);
void solve_linear_system(const double A[N][N], const double b[N], double x[N]); // Simple substitution for symmetric positive matrices

double inverse_power_method(const double A[N][N], double eigenvector[N], int max_iters, double tol) {
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
        // Solve (A - σI)x = v (σ = 0, so just solve Ax = v)
        solve_linear_system(A, v, Av);

        // Normalize the result to get the next eigenvector approximation
        normalize(Av);

        // Approximate the eigenvalue λ = (v^T A v) / (v^T v)
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

// Matrix-vector multiplication
void mat_vec_mult(const double A[N][N], const double vec[N], double result[N]) {
    for (int i = 0; i < N; i++) {
        result[i] = 0.0;
        for (int j = 0; j < N; j++) {
            result[i] += A[i][j] * vec[j];
        }
    }
}

// Normalize a vector
void normalize(double vec[N]) {
    double norm = vector_norm(vec);
    for (int i = 0; i < N; i++) {
        vec[i] /= norm;
    }
}

// Calculate the 2-norm of a vector
double vector_norm(const double vec[N]) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += vec[i] * vec[i];
    }
    return sqrt(sum);
}

// Solve a linear system A * x = b for a symmetric positive definite matrix
void solve_linear_system(const double A[N][N], const double b[N], double x[N]) {
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
