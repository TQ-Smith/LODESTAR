
#include "Matrix.h"

#include <stdlib.h>

double** create_matrix(int n, int m) {

    double** matrix = (double**) calloc(n, sizeof(double**));

    for (int i = 0; i < n; i++)
        matrix[i] = (double*) calloc(m, sizeof(double*));
    
    return matrix;

}

void destroy_matrix(double** matrix, int n) {

    if (matrix == NULL)
        return;

    for (int i = 0; i < n; i++)
        free(matrix[i]);

    free(matrix);
    
}