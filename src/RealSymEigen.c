
#include "RealSymEigen.h"

RealSymEigen* init_real_sym_eigen(int N) {
    RealSymEigen* eigen = malloc(sizeof(RealSymEigen));
    eigen -> A = malloc(N * N * sizeof(double));
    eigen -> W = malloc(N * sizeof(double));
    eigen -> Z = malloc(N * N * sizeof(double));
    eigen -> ISUPPZ = malloc(2 * N * sizeof(int));
    eigen -> WORK = malloc(26 * N * sizeof(double));
    eigen -> IWORK = malloc(10 * N * sizeof(int));
    eigen -> N = N;
    return eigen;
}

int compute_k_eigenpairs(RealSymEigen* eigen, int k) {
    return dsyevr('V', 'I', 'U', eigen -> N, eigen -> A, eigen -> N, 
                    0, 0, eigen -> N - k + 1, eigen -> N, dlamch('S'), 
                    &(eigen -> M), eigen -> W, eigen -> Z, eigen -> N, 
                    eigen -> ISUPPZ, eigen -> WORK, 26 * eigen -> N, 
                    eigen -> IWORK, 10 * eigen -> N);
}

int compute_k_eigenvalues(RealSymEigen* eigen, int k) {
    return dsyevr('N', 'I', 'U', eigen -> N, eigen -> A, eigen -> N, 
                    0, 0, eigen -> N - k + 1, eigen -> N, dlamch('S'), 
                    &(eigen -> M), eigen -> W, eigen -> Z, eigen -> N, 
                    eigen -> ISUPPZ, eigen -> WORK, 26 * eigen -> N, 
                    eigen -> IWORK, 10 * eigen -> N);
}

void destroy_real_sym_eigen(RealSymEigen* eigen) {
    if (eigen == NULL)
        return;
    free(eigen -> A);
    free(eigen -> W);
    free(eigen -> Z);
    free(eigen -> ISUPPZ);
    free(eigen -> WORK);
    free(eigen -> IWORK);
    free(eigen);
}