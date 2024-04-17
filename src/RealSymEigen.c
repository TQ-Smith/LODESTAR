
#include "RealSymEigen.h"

extern void dsyevr_(char *JOBZp, char *RANGEp, char *UPLOp, 
                        int *Np, double *A, int *LDAp, double *VLp, 
                        double *VUp, int *ILp, int *IUp, double *ABSTOLp, 
                        int *Mp, double *W, double *Z, int *LDZp, int *ISUPPZ, 
                        double *WORK, int *LWORKp, int *IWORK, int *LIWORKp, int *INFOp);

extern double dlamch_(char *CMACHp);

int dsyevr(
    char JOBZ, 
    char RANGE, 
    char UPLO, 
    int N,
    double *A, 
    int LDA, 
    double VL, 
    double VU,
    int IL, 
    int IU, 
    double ABSTOL, 
    int *M,
    double *W, 
    double *Z, 
    int LDZ, 
    int *ISUPPZ,
    double *WORK, 
    int LWORK, 
    int *IWORK, 
    int LIWORK
) {

    int INFO;

    dsyevr_(&JOBZ, &RANGE, &UPLO, 
            &N, A, &LDA, &VL, 
            &VU, &IL, &IU, &ABSTOL,
            M, W, Z, &LDZ, ISUPPZ, 
            WORK, &LWORK, IWORK, &LIWORK, &INFO);
    
    return INFO;

}

double dlamch(char CMACH) {
  return dlamch_(&CMACH);
}


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
    if (k == eigen -> N)
        return dsyevr('V', 'A', 'U', eigen -> N, eigen -> A, eigen -> N, 
                    0, 0, 0, eigen -> N, dlamch('S'), 
                    &(eigen -> M), eigen -> W, eigen -> Z, eigen -> N, 
                    eigen -> ISUPPZ, eigen -> WORK, 26 * eigen -> N, 
                    eigen -> IWORK, 10 * eigen -> N);

    return dsyevr('V', 'I', 'U', eigen -> N, eigen -> A, eigen -> N, 
                    0, 0, eigen -> N - k + 1, eigen -> N, dlamch('S'), 
                    &(eigen -> M), eigen -> W, eigen -> Z, eigen -> N, 
                    eigen -> ISUPPZ, eigen -> WORK, 26 * eigen -> N, 
                    eigen -> IWORK, 10 * eigen -> N);
}

int compute_k_eigenvalues(RealSymEigen* eigen, int k) {
    if (k == eigen -> N)
        return dsyevr('N', 'A', 'U', eigen -> N, eigen -> A, eigen -> N, 
                    0, 0, 0, eigen -> N, dlamch('S'), 
                    &(eigen -> M), eigen -> W, eigen -> Z, eigen -> N, 
                    eigen -> ISUPPZ, eigen -> WORK, 26 * eigen -> N, 
                    eigen -> IWORK, 10 * eigen -> N);
    
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