
#include "MultidimensionalScaling.h"

/*
int
main()
{
  double *A, *B, *W, *Z, *WORK;
  int *ISUPPZ, *IWORK;
  int  i, j;
  int  M;

  A = malloc(N*N*sizeof(double));
  A[0] = 1; A[1] = 2; A[2] = 3; A[3] = 2; A[4] = 1; A[5] = 2; A[6] = 3; A[7] = 2; A[8] = 1;

  printf("A:\n");
  for (i=0; i<N; ++i) {
    for (j= i; j< N; ++j) {
      printf("%lf\t", A[INDEX(i, j)]);
    }
    printf("\n");
  }

  W = malloc(N*sizeof(double));
  Z = malloc(N*N*sizeof(double));
  ISUPPZ = malloc(2*N*sizeof(int));
  WORK = malloc(26*N*sizeof(double));
  IWORK = malloc(10*N*sizeof(int));

  dsyevr('V', 'I', 'U', N, A, N, 0, 0, 2, 3, dlamch('S'), &M,
         W, Z, N, ISUPPZ, WORK, 26*N, IWORK, 10*N);

  
  printf("%6.2f\t", Z[0]); printf("%6.2f\t", Z[3]); printf("%6.2f\t\n", Z[6]);
  printf("%6.2f\t", Z[1]); printf("%6.2f\t", Z[4]); printf("%6.2f\t\n", Z[7]);
  printf("%6.2f\t", Z[2]); printf("%6.2f\t", Z[5]); printf("%6.2f\t\n", Z[8]);

  printf("Eigenvalues:\n");
  for (i = 0; i < N; i++)
    printf("%6.2f\t", W[i]);
  printf("\n");

  return 0;
}
*/