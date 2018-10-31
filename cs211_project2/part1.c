#include <stdio.h>                                                   
#include "cblas.h"                                                         
#include "lapacke.h"                                                       
#include <string.h> 
#include <iostream>                                                                           
using namespace std;                                                                           

int main (int argc, const char * argv[]) {                                                     
                                                                                               

  int m = 3;    
  int n = 3;    
  int lda = 3; 
  int ldb = 1;
  int sizes_for_m [] = {1000,2000,3000,4000,5000};


/* used to test for correctness in  algorithms

  A[0] = 4; A[1] = 6; A[2] = 3;  B[0] = 2;
  A[3] = 3; A[4] = 3; A[5] = 6;  B[1] = 10;                            
  A[6] = 1; A[7] = 4; A[8] = 1;  B[2] = 2;

*/


for(int i = 0; i <5; i++){
  m = sizes_for_m[i];

  n = m;    
  lda = m; 
  ldb = 1;

  double * A = new double[m*n]; 
  double * A_t = new double[m*n]; 
  double * B = new double[m];   
  int * ipiv = new int[m]; 

    
//Filling in A and B with random doubles
 for(int i=0; i<m; i++){
     for(int j=0; j<m; j++){
        A[i*m+j] = ((double)rand()*100)/((double)RAND_MAX + 0); //generate random double 0.00-100.0
     }
        B[i] = ((double)rand()*100)/((double)RAND_MAX + 0);
  }
  

clock_t time0;
time0 = clock();
    
  LAPACKE_dgetrf( LAPACK_ROW_MAJOR, m, n, A, lda, ipiv );   
/* used to print contents of matrices and vectors for debugging purposes:


    
      printf("\nfactorization:\n");
  for (int i = 0; i < m; i++) {  
          for (int j = 0; j < n; j++){ 
              printf("  %lf ", A[i*n+j]);                 
      }                                                                
      printf("\n");
  }  

 printf("\nipiv contents:\n");
for(int i = 0; i < m; i++)
  //  cout<<ipiv[i]<<" ";
  printf("  %d\n",ipiv[i]); 
  printf("\n");

*/

//doing operations that were one on pivot vector on the B vector. We have to do this because since rows were swapped in A and not B in dgetrf

for(int i = 0; i<m; i++)
    ipiv[i] -=1;

for(int i = 0; i<m; i++){
  double temp=0; 
  temp = B[i];
  B[i] = B[ipiv[i]];
  B[ipiv[i]] = temp;  
}


cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, m, 1, 1.0, A, lda, B, ldb);      

/*
  for(int i=0; i < n; i++){
     printf("  %lf ", B[i]);
     printf("\n");
  }
*/

 cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, m, 1, 1.0, A, lda, B, ldb);      

/*
  for(int i=0; i < n; i++){
     printf("  %lf ", B[i]);
     printf("\n");
  }
*/
time0 = clock()-time0;

printf("LAPACK library  went through  (%.15f seconds) to complete running for a matrix size of %dx%d. \n",((float)time0)/CLOCKS_PER_SEC, m,n);

}

  return 0;
}
