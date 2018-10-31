#include <iostream>
#include "lapacke.h"
#include "cblas.h"
#include <unistd.h>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>

//declearation
int mydgetrf(int row, int col, double *a, int lda, int *ipiv);
int Blocked_dgetrf(int row, int col, double *a, int lda, int *ipiv, int block_size);
int mydtrsm(char trans, int n, int nrhs, double *a, int lda, int* ipiv, double *b, int ldb);

//HDdiff is from Stackoverflow
struct timespec HDdiff(struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec - start.tv_nsec) < 0) {
        temp.tv_sec = end.tv_sec - start.tv_sec - 1;
        temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
    }
    else {
        temp.tv_sec = end.tv_sec - start.tv_sec;
        temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    }
    return temp;
}

struct timespec HDadd(struct timespec a, struct timespec b)
{
    struct timespec temp;
    temp.tv_nsec = (a.tv_nsec + b.tv_nsec) % 1000000000;
    temp.tv_sec = a.tv_sec + b.tv_sec + (a.tv_nsec + b.tv_nsec) / 1000000000;
    return temp;
}

//Sqeuence to time the LAPACK perfermance
void testLapack(double *a, double *b, int n) {
    struct timespec begin, end, diff;
    int lda = n, info = 3;
    int *ipiv = new int[n];

    clock_gettime(CLOCK_MONOTONIC, &begin);

    info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, a, lda, ipiv);
    if (info != 0) {
        std::cout << "LAPACKE_dgetrf FAILED" << std::endl;
        return;
    }

    char TRANS = 'N';
    int m = 1;
    info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, TRANS, n, m, a, n, ipiv, b, n);
    if (info != 0) {
        std::cout << "LAPACKE_dgetrs FAILED" << std::endl;
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    diff = HDdiff(begin, end);
    printf("LAPACK, n=%d, Time:%ld seconds and %ld nanoseconds.\n", n, diff.tv_sec, diff.tv_nsec);
}

//Sqeuence to time the my algorithm perfermance
void testMine(double *a, double *b, int n) {
    struct timespec begin, end, diff;
    int lda = n, info = 3;
    int *ipiv = new int[n];

    clock_gettime(CLOCK_MONOTONIC, &begin);

    info = mydgetrf(n, n, a, lda, ipiv);
    if (info != 0) {
        std::cout << "mydgetrf FAILED" << std::endl;
        return;
    }

    char TRANS = 'N';
    int m = 1;
    info = mydtrsm(TRANS, n, m, a, n, ipiv, b, n);
    if (info != 0) {
        std::cout << "mydtrsm FAILED" << std::endl;
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    diff = HDdiff(begin, end);
    printf("My algorithm, n=%d, Time:%ld seconds and %ld nanoseconds.\n", n, diff.tv_sec, diff.tv_nsec);
}

//Sqeuence to time the my GEPP algorithm perfermance
void testBlcoked(double *a, double *b, int n, int block_size) {
    struct timespec begin, end, mid, diff;
    int lda = n, info = 3;
    int *ipiv = new int[n];

    clock_gettime(CLOCK_MONOTONIC, &begin);
    //=============================================================================================================LOOK<-HERE
    info = Blocked_dgetrf(n, n, a, lda, ipiv, block_size);
    if (info != 0) {
        std::cout << "mydgetrf FAILED" << std::endl;
        return;
    }
    clock_gettime(CLOCK_MONOTONIC, &mid);
    char TRANS = 'N';
    int m = 1;
    info = mydtrsm(TRANS, n, m, a, n, ipiv, b, n);
    if (info != 0) {
        std::cout << "mydtrsm FAILED" << std::endl;
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    diff = HDdiff(begin, mid);
    printf("My Blocked GEPP, First Half, n=%d, Time:%ld seconds and %ld nanoseconds.\n", n, diff.tv_sec, diff.tv_nsec);
    diff = HDdiff(mid, end);
    printf("My Blocked GEPP, Second Half, n=%d, Time:%ld seconds and %ld nanoseconds.\n", n, diff.tv_sec, diff.tv_nsec);
    diff = HDdiff(begin, end);
    printf("My Blocked GEPP, n=%d, Time:%ld seconds and %ld nanoseconds.\n", n, diff.tv_sec, diff.tv_nsec);
}

//Unblocked dgetrf
int mydgetrf(int row, int col, double *a, int lda, int *ipiv) {
    int n = row;
    if (n != col) {
        std::cout << "ERROR, ONLY SUPPORT REGTANGLE MATRIX" << std::endl;
        return -1;
    }
    ipiv[n - 1] = n;
    for (int i = 0; i < n - 1; ++i) {
        int maxp = i;
        double max = std::abs(a[i*n + i]);
        for (int t = i + 1; t < n; ++t) {
            if (std::abs(a[i*n + t]) > max) {
                maxp = t;
                max = std::abs(a[i*n + t]);
            }
        }
        if (max == 0) {
            std::cout << "LUfactoration failed: coefficient matrix is singular" << std::endl;
            return -1;
        }
        else {
            //save pivoting infomation in LAPACK format
            ipiv[i] = maxp + 1;
            if (maxp != i) {
                //swap rows
                for (int j = 0; j < n; ++j) {
                    double tmp = a[j*n + i];
                    a[j*n + i] = a[j*n + maxp];
                    a[j*n + maxp] = tmp;
                }
            }
        }
        for (int j = i + 1; j < n; ++j) {
            a[i*n + j] /= a[i*n + i];
            for (int k = i + 1; k < n; k++) {
                a[k*n + j] -= a[i*n + j] * a[k*n + i];
            }
        }
    }
    return 0;
}

//Blocked dgetrf
int Blocked_dgetrf(int row, int col, double *a, int lda, int *ipiv, int block_size) {
    int info = 0;
    if (row != col) {
        std::cout << "ERROR, ONLY SUPPORT REGTANGLE MATRIX" << std::endl;
        return -1;
    }

    struct timespec st1, st2, st3, st4, st5;
    st1.tv_nsec = 0;
    st1.tv_sec = 0;
    st2 = st1;
    st3 = st1;
    st4 = st1;
    st5 = st1;

    for (int p = 0; p < row; p += block_size) {
        int pb = std::min(row - p, block_size);
        //struct timespec begin, end, cp1, cp2, cp3, cp4, diff;
        //clock_gettime(CLOCK_MONOTONIC, &begin);
        ////DGETRF2
        int rowToGo = row - p;
        int colToGo = pb;
        ipiv[row - 1] = row;//p+row-1=row-1
        for (int i = p; i < std::min(row - 1, p + pb); ++i) {
            int maxp = i;
            double max = std::abs(a[i*lda + i]);
            for (int t = i + 1; t < row; ++t) {
                if (std::abs(a[i*lda + t]) > max) {
                    maxp = t;
                    max = std::abs(a[i*lda + t]);
                }
            }
            if (max == 0) {
                std::cout << "LUfactoration failed: coefficient matrix is singular" << std::endl;
                return -1;
            }
            else {
                //save pivoting infomation in LAPACK format
                ipiv[i] = maxp + 1;
                if (maxp != i) {
                    //swap rows
                    for (int j = p; j < p + colToGo; ++j) {
                        double tmp = a[j*lda + i];
                        a[j*lda + i] = a[j*lda + maxp];
                        a[j*lda + maxp] = tmp;
                    }
                }
            }
            for (int j = i + 1; j < row; ++j) {
                a[i*lda + j] /= a[i*lda + i];
                for (int k = i + 1; k < p + colToGo; k++) {
                    a[k*lda + j] -= a[i*lda + j] * a[k*lda + i];
                }
            }
        }
        ////END DGETRF2
        //clock_gettime(CLOCK_MONOTONIC, &cp1);
        //Pivot indices are correct, no need to do a correction

        //Apply interchanges to columns 1:p-1

        //dlawsp
        //           K1      K2
        for (int i = p; i < p + pb; ++i) {
            if (ipiv[i] != i + 1) {
                //swap rows
                //                 col of A
                for (int j = 0; j < p; ++j) {
                    double tmp = a[j*lda + i];
                    a[j*lda + i] = a[j*lda + ipiv[i] - 1];
                    a[j*lda + ipiv[i] - 1] = tmp;
                }
            }
        }
        //clock_gettime(CLOCK_MONOTONIC, &cp2);
        bool inthere = false;
        //Line 197
        if (p + pb < col) {
            inthere = true;
            //Apply interchanges to columns p+pb:n
            //dlawsp
            //          K1        K2
            for (int i = p; i < p + pb; ++i) {
                if (ipiv[i] != i + 1) {
                    //swap rows
                    for (int j = p + pb; j < col; ++j) {
                        double tmp = a[j*lda + i];
                        a[j*lda + i] = a[j*lda + ipiv[i] - 1];
                        a[j*lda + ipiv[i] - 1] = tmp;
                    }
                }
            }
            //clock_gettime(CLOCK_MONOTONIC, &cp3);
            //Compute block row of U
            //Lower triangular dtrsm
            for (int j = p + pb; j < col; ++j) { //col of X b
                for (int i = p; i < p + pb; ++i) {  //row of X b
                    for (int k = p; k < i; ++k) {
                        a[j*lda + i] -= a[j*lda + k] * a[k*lda + i];
                    }
                }
            }
            //clock_gettime(CLOCK_MONOTONIC, &cp4);
            if (p + pb < row) {
                //Update trailing submatrix
                //DEGMM
                //BLOCKED MM

                int B = 64;
                for (int j = p + pb; j < col; j += B) {
                    int j2 = std::min(j + B, col);
                    for (int k = p; k < p + pb; k += B) {
                        int k2 = std::min(k + B, p + pb);
                        for (int i = p + pb; i < row; i += B) {
                            int i2 = std::min(i + B, row);
                            for (int j1 = j; j1 < j2; ++j1)
                                for (int k1 = k; k1 < k2; ++k1)
                                    for (int i1 = i; i1 < i2; ++i1)
                                        //a[i][j] -= a[i][k]*a[k][j]
                                        a[j1*lda + i1] -= a[k1*lda + i1] * a[j1*lda + k1];
                        }
                    }
                }



            }
        }//if (p + pb < col)
    }


    return 0;
}



//My dtrsm
//signiture based on LAPACKE_dgetrs
int mydtrsm(char trans, int n, int nrhs, double *a, int lda, int* ipiv, double *b, int ldb) {
    if (trans != 'N') {
        std::cout << "ERROR, ONLY ACCEPT N TYPE MATRIX" << std::endl;
        return -1;
    }
    if ((nrhs != 1) || (lda != ldb) || (lda != n)) {
        std::cout << "ERROR, NOT SUPPORTED." << std::endl;
        return -1;
    }

    //Forward Substitution
    //Preprocess
    for (int i = 0; i < n; ++i) {
        double temp = b[ipiv[i] - 1];
        b[ipiv[i] - 1] = b[i];
        b[i] = temp;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            b[i] -= b[j] * a[j*n + i];
        }
    }

    //Backward Substitution
    for (int i = n - 1; i >= 0; --i) {
        for (int j = i + 1; j < n; j++) {
            b[i] -= b[j] * a[j*n + i];
        }
        b[i] /= a[i*n + i];
    }
    return 0;
}


int main(int argc, char *argv[]) {
    double *a, *b;
    int n = 3, bs =1, opt;
/*
    while ((opt = getopt(argc, argv, "n:b:")) != EOF) {
        switch (opt) {
        case 'n':
            n = atoi(optarg);
            break;
        case 'b':
            bs = atoi(optarg);
            break;
        case '?':
        default:
            std::cerr << "Usage run -n <size> [-b <block size>]" << std::endl;
            return -1;
        }
    }
    if (n == 0) {
        std::cerr << "Usage run -n <size> [-b <block size>]" << std::endl;
        return -1;
    }
*/
    srand(419);

    a = new double[n*n];
    b = new double[n];
    
  a[0] = 4; a[1] = 6; a[2] = 3; b[0] = 2;                                                               
  a[3] = 3; a[4] = 3; a[5] = 6; b[1] = 10; 
  a[6] = 1; a[7] = 4; a[8] = 1; b[2] = 2;
/*
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a[i*n + j] = (double)((rand() << 15) | rand()) / (double)rand();
        }
        b[i] = (double)((rand() << 15) | rand()) / (double)rand();
    }
  */
  


 
    double *al, *bl, *ag, *bg;
    al = new double[n*n];
    bl = new double[n];
    ag = new double[n*n];
    bg = new double[n];
    memcpy(al, a, n*n * sizeof(double));
    memcpy(bl, b, n * sizeof(double));
    memcpy(ag, a, n*n * sizeof(double));
    memcpy(bg, b, n * sizeof(double));

    testLapack(al, bl, n);
    
    //printing out a
    printf("\n Lapack a: \n");
    for (int i = 0; i < n; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", al[j+i]);                                                        
      }
                                                                                                
      printf("\n");          
    }   
  


    testMine(a, b, n);

    
    //printing out a
    printf("\n My  a: \n");
    for (int i = 0; i < n; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", a[j+i]);                                                        
      }
                                                                                            
      printf("\n");          
    }   
  
    testBlcoked(ag, bg, n, bs);
    double sumOfSquare = 0;
    double norm;
    
    for (int i = 0; i < n; ++i) {
        sumOfSquare += (b[i] - bl[i])*(b[i] - bl[i]);
    }
    norm = sqrt(sumOfSquare);
    std::cout << "The norm of difference between LAPACK and My unoptimized algorithm is " << std::scientific << norm << std::endl;
    
    sumOfSquare = 0;
    for (int i = 0; i < n; ++i) {
        sumOfSquare += (bg[i] - bl[i])*(bg[i] - bl[i]);
    }
    norm = sqrt(sumOfSquare);
    std::cout << "The norm of difference between LAPACK and My OPTIMIZED algorithm is " << std::scientific << norm << std::endl;
    
     
    //printing out a
    printf("\n altered a: \n");
    for (int i = 0; i < n; i++) {                                                                
      for (int j = 0; j < n; j++){                                                             
          printf("  %lf ", a[j+i]);                                                        
      }                                                                                        
    }   
      printf("\n");          
  

    //printing out b
    printf("\naltered b:\n");
    for(int i=0; i < n; i++){
     printf("  %lf ", b[i]);
     printf("\n");

    }    


    delete[] a;
    delete[] b;
    delete[] al;
    delete[] bl;
    delete[] ag;
    delete[] bg;

    return 0;
}

