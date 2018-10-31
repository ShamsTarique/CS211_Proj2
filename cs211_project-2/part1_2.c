#include<stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include "cblas.h"
#include"lapacke.h"

using namespace std;

//TODO: take the code from the main() and make them into their own functions for mor eencapsulation -- estimated time to do: <10 minutes
//TODO: add user input and error checking so code does not break before uploading to Github <--good practice
//TODO: comment better before adding to GitHub 

int main(){
    
    int n = 3;
    int m =3;

    double* A = new double[m*n]; 
    double* B =new double[n];
    int* pivot = new int[n];
    double* Abk = new double[m*n];         
    
 A[0] = 4; A[1] = 6; A[2] = 3;  B[0] = 2;
 A[3] = 3; A[4] = 3; A[5] = 6;  B[1] = 10;                            
 A[6] = 1; A[7] = 4; A[8] = 1;  B[2] = 2;

 

//fill up A, B,Abk, and pivot
    for(int  i=0; i<n; i++){
       for(int j=0; j<n; j++){
        //  A[i*n+j] = Abk[i*n+j] = ((double)rand()*100)/((double)RAND_MAX + 0); //generate random double 0.00-100.00
       }
      // B[n] = ((double)rand()*100)/((double)RAND_MAX + 0);
       pivot[i] = i;
     }
    
printf("1");
    for(int i = 0; i<m*m; i++)
        Abk[i] = A[i]; 

//factorize --mydgetrf
printf("2");
    // pivot[n - 1] = n;
     for (int i = 0; i < n-1 ; ++i) {
         int maxp = i;
         double max = abs(A[i*n + i]);
         for (int t = i + 1; t < n; ++t) {
          //  if (abs(A[i*n + t]) > max) {
              if (abs(A[n*t + i]) > max) {
                 maxp = t;
                 max = abs(A[n*t + i]);
             }
         }
             //save pivoting infomation in LAPACK format
            // pivot[i] = maxp + 1;
            
            
              if (maxp != i) {
                //save pivot information
                int temps = pivot[i];
                pivot[i] = pivot[maxp];
                pivot[maxp] = temps;              
                //swap rows
                
                double* tempv = new double [m];
                for(int z = 0; z < m; z++){
                    tempv[z] = A[n*i +z];    
                    A[n*i+z] = A[maxp*n+z];
                    A[maxp*n+z] = tempv[z];
                }


//                 for (int j = 0; j < n; ++j) {
  //                   double tmp = A[j*n + i];
   //                  A[j*n + i] = A[j*n + maxp];
    //                 A[j*n + maxp] = tmp;
     //            }


printf("4");

             }
     
         for (int j = i + 1; j < n; ++j) {
             A[n*j + i] /= A[i*n + i];
             for (int k = i + 1; k < n; k++) {
                 A[n*j +k] -= A[n*j +i] * A[n*i +k];
             }
         }

    }

   printf("\nfactorization:\n");
   for (int i = 0; i < m; i++) {  
            for (int j = 0; j < n; j++){ 
                printf("  %lf ", A[i*n+j]);                 
       }                                                                
       printf("\n");
    }  


  printf("\nipiv contents:\n");
  for(int i = 0; i < m; i++)
    printf("  %d\n",pivot[i]); 
  printf("\n");
  



//TODO: fix both of these lower blocks V
//forward substitution --lower  --mydtrsm
      //make sure B matches the row swaps that happened in A  
     
    
    
printf("B before pivot swapping\n");
    for(int i = 0; i < n; i++)
        printf("%.15f \n", B[i]);  
 
      for(int i = 0; i <m; i++){
        double temp = 0;
        temp =  B[i];
        B[i] = B[pivot[i]];
        B[pivot[i]] = temp;
      }


printf("B after pivot swapping\n");
    for(int i = 0; i < n; i++)
        printf("%.15f \n", B[i]);  

    double *y = new double [m];
    for(int i =0;i<m; i++){
        y[i] =0;
    }    

    y[0] = B[pivot[0]];   
    for(int i = 1; i <n; i++){
        double sum = 0;
        for(int j = 1; j<i-1; j++)
            sum += y[j]*A[i*m+j];
        y[i] = B[pivot[i]]-sum;
        

    }

//back substitution -- upper --mydtrsm



double*x = new double [m];
for(int i =0; i <m; i++){
    x[i] = y[i]/A[m*i+i];
}

for( int i = n-1; i >0; i--){
    double sum = 0;
    for(int j = i; j<n; j++)
        sum += x[j]*A[m*i+j];
    x[i] = (y[i]-sum)/A[m*i+i];
}


printf("\nMyLU solution for x:\n");
for(int k =0; k<m; k++)
    printf("%.15f \n", x[k]);

 
return 0;
}
