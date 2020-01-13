/******************************************************************************

                            Online C Compiler.
                Code, Compile, Run and Debug C program online.
Write your code in this editor and press "Run" button to compile and execute it.

*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>



int main()
{
    int j,k,i,n=3;
    double A[3][3]={1,7,1,2,1,2,4,1,5},L[3][3]={1,0,0,0,1,0,0,0,1},P[3][3]={1,0,0,0,1,0,0,0,1};
    for(k=0;k<=n-2;k++){
        int max=A[k][k];
        int index=k;
        for(j=k+1;j<=n-1;j++){
            if(A[j][k]>max){
                max=A[j][k];
                index=j;
            }
        }
        if(index!=k){
            float temp;
            for(int l=0;l<=n-1;l++){
                temp=A[k][l];
                A[k][l]=A[index][l];
                A[index][l]=temp;
                temp=P[k][l];
                P[k][l]=P[index][l];
                P[index][l]=temp;
               
            }
        }
    }
    for(k=0;k<=n-2;k++){
                for(i=k+1;i<=n-1;i++){
          L[i][k]=A[i][k]/A[k][k];
          for(j=k;j<=n-1;j++){
              A[i][j]=A[i][j]-L[i][k]*A[k][j];
          }
        }
    }
    printf("U\n");
    for(k=0;k<=n-1;k++){
            printf("%f %f %f\n",A[k][0],A[k][1],A[k][2]);
    }
    printf("L\n");
    for(k=0;k<=n-1;k++){
            printf("%f %f %f\n",L[k][0],L[k][1],L[k][2]);
    }
        printf("P\n");
    for(k=0;k<=n-1;k++){
            printf("%f %f %f\n",P[k][0],P[k][1],P[k][2]);
    }
    
    
        for(k=0;k<=n-2;k++){
            for(i=k+1;i<=n-1;i++){
                for(j=0;j<n;j++){
                P[i][j]-=L[i][k]*P[k][j];
                }
            }
        }
    
    printf("P\n");
    for(k=0;k<=n-1;k++){
            printf("%f %f %f\n",P[k][0],P[k][1],P[k][2]);
    }
 
  for(j=0;j<n;j++){
       P[n-1][j]/=A[n-1][n-1];
      for(k=n-2;k>=0;k--){
            for(i=k+1;i<=n-1;i++){
                P[k][j]-=A[k][i]*P[i][j];
                }
                P[k][j]/=A[k][k];
            }
    }
        printf("P\n");
    for(k=0;k<=n-1;k++){
            printf("%f %f %f\n",P[k][0],P[k][1],P[k][2]);
    }
    
    
    
    
    
    
    
    
    
    
}
