#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#define GG 20
#define DG 0
#include<math.h>
void input_matrix( FILE *pf, double **M, int rows,int colums)
{
  int i,j;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < colums; j++) {
    	#pragma omp ordered 
      fscanf(pf, "%lf", &M[i][j]);
    }
  }  
}

void print_matrix(double **M, int rows,int colums) 
{
  int i, j;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < colums; j++) {
    	#pragma omp ordered
      printf("%lf\t", M[i][j]);
      
    }
    printf("\n");
  }   
}

double** allocate_matrix(int rows,int colums) 
{  
  double **M;
  int i;
  
  M = (double**) malloc ( rows* sizeof(double*) );
  for (i = 0; i < rows; i++) {
    M[i] = (double*) malloc ( colums * sizeof(double));
  }
  return M;
}

void make_diag_matrix(double **M,int n)
{
int i,j;

for(i=0;i<n;i++){
	for(j=0;j<n;j++){
		if(i==j){
		M[i][j]=1;
		}
		else M[i][j]=0;
	}
}	
}

void gem(double **A,double **L,int n,int k)
{
    /*Gaussian elimination
    * Computation of coefficient ratio between rows stored in matrix L
    * Reduction of rows starting from row=k+1,column=k in matrix A
    */
    int i,j;
    #pragma omp parallel for shared(L,A) private(j,i)
    for(i=k+1;i<=n-1;i++){
              L[i][k]=A[i][k]/A[k][k]; 
              for(j=k;j<=n-1;j++){
                  A[i][j]=A[i][j]-L[i][k]*A[k][j];
              }
            }
}

void swap(double **A,int k,int l,int index)
{
    double temp;
	temp=A[k][l];
    A[k][l]=A[index][l];
    A[index][l]=temp;
	
}
void pivoting(double **A,double **P,double **L,int n)
{
    /*Pivoting and Gaussian elimination
    *Searching for j>k such that |a(j)(k)|>|a(k)(k)|
    *If j exists, switch A(j)(:)<->A(k)(:), P(j)(:)<->P(k)(:) and L(j)(0:k-1)<->L(k)(0:k-1)
    *Then we call the function gem that does the gauss elimination on column k
    */
    int j,k,index,l;
    double max;
    for(k=0;k<n-1;k++){
         max= fabs(A[k][k]);
         index=k;
	    #pragma omp parallel for shared(A,k) private(j) reduction(max:max,index)
        for(j=k+1;j<=n-1;j++){
            if(fabs(A[j][k])>max){
                max=fabs(A[j][k]);
                index=j;
            }
           
           
        }
        if(index!=k){
            #pragma omp parallel for shared(A,P,L,n,k,index) schedule(static,2)
            for( l=0;l<=n-1;l++){
            	swap(A,k,l,index);
                swap(P,k,l,index);
                if(l<k){
                	swap(L,k,l,index);
                }
            }
		 }
		 gem(A,L,n,k);
	}
}
void forwardsub(double **P,double **L,int n)
{
    /*Algorithm for solving linear sistem that involves a lower triangular matrix-> Ly=P
    *The computation starts from the row 0 and uses L(0)(:),P(0)(:) to obtain y(0)(:)
    *Iterative variable k goes from 0 to n-2
    *Then iteration k uses L(k)(:),P(k)(:),y(0:k-1) to obtain y(k)(:)
    */
	int i,j,k;
	//#pragma omp parallel for shared(L,P,n) schedule(static,2) private(i,j,k)
    //for(j=0;j<n;j++){
	 for(k=0;k<=n-2;k++){
            #pragma omp parallel for shared(L,P,n,k) schedule(static,2) private(i,j)
	 	     for(i=k+1;i<=n-1;i++){
                 //printf("Thread %d  iteration %d\n", omp_get_thread_num(),k);
                for(j=0;j<n;j++){	             
                    P[i][j]-=L[i][k]*P[k][j];
                }
            }
        }
}
void backwardsub(double **A,double **P,int n)
{
    /*Algorithm for solving linear sistem that involves an upper triangular matrix-> Ux=P
    *The computation starts from the row n and uses L(n)(:),P(n)(:) to obtain y(n)(:)
    *Iterative variable k goes from n-1 to 0
    *Then iteration k uses L(k)(:),P(k)(:),y(k+1:n) to obtain y(k)(:)
    */
	int i,j,k;
	#pragma omp parallel for shared(A,P,n) schedule(static,2) private(j,k,i)
	for(j=0;j<n;j++){
       P[n-1][j]/=A[n-1][n-1]; 
       //printf("Thread %d  iteration %d\n", omp_get_thread_num(),j);
            for(k=n-2;k>=0;k--){
                for(i=k+1;i<=n-1;i++){
                    P[k][j]-=A[k][i]*P[i][j];
                   // printf("P[%d]=P[%d]\n", k,i);
                }
                P[k][j]/=A[k][k];  
            }
    }
}
void release_mem(double **M,int rows)
{
int i;

#pragma omp parallel for if(rows*rows>7000) private(i)
for(i=0;i<rows;i++)
free(M[i]);
free(M);
}

void random(double **M,int rows,int columns)
{
	int i,j;
	
	//srand ((unsigned)time(NULL ));
	for(i=0;i<rows;i++){
		for(j=0;j<columns;j++){
			M[i][j]=DG+ (float)rand () / RAND_MAX * (GG - DG);

		}
	}
	
}
void process(int rows)
{
  double **L, **P,**A;
  
  int i, cols;
  cols=rows;
  A=allocate_matrix(rows,cols);
  P=allocate_matrix(rows,cols);
  L=allocate_matrix(rows,cols);
  random(A,rows,cols);
  make_diag_matrix(P,rows);
  make_diag_matrix(L,rows);
  //printf("\n%d. matrix\n\n",i);
  //printf("P:\n");
  //print_matrix(P,rows,cols);
  //printf("A:\n");
  //print_matrix(A,rows,cols);
  //printf("After pivoting:A\n");
  clock_t t; 
  t = clock(); 
  pivoting(A,P,L,rows);
  //print_matrix(A,rows,cols);
  //printf("After pivoting:P\n");
  //print_matrix(P,rows,cols);
  //printf("After GEM :A:\n");
  //print_matrix(A,rows,cols);
  //printf("After GEM :L:\n");
  //print_matrix(L,rows,cols);
  forwardsub(P,L,rows);
  //printf("After GEM:P\n");
  //print_matrix(P,rows,cols);
  backwardsub(A,P,rows);
  t = clock() - t;
  double time_taken = ((double)t)/CLOCKS_PER_SEC; 
  //printf("Inverse of A:\n");
 //print_matrix(P,rows,cols);
  printf("Time elapsed: %lf\n",time_taken);
  release_mem(P,rows);
  release_mem(A,rows);
  release_mem(L,rows);
}

int main(void)
{
    for(int dim=1000;dim<1001;dim+=10){
        printf("Dimension: %d\n",dim);
            process(dim);
    }
}



