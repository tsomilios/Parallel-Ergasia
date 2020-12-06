#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include "mmio.h"
#include "coo2csc.h"

#include <omp.h>

#define CHUNKSIZE   100

void print1DMatrix(int* matrix, int size){
    int i = 0;
    for(i = 0; i < size; i++){
        printf("%d: %d \n",i, matrix[i]);
    }
}

int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    uint32_t M, N, nz;   
    int i, *I, *J;
    double *val;
    int binary = atoi(argv[2]);
    int num_of_threads = atoi(argv[3]);
    struct timeval start, end;

    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename] [0 for non binary 1 for binary matrix]\n", argv[0]);
		exit(1);
	}
    else    
    { 
        if ((f = fopen(argv[1], "r")) == NULL) 
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    /* reseve memory for matrices */
    /* For the COO */
    I = (uint32_t *) malloc(nz * sizeof(uint32_t));
    J = (uint32_t *) malloc(nz * sizeof(uint32_t));
    val = (double *) malloc(nz * sizeof(double));

    /* For the CSC */
    uint32_t* cscRow = (uint32_t *) malloc(nz * sizeof(uint32_t));
    uint32_t* cscColumn = (uint32_t *) malloc((N + 1) * sizeof(uint32_t));

    /* Depending on the second argument of the main call our original matrix may be binary or non binary so we read the file accordingly */
   
    if (!mm_is_pattern(matcode))
    {
    for (uint32_t i=0; i<nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }
    }
    else
    {
    for (uint32_t i=0; i<nz; i++)
    {
        fscanf(f, "%d %d\n", &I[i], &J[i]);
        val[i]=1;
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }
    }




    if (f !=stdin) fclose(f);

    if(M != N) {
        printf("COO matrix' columns and rows are not the same");
    }

    /*
   
    Code that converts any symmetric matrix in upper/lower triangular
    */
    int flag = 0;
    if(I[0] > J[0]) {
        flag = 1;
    }
    if(flag == 0){
        printf("Ypper trianglular I,J \n");
        coo2csc(cscRow, cscColumn, I, J, nz, M, 0);
    }
    else if(flag == 1){
        printf("Lower triangle J,L \n");
        coo2csc(cscRow, cscColumn, J, I, nz, N, 0);
    }
    else{
        exit;
    } 

    /* Initialize c3 with zeros*/
    int* c3;
    c3 = malloc(N * sizeof c3);    
    for(int i = 0; i < N; i++){
        c3[i] = 0;
    }

    printf("Matrix Loaded, now Searching!\n");

   

    int sum = 0;
    omp_set_dynamic(0);     // Disabling dynamic teams
    omp_set_num_threads(num_of_threads); // Use the same number of threads for all parallel regions
    
     /* We measure time from this point */
    gettimeofday(&start,NULL);
    
    #pragma omp parallel for shared(sum, c3)
    
    for(int i = 1; i < N; i++) {
        for(int j = 0; j < cscColumn[i+1] - cscColumn[i]; j++) {
            int row1 = cscRow[cscColumn[i] + j];
            int col1 = i;
            for(int k = 0; k < cscColumn[row1+1] - cscColumn[row1]; k++) {
                int row2 = cscRow[cscColumn[row1] + k];
                int col2 = row1;                
                if(row2>col1) {
                    for (int l = 0; l < cscColumn[row2+1] -cscColumn[row2]; l++) {
                        int temp = cscRow[cscColumn[row2] + l];
                        if(temp == col1) {
                            #pragma omp critical
                            sum++;
                            c3[col1]++;
                            c3[row2]++;
                            c3[col2]++;
                        }
                    }
                }
                else {
                    for (int l = 0; l < cscColumn[col1+1] - cscColumn[col1]; l++) {
                        int temp = cscRow[cscColumn[col1] + l];
                        if(temp == row2) {
                            #pragma omp critical
                            sum++;
                            c3[col1]++;
                            c3[row2]++;
                            c3[col2]++;
                        }
                    }
                }
            }
        }
    }

    /* We stop measuring time at this point */
    gettimeofday(&end,NULL);
    double duration = (end.tv_sec+(double)end.tv_usec/1000000) - (start.tv_sec+(double)start.tv_usec/1000000);

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    /*for (i=0; i<nz; i++){
        fprintf(stdout, "%d %d %20.19g\n", I[i]+1, J[i]+1, val[i]);
    }*/
    printf("Threads: %d \n",num_of_threads);
    printf("Sum: %d \n", sum);
    printf("Duration: %f \n", duration);



    /* Deallocate the arrays */
    free(I);
    free(J);
    free(c3);

	return 0;
}

