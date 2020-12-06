#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include "mmio.h"
#include "coo2csc.h"

#include <pthread.h>

#define MAX_THREAD 1000

 struct matrix{
    uint32_t* cscRow;
    uint32_t* cscColumn;
    uint32_t* c_values;
    int nz;
    int start;
    int end;
    int id;
 };

void *multiplication(void* arg) {
    struct matrix* mul_matrix = arg; 

    for(int i = mul_matrix->start; i < mul_matrix->end; i++) {
      for(int j = 0; j < mul_matrix->cscColumn[i+1] - mul_matrix->cscColumn[i]; j++) {
          int index_p = mul_matrix->cscRow[mul_matrix->cscColumn[i] + j];
          int index_col = i;
          int k_size = mul_matrix->cscColumn[index_p+1] - mul_matrix->cscColumn[index_p];  
          int l_size = mul_matrix->cscColumn[index_col+1] - mul_matrix->cscColumn[index_col];    
          int *indice_l = malloc((l_size) * sizeof(int));
          int *indice_k = malloc((k_size) * sizeof(int));
          for(int x = 0; x < k_size; x++) {
              indice_k[x] = mul_matrix->cscRow[mul_matrix->cscColumn[index_p] + x];
          }
          for(int x = 0; x < l_size; x++) {
              indice_l[x] = mul_matrix->cscRow[mul_matrix->cscColumn[index_col] + x];
          }
          int k_pointer = 0;
          int l_pointer = 0;
          int value = 0;
          while(k_pointer != k_size && l_pointer != l_size) {
              if(indice_k[k_pointer] == indice_l[l_pointer]) {
                  value++;
                  k_pointer++;
                  l_pointer++;
              }
              else if(indice_k[k_pointer] > indice_l[l_pointer]) {
                  l_pointer++;
              }
              else
              {
                  k_pointer++;
              }                
          }        
          if(value) {
              mul_matrix->c_values[mul_matrix->cscColumn[i] + j] = value;
          }
      }
    } 

    pthread_exit(NULL);
}

int main(int argc, char* argv[]) {
  
  int ret_code;
    MM_typecode matcode;
    FILE *f;
    uint32_t M, N, nz;   
    int i,*I, *J;
    double *val;
    int binary = atoi(argv[2]);
    int num_of_threads = atoi(argv[3]);
    struct timeval start, end;

    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename] [0 for binary or 1 for non binary] [num of threads]\n", argv[0]);
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

    /* Reseve memory for matrices */

    /* For the COO */
    /* Double the memory to store the full matrix */
    I = (uint32_t *) malloc(2 * nz * sizeof(uint32_t));
    J = (uint32_t *) malloc(2 * nz * sizeof(uint32_t));
    val = (double *) malloc(nz * sizeof(double));

    
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
    
  
    
    /* Stop reading coo file */
    if (f !=stdin) fclose(f);

    if(M != N) {
        printf("COO matrix' columns and rows are not the same");
    }

    /* Add each value once more so we get the full symmetric matrix */
    /* Requires every element in the diagonal to be zero */
    for(uint32_t i = 0; i < nz; i++) {
        I[nz + i] = J[i];
        J[nz + i] = I[i];
    }

    /* For the CSC */
    uint32_t* cscRow = (uint32_t *) malloc(2 * nz * sizeof(uint32_t));
    uint32_t* cscColumn = (uint32_t *) malloc((N + 1) * sizeof(uint32_t));

   
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

    /* For the C CSC */
    uint32_t* c_cscRow = (uint32_t *) malloc(2 * nz * sizeof(uint32_t));
    uint32_t* c_values = (uint32_t *) malloc(2 * nz * sizeof(uint32_t));
    uint32_t* c_cscColumn = (uint32_t *) malloc((N + 1) * sizeof(uint32_t));

    printf("Matrix Loaded, now Searching!\n");

    /* Initialize c3 with zeros*/
    int* c3;
    c3 = malloc(N * sizeof c3);    
    for(int i = 0; i < N; i++){
        c3[i] = 0;
    }

    /* Initialize t with ones*/
    int* t;
    t = malloc(N * sizeof t);    
    for(int i = 0; i < N; i++){
        t[i] = 1;
    }

    /* Initialize result with zeros*/
    int* result_vector;
    result_vector = malloc(N * sizeof result_vector);    
    for(int i = 0; i < N; i++){
        result_vector[i] = 0;
    }

    /* We measure time from this point */
    //gettimeofday(&start,NULL);

    //Assign matrix atributes
    struct matrix matrix[num_of_threads];

    pthread_t *threads;
    threads = (pthread_t *)malloc(num_of_threads*sizeof(pthread_t));

    /* We measure time from this point */
    gettimeofday(&start,NULL); 

    //Parallelize the for loop by breaking it into chunks
    int chunk = 1;
    if(num_of_threads > 0) {
        chunk = N / (num_of_threads);
    }

 


    for(int i = 0; i < num_of_threads-1; i++) {
      matrix[i].cscRow = cscRow;
      matrix[i].cscColumn = cscColumn;
      matrix[i].c_values = c_values;
      matrix[i].start = i * chunk;
      matrix[i].end = matrix[i].start + chunk;
      matrix[i].nz = nz;
      matrix[i].id = i;
      //multiplication(&matrix[i]);
      pthread_create(&threads[i], NULL, multiplication, &matrix[i]);
      
    }
    // The last thread is left out so as to calculate the mod of the chunk division!
    matrix[num_of_threads - 1].cscRow = cscRow;
    matrix[num_of_threads - 1].cscColumn = cscColumn;
    matrix[num_of_threads - 1].c_values = c_values;
    matrix[num_of_threads - 1].start = (num_of_threads - 1) * chunk;
    matrix[num_of_threads - 1].end = matrix[num_of_threads - 1].start + chunk + (N % num_of_threads);
    matrix[num_of_threads - 1].nz = nz;
    matrix[num_of_threads - 1].id = num_of_threads - 1;    
    pthread_create(&threads[num_of_threads - 1], NULL, multiplication, &matrix[num_of_threads - 1]);

    for(int i = 0; i < num_of_threads; i++) {
      pthread_join(threads[i], NULL);
    }

    c_cscColumn = cscColumn;
    c_cscRow = cscRow;


  
    /* Multiplication of a NxN matrix with a Nx1 vector*/
    for(int i = 0; i < N; i++) {
        // printf("i: %d \n", i);
        for(int j = 0; j < c_cscColumn[i+1] - c_cscColumn[i]; j++) {
            int row = c_cscRow[c_cscColumn[i] + j];
            int col = i;
            int value = c_values[c_cscColumn[i] + j];
            result_vector[row] += value * t[col];
        }
    }
    
    int triangle_sum = 0;
    for(int i = 0; i < N; i++) {
        c3[i] = result_vector[i] / 2;
        triangle_sum += c3[i];
    }

    triangle_sum = triangle_sum / 3;

    /* We stop measuring time at this point */
    gettimeofday(&end,NULL);
    double duration = (end.tv_sec+(double)end.tv_usec/1000000) - (start.tv_sec+(double)start.tv_usec/1000000);
    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    //for (i=0; i<nz; i++){
        //fprintf(stdout, "%d %d %20.19g\n", I[i]+1, J[i]+1, val[i]);
    //}
    printf("\nNum p threads: %d",  num_of_threads);
    //printf("\nTriangle Sum: %d",  triangle_sum);
    printf("\nDuration: %f\n",  duration);
  
    free(I);
    free(J);
    free(c3);
    free(t);
    free(result_vector);

	return 0;
}
