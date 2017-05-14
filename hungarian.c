/* The following algorithm is an implementation of the Munkres' Assignment 
 * Algorithm (sometimes referred to as the Hungarian Algorithm). Munkres' 
 * algorithm is O(n^3) and is a modification from the original Kuhn's algorithm.
 * This implementation is for rectangular matrix (nxm) where n <= m. */

#include "global.h"
#define eps 0.00001

double **C;     /*cost matrix*/
int **M, **path;
int *RowCover, *ColCover, *Selection;
int nrow;       /* number of weight vectors */
int ncol;       /* number of solution vectors */
int nfunc;      /* number of functions */
int path_count;
int path_row_0;
int path_col_0;
int step;

/* Main function to apply the Kuhn-Munkres' algorithm. The T matrix have the 
 * form nxk and contain n solution vectors with k objective functions. The W 
 * matrix have the form kxm and contain m weight vectors uniformly scattered
 * in the k-dimensional objective space. 
 * Return who solution was selected and was not. */
int* hungarian(double **T, double **W)
{
  InitMunkres(T, W);
  step = 1;
  RunMunkres();  
  return Selection;  
}

/*            *** Construct the linear assignment problem ***          */

/* Create an nxm  matrix called the cost matrix in which each element 
 * represents the cost of assigning one of n solutions to one of m weights.  
 * Rotate the matrix so that there are at least as many columns as rows.*/
void InitMunkres(double **T, double **W)
{    
    path_count = 0;
    int r,c;
    for (r = 0; r < nrow; r++)
    {
        RowCover[r] = 0;        
        for (c = 0; c < ncol; c++)
        {
            C[r][c] = 0.0;
            M[r][c] = 0;
        }
    }
    for (c = 0; c < ncol; c++)
    {
        ColCover[c]  = 0;
        Selection[c] = 0;
    }
    
    for (r = 0; r < 2*ncol+1; r++)
        for (c = 0; c < 2; c++)
            path[r][c] = 0;
    
    Normalize(T);    
    CostMatrix(T, W);
}

/* Normalize the objective values of all current solutions */
void Normalize(double **T)
{    
    int r, c;
    double maxval, minval;
    for (c = 0; c < nfunc; c++)
    {
        maxval = T[0][c];
        minval = T[0][c];
        for (r = 1; r < ncol; r++)
        {            
            if (maxval < T[r][c])
                maxval = T[r][c];
            if (minval > T[r][c])
                minval = T[r][c];
        }
        if (minval == maxval)
            minval = 0.0;
        if (maxval == 0.0)
            maxval = 1.0;        
        for (r = 0; r < ncol; r++)
            T[r][c] = (T[r][c] - minval)/(maxval - minval);
    }
}

void CostMatrix(double **T, double **W)
{
    /*the cost matrix C is the transpose of Utility(Tk|Wk), k is the objective*/        
    int j, i, k;
    for (i = 0; i < ncol; i++) 
       for (j = 0; j < nrow; j++)           
           for (k = 0; k < nfunc; k++)  /* modified Tchebycheff decomposition */
                if (C[j][i] < T[i][k]/W[k][j])
                    C[j][i] = T[i][k]/W[k][j];
}

/*              *** Kuhn-Munkres' algorithm implementation ***           */

void RunMunkres() 
{
    int done = 0;
    while (done == 0) 
    {
        switch (step) {
            case 1:
                step_one();
                break;
            case 2:
                step_two();
                break;
            case 3:
                step_three();
                break;
            case 4:
                step_four();
                break;
            case 5:
                step_five();
                break;
            case 6:
                step_six();
                break;
            case 7:
                step_seven();
                done = 1;
                break;
        }
    }    
}
               /*** STEP 1 ***/
/* For each row of the matrix, find the smallest element and subtract it from 
 * every element in its row.  Go to Step 2.*/
void step_one() 
{
    double min_in_row;
    int r,c;
    for (r = 0; r < nrow; r++) 
    {
        min_in_row = C[r][0];
        for (c = 1; c < ncol; c++)
            if (C[r][c] < min_in_row)
                min_in_row = C[r][c];
        for (c = 0; c < ncol; c++)
            C[r][c] -= min_in_row;
    }
    step = 2;
}
               /*** STEP 2 ***/
/* Find a zero (Z) in the resulting matrix. If there is no starred zero in 
 * its row or column, star Z. Repeat for each element in the matrix. 
 * Go to Step 3.*/
void step_two() 
{
    int r,c;
    for (r = 0; r < nrow; r++)
        for (c = 0; c < ncol; c++) 
            if (C[r][c] < eps && ColCover[c] == 0)
            {
                M[r][c] = 1;
                ColCover[c] = 1;
                break;
            }
    step = 3;
}
               /*** STEP 3 ***/
/* Cover each column containing a starred zero. If K columns are covered, the 
 * starred zeros describe a complete set of unique assignments. In this case, 
 * go to Step 7, otherwise, go to Step 4*/
void step_three() 
{
    int colcount = 0, r, c;
    for (c = 0; c < ncol; c++)
        for (r = 0; r < nrow; r++)
            if (M[r][c] == 1)
            {
                ColCover[c] = 1;
                colcount ++;
                break;
            }    
    if (colcount >= ncol || colcount >= nrow)
        step = 7;
    else
        step = 4;
}
               /*** STEP 4 ***/
/* Find an uncovered zero and prime it. If there is no starred zero in the row
 * containing this primed zero, go to Step 5. Otherwise, cover this row and 
 * uncover the column containing the starred zero. Repeat this process 
 * until there are no uncovered zeros left and go to Step 6.*/
void find_a_zero(int *row, int *col)
{
    int r,c;
    *row = -1;
    *col = -1;
    for (r = 0; r < nrow; r++)
        if (RowCover[r] == 0)
            for (c = 0; c < ncol; c++) 
                if (C[r][c] < eps && ColCover[c] == 0)
                {
                    *row = r;
                    *col = c;
                    return;
                }
    return;        
}

int find_star_in_row(int row) 
{
    int col = -1, c;
    for (c = 0; c < ncol; c++)
        if (M[row][c] == 1)
        {
            col = c;
            break;
        }
    return col;
}

void step_four() 
{
    int row = -1;
    int col = -1;
    int done = 0, c;
    while (done == 0) 
    {
        find_a_zero(&row, &col);
        if (row == -1) 
        {
            done = 1;
            step = 6;
        }
        else 
        {
            M[row][col] = 2;
            c = find_star_in_row(row);
            if (c != -1) 
            {
                RowCover[row] = 1;
                ColCover[c] = 0;
            }
            else 
            {
                done = 1;
                step = 5;
                path_row_0 = row;
                path_col_0 = col;
            }
        }
    }
}
               /*** STEP 5 ***/
/* Construct a serie of alternating primed and starred zeros as follows.  
 * Let Z0 represent the uncovered primed zero found in Step 4. Let Z1 denote 
 * the starred zero in the column of Z0 (if any). Let Z2 denote the primed zero 
 * in the row of Z1 (there will always be one). Continue until the serie 
 * terminates at a primed zero that has no starred zero in its column. Unstar 
 * each starred zero of the serie, star each primed zero of the serie, erase 
 * all primes and uncover every line in the matrix.  Return to Step 3.*/
void find_star_in_col(int c, int* r) 
{
    *r = -1;
    int row;
    for (row = 0; row < nrow; row++)
        if (M[row][c] == 1)
        {
            *r = row;
            break;
        }
}

void find_prime_in_row(int r, int* c) 
{
    int col;
    for (col = 0; col < ncol; col++)
        if (M[r][col] == 2)
        {
            *c = col;
            break;
        }
}

void augment_path() 
{
    int p;
    for (p = 0; p < path_count; p++)
        if (M[path[p][0]][path[p][1]] == 1)
            M[path[p][0]][path[p][1]] = 0;
        else
            M[path[p][0]][path[p][1]] = 1;
}

void clear_covers() 
{
    int r,c;
    for (r = 0; r < nrow; r++)
        RowCover[r] = 0;
    for (c = 0; c < ncol; c++)
        ColCover[c] = 0;
}

void erase_primes() 
{
    int r,c;
    for (r = 0; r < nrow; r++)
        for (c = 0; c < ncol; c++)
            if (M[r][c] == 2)
                M[r][c] = 0;
}

void step_five() 
{
    int done = 0;
    int r = -1;
    int c = -1;
    path_count = 1;
    path[path_count - 1][0] = path_row_0;
    path[path_count - 1][1] = path_col_0;
    while (done == 0) 
    {
        find_star_in_col(path[path_count - 1][1], &r);
        if (r > -1) 
        {
            path_count += 1;
            path[path_count - 1][0] = r;
            path[path_count - 1][1] = path[path_count - 2][1];
        }
        else
            done = 1;
        if (done == 0) 
        {
            find_prime_in_row(path[path_count - 1][0], &c);
            path_count += 1;
            path[path_count - 1][0] = path[path_count - 2][0];
            path[path_count - 1][1] = c;
        }
    }
    augment_path();
    clear_covers();
    erase_primes();
    step = 3;
}
               /*** STEP 6 ***/
/* Add the smallest uncovered value to every element of each covered row, and 
 * subtract it from every element of each uncovered column. Return to Step 4 
 * without altering any stars, primes, or covered lines. */
double find_smallest() 
{
    double minval = 100000.0;
    int r,c;
    for (c = 0; c < ncol; c++)
        if (ColCover[c] == 0)
            for (r = 0; r < nrow; r++) 
                if (minval > C[r][c] && RowCover[r] == 0)
                    minval = C[r][c];
    return minval;
}

void step_six() 
{
    int r, c;
    double minval;
    minval = find_smallest();
    for (r = 0; r < nrow; r++)
        for (c = 0; c < ncol; c++) 
        {
            if (RowCover[r] == 1)
                C[r][c] += minval;
            if (ColCover[c] == 0)
                C[r][c] -= minval;
        }
    step = 4;
}
               /*** STEP 7 ***/
/* Assignment pairs are indicated by the positions of the starred zeros in the
 * cost matrix. If C(i,j) is a starred zero, then the element associated with 
 * row i is assigned to the element associated with column j.*/
void step_seven()
{
    int r,c;
    for (r = 0; r < nrow; r++)
        for (c = 0; c < ncol; c++)
            if(M[r][c] == 1)
                Selection[c] = 1;
}

void InitMemoryKM (int nsol, int nfun, int nweight)
{
    nrow  = nweight;
    ncol  = nsol;    
    nfunc = nfun;
    int r;
    C = (double **)malloc(nrow*sizeof(double*));   /* cost matrix */
    M = (int **)malloc(nrow*sizeof(int*));         /* mask matrix */
    for (r = 0; r < nrow; r++)
    {
        C[r] = (double *)malloc(ncol*sizeof(double));
        M[r] = (int*)malloc(ncol*sizeof(int));        
    }
    RowCover  = (int*)malloc(nrow*sizeof(int)); 
    ColCover  = (int*)malloc(ncol*sizeof(int));   
    Selection = (int*)malloc(ncol*sizeof(int)); 
    path = (int **)malloc((2*ncol+1)*sizeof(int*)); /* path matrix */
    for (r = 0; r < 2*ncol+1; r++)
        path[r] = (int*)malloc(2*sizeof(int));
}

void FreeMemoryKM()
{
    int r;    
    for (r = 0; r < nrow; r++)
    {
        free(C[r]);
        //free(M[r]);         
    }
    for (r = 0; r < 2*ncol+1; r++)
        free(path[r]);
    
    free(RowCover); 
    free(ColCover);
    free(Selection);
    free(path);
    free(C);
    free(M);
}