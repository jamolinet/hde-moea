/* Generate a set of weight vectors using one of two methods. 
 * The first is the Simplex Lattice Design that obtain a discrete number of vectors.
 * The second is a Hammersley method join to Uniform Design of Mixture that 
 * obtain an arbitrary number of vectors well distributed. */

#include "global.h"

/************************
 * Simplex Lattice Points
 ************************/

/* calculates the number of weight vectors for a number of objective functions 
 * (nobj) and subdivision of each space (parts)*/
int NumberOfWeights()
{
    int n = parts + nobj - 1;
    int r = nobj - 1;
    int x = 1, y = 1, i;
    for (i = n; i > (n - r); i--)
        x *= i;
    for (i = r; i > 1; i--)
        y *= i;
    return x/y;
}

/* generate weight vectors */
void GetWeights()
{    
    int i, j, depth, current = 0;
    for (i = 0; i < parts+1; i++)
    {
        weights[0][current] = i*1.0/parts;        
        aux[0] = i;
        depth = 1;
        current = AuxWeights(current,depth);
    }
    for (i = 0; i < nobj; i++)
        for (j = 0; j < popsize; j++)
            if (weights[i][j] == 0.0)
                weights[i][j] = pow(10,-3);
}

int AuxWeights(int current, int depth)
{
    int i, j, sum = 0;
    for (i = 0; i < depth; i++)  
        sum += aux[i]; 
    
    if (depth == nobj - 1) 
    {
        weights[nobj-1][current] = 1 - sum * 1.0 / parts;
        current++;
    }
    else
        for (i = 0; i < parts - sum + 1; i++)
        {            
            weights[depth][current] = i*1.0/parts;
            aux[depth] = i;            
            current = AuxWeights(current, depth + 1);
            if (i < parts - sum + 1)
                for (j = 0; j < depth; j++)
                    weights[j][current] = aux[j]*1.0/parts;
        }
    return current;
}


/***************************************************
 * Uniform Design of Mixture with Hammersley Method 
 ***************************************************/

int primes[] = {2,3,5,7,11,13,17,19,23,29}; /*first 10 prime*/

/* Radical inverse of a number (index) in a prime base */
double RadicalInverse(int index, int base)
{
    double result = 0;
    double f = 1.0/base;
    int i = index;
    while(i > 0)
    {
        result += f*(i % base);
        i = (int)(i*1.0/base);
        f = f/base;
    }
    return result;
}

/* Hammersley Method to construct a low discrepancy set of n design points 
 * in a k-dimensional space */
void Hammersley(double **set, int n, int k)
{
    int i,j;
    for (i = 0; i < n; i++)
    {
        set[i][0] = (2*(i+1) - 1.0)/(2*n);
        for (j = 1; j < k; j++)
            set[i][j] = RadicalInverse(i+1, primes[j-1]);            
    }                   
}

/* Generate a set of weight vectors using Uniform Design of Experiments with Mixture. 
 * Transformation method of Wang and Fang (1990) */
void UDM()
{   
    int i,j,h;
    int n = popsize - nobj;
    int k = nobj;    
    double **set = (double**)malloc(n*sizeof(double*));
    double **w   = (double**)malloc(n*sizeof(double*));
    for (i = 0; i < n; i++)
    {
        set[i] = (double*)malloc((k-1)*sizeof(double));
        w[i]   = (double*)malloc(k*sizeof(double));
    }
    
    Hammersley(set, n, k-1);
    
    for (i = 0; i < n; i++)
        for (j = 0; j < k; j++)
        {
            if (j != k-1)
                w[i][j] = 1 - pow(set[i][j],1.0/(k-j-1));
            else
                w[i][j] = 1.0;
            for (h = 0; h < j; h++)
                w[i][j] *= pow(set[i][h],1.0/(k-h-1));
        }
        
    for (i = 0; i < nobj; i++)
    {   
        /* the first k weights vectors are extreme points */
        for (j = 0; j < nobj; j++)
            if (i == j)
                weights[i][j] = 1.0;
            else
                weights[i][j] = 0.001;  
        
        /* the rest of weights from the mixture */
        for (j = nobj; j < popsize; j++) 
            weights[i][j] = w[j-nobj][i];
    }
    
    for (i = 0; i < n; i++)
    {
        free(set[i]);
        free(w[i]);
    }
    free(set); free(w);
}