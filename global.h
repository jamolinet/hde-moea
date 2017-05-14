/* 
 * File:   global.h
 * Author: Jose Antonio Molinet Berenguer
 * Advisor: Carlos A. Coello Coello
 *
 * Created on 15 de septiembre de 2013, 10:55 PM
 */
#ifndef GLOBAL_H
#define	GLOBAL_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/unistd.h>
#include "rand.h"

# define GNUPLOT_COMMAND "gnuplot -persist"

#ifdef	__cplusplus
extern "C" {
#endif
    
    /* structures */
    typedef struct
    {
        double *vars;   /* real variables of the problem */
        double *objs;   /* objective functions */
    }
    individual;
    
    individual *oldpop; /* population of parents */
    individual *newpop; /* population of childs */
    double **weights;   /* weights vectors */
    double **solutions; /* solutions vectors */
    int *aux;           /* auxiliary structure */
    int *sel;           /* selection from old and new population */
    
    /* parameters */
    int popsize;        /* population size */
    int generations;    /* number of generations */
    int nvar;           /* number of variables */
    int nobj;           /* number of objective functions */
    double *upBound;    /* upper bound for each variable */
    double *lowBound;   /* lower bound for each variable */    
    double F;           /* mutation factor in DE */
    double CR;          /* crossover factor in DE */
    unsigned long seed; /* seed for the random generator */        
    int parts;          /* divisions of each dimension of weights vectors */
    void (*functions)(double*,double*,int,int);  /* objective function */
       
    /* Main functions */
    void InitMemory();
    void FreeMemory();
    double rnd(double, double);
    int rndint(int, int);
    void InitialPopulation();
    void GenerateOffspring();
    void Evaluate(individual*);
    void Select();
    void ReadParameters(char*);
    void SaveSolutions(char*);
    void (*getFunction(char *))(double*, double*,int,int);
    
    /* Generate weight vectors using Uniform Design of Mixture */
    double RadicalInverse(int, int);
    void Hammersley(double **, int, int);
    void UDM();
    
    /* Generate weight vectors using Simplex Lattice Points */
    int NumberOfWeights();
    int AuxWeights(int, int);
    void GetWeights();
    
    /* Functions to construct the assignment problem and 
     * to apply the Kuhn-Munkres' algorithm */
    void InitMemoryKM(int, int, int);
    void FreeMemoryKM();
    void Normalize(double **);
    void CostMatrix(double **, double **);
    void InitMunkres(double **, double **);
    void step_one();
    void step_two();
    void step_three();
    void find_a_zero(int *, int *);    
    int find_star_in_row(int );
    void step_four();
    void find_star_in_col(int , int *);
    void find_prime_in_row(int , int *);
    void augment_path();
    void clear_covers();
    void erase_primes();
    void step_five();
    double find_smallest();
    void step_six();
    void step_seven();
    void RunMunkres();
    int* hungarian(double **, double **);
    
     /* Functions to display the current population using GNUPLOT*/
    FILE *gp;
    int choice;         /* choose if the solutions will be plotted */
    void onthefly_display(int);

#ifdef	__cplusplus
}
#endif

#endif	/* GLOBAL_H */
