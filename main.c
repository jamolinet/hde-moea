/* 
 * File:   main.c
 * Author: Jose Antonio Molinet Berenguer
 * Advisor: Carlos A. Coello Coello
 *
 * Created on 15 de septiembre de 2013, 10:48 PM
 */

#include "global.h"

int main(int argc, char** argv) 
{
    if (argc < 4)
    {
        printf("Specify input and output files\n");
        exit(1);
    }
    
    int i; 
    char *input  = argv[1];
    char *output = argv[2];
    if (strcmp(argv[3],"1") == 0) 
        choice = 1;         /* the solutions will be plotted */
    else  
        choice = 0;
    seed = time(NULL);
    
    ReadParameters(input);    
    UDM();                  /* obtain weights vector using Uniform Design */   
    InitialPopulation();
    Evaluate(oldpop);    
                   
    for (i = 0; i < generations; i++)
    { 
        GenerateOffspring();
        Evaluate(newpop);
        Select(); 
        
        /* GNUPLOT */
        if (choice == 1)
            onthefly_display(i+1);  
    }       
    
    SaveSolutions(output);
    FreeMemory();      
    
    pclose(gp);  /* GNUPLOT */
    
    return (EXIT_SUCCESS);
}


/* Selection scheme using Kuhn-Munkres' algorithm */
void Select()
{
    int i,j,k;
    for (i = 0; i < popsize; i++)
        for (j = 0; j < nobj; j++)
        {
            solutions[i][j]         = oldpop[i].objs[j];
            solutions[i+popsize][j] = newpop[i].objs[j];
        }
    
    sel = hungarian(solutions, weights);
    
    k = 0;
    for (i = 0; i < popsize; i++)
        if (sel[i] == 0)   // individuals from oldpop to be replaced /
        {
            aux[k] = i;
            k++;
        }
    
    for (i = popsize; i < 2*popsize; i++)        
        if (sel[i] == 1)
        {
            k--;
            for (j = 0; j < nobj; j++)  
                oldpop[aux[k]].objs[j] = newpop[i-popsize].objs[j];
            for (j = 0; j < nvar; j++)
                oldpop[aux[k]].vars[j] = newpop[i-popsize].vars[j];     
        }  
}

/* Generate offspring using Differential Evolution algorithm */
void GenerateOffspring()
{
    int i, a, b, c, j, jrand;
    for (i = 0; i < popsize; i++) 
    {   
        /* select three different individuals */
        do { a = rndint(0,popsize-1); } while(a == i);
        do { b = rndint(0,popsize-1); } while(b == i || b == a);
        do { c = rndint(0,popsize-1); } while(c == i || c == a || c == b);
        
        /* apply crossover and mutation to obtain a new individual */
        jrand = rndint(0,nvar-1);
        for (j = 0; j < nvar; j++)
        {
            if (rnd(0,1) < CR || j == jrand)
            {
                newpop[i].vars[j] = oldpop[a].vars[j] + F*(oldpop[b].vars[j] - oldpop[c].vars[j]);
                
                if (newpop[i].vars[j] < lowBound[j])
                    newpop[i].vars[j] = lowBound[j];
                else 
                    if (newpop[i].vars[j] > upBound[j])
                        newpop[i].vars[j] = upBound[j];
            }
            else
                newpop[i].vars[j] = oldpop[i].vars[j];           
        }
    }
}

/* Evaluate each individual in the population */
void Evaluate(individual *pop)
{
    int i;
    for (i = 0; i < popsize; i++)
        (*functions)(pop[i].vars, pop[i].objs, nvar, nobj);
}


/* Create initial population with random values */
void InitialPopulation()
{
    int i,j;
    init_genrand(seed);
    for (i = 0; i < popsize; i++)
        for (j = 0; j < nvar; j++)
            oldpop[i].vars[j] = rnd(lowBound[j], upBound[j]);
}

/* Read the input file with the parameters */
void ReadParameters(char *input)
{
    int i;
    FILE *file;
    if ((file = fopen(input, "r")) == NULL)
    {
        printf("Error opening input file\n");
        exit(1);
    }
    char idfunc[20];
    fscanf(file, "%s", idfunc);    
    fscanf(file, "%d", &popsize);
    fscanf(file, "%d", &generations);
    fscanf(file, "%d", &nobj);
    fscanf(file, "%d", &nvar);
    
    InitMemory();
    
    for (i = 0; i < nvar; i++)
        fscanf(file, "%lf %lf", &lowBound[i], &upBound[i]);
    
    fscanf(file, "%lf", &F);
    fscanf(file, "%lf", &CR);
    fclose(file);
    
    functions = getFunction(idfunc); 
    
    /* GNUPLOT */
    gp = popen(GNUPLOT_COMMAND,"w");
    if (gp==NULL)
    {
        printf("\n Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND in file global.h\n");
        printf("\n Edit the string to suit your system configuration and re-run the program\n");
        exit(1);
    }
}

/* Store the final generation (Pareto Front and Pareto Set) */
void SaveSolutions(char *output)
{
    int i,j;
    FILE *file;
    if ((file = fopen(output, "w")) == NULL)
    {
        printf("Error opening output file\n");
        exit(1);
    }
    for (i = 0; i < popsize; i++)
    {
        for (j = 0; j < nobj; j++)
            fprintf(file, "%lf\t", oldpop[i].objs[j]);
        for (j = 0; j < nvar; j++)
            fprintf(file, "%lf\t", oldpop[i].vars[j]);
        fprintf(file, "\n");
    }
    fclose(file);
}

void InitMemory()
{
    int i;
    oldpop  = (individual*)malloc(popsize*sizeof(individual));
    newpop  = (individual*)malloc(popsize*sizeof(individual));    
    for (i = 0; i < popsize; i++)
    {
        oldpop[i].vars = (double*)malloc(nvar*sizeof(double));
        oldpop[i].objs = (double*)malloc(nobj*sizeof(double));
        newpop[i].vars = (double*)malloc(nvar*sizeof(double));
        newpop[i].objs = (double*)malloc(nobj*sizeof(double));        
    }    
    solutions = (double**)malloc((2*popsize)*sizeof(double*));
    for(i = 0; i < 2*popsize; i++)
        solutions[i] = (double*)malloc(nobj*sizeof(double));
    
    weights = (double**)malloc(nobj*sizeof(double*));
    for(i = 0; i < nobj; i++)
        weights[i] = (double*)malloc(popsize*sizeof(double));
    
    lowBound = (double*)malloc(nvar*sizeof(double));
    upBound  = (double*)malloc(nvar*sizeof(double));
    aux = (int*)malloc(popsize*sizeof(int));
    
    /* allocate memory for Kuhn-Munkres procedure */
    InitMemoryKM(2*popsize, nobj, popsize);           
}

void FreeMemory()
{
    int i;    
    for (i = 0; i < popsize; i++)
    {
        free(oldpop[i].vars);
        free(oldpop[i].objs);
        free(newpop[i].vars);
        free(newpop[i].objs);
    }    
    for(i = 0; i < 2*popsize; i++) 
        free(solutions[i]);
    /*for(i = 0; i < nobj; i++)
        free(weights[i]);*/
        
    free(newpop);
    free(oldpop);
    free(lowBound);
    free(upBound);
    free(weights);
    free(solutions);
    free(aux);
    
    /* deallocate memory for Kuhn-Munkres procedure */
    FreeMemoryKM(); 
}

/* generate a real random value in [low, up] */
double rnd(double low, double up)
{
    return genrand_real1()*(up - low) + low;
}

/* generate an integer random value in [low, up] */
int rndint(int low, int up)
{
    return (int)rnd(low, up);
}

/* Function to display the current population */
void onthefly_display(int generation)
{
    int i;
    FILE *fpt;
    fpt = fopen("plot.out","w");
    for (i = 0; i < popsize; i++)
    {
        if (nobj == 2)
            fprintf(fpt,"%e\t%e\n",oldpop[i].objs[0],oldpop[i].objs[1]);
        else
            fprintf(fpt,"%e\t%e\t%e\n",oldpop[i].objs[0],oldpop[i].objs[1],oldpop[i].objs[2]);
        fflush(fpt);               
    }
    
    if (nobj == 2)
        fprintf(gp,"set title 'Generation #%d'\n unset key\n plot 'plot.out' w points pointtype 7 pointsize 0.8\n",generation);
    else
        fprintf(gp,"set title 'Generation #%d'\n set ticslevel 0\n set view %d,%d\n unset key\n splot 'plot.out' w points pointtype 7 pointsize 0.8\n",generation,55,112);
    fflush(gp);
    
    fclose(fpt);
    return;
}
