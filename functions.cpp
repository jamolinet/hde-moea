
#include "global.h"

#include "wfg/ExampleProblems.h"
#include "wfg/TransFunctions.h"

# define PI 3.14159265358979
int nwfg;   /* number WFG problem */

using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;
using std::vector;

void zdt1 (double *vars, double *objs, int nvar, int nobj)
{
    double f1, f2, g, h;
    int i;
    f1 = vars[0];
    g = 0.0;
    for (i = 1; i < nvar; i++)
        g += vars[i];
    
    g  = 9.0*g/(nvar - 1);
    g += 1.0;
    h  = 1.0 - sqrt(f1/g);
    f2 = g*h;
    objs[0] = f1;
    objs[1] = f2;
}

void zdt2 (double *vars, double *objs, int nvar, int nobj)
{
    double f1, f2, g, h;
    int i;
    f1 = vars[0];
    g = 0.0;
    for (i = 1; i < nvar; i++)
        g += vars[i];
    
    g  = 9.0*g/(nvar - 1);
    g += 1.0;
    h  = 1.0 - pow((f1/g),2.0);
    f2 = g*h;
    objs[0] = f1;
    objs[1] = f2;
}

void zdt3 (double *vars, double *objs, int nvar, int nobj)
{
    double f1, f2, g, h;
    int i;
    f1 = vars[0];
    g = 0.0;
    for (i = 1; i < nvar; i++)
        g += vars[i];
    
    g  = 9.0*g/(nvar - 1);
    g += 1.0;
    h  = 1.0 - sqrt(f1/g) - (f1/g)*sin(10.0*PI*f1);
    f2 = g*h;
    objs[0] = f1;
    objs[1] = f2;
}

void zdt4 (double *vars, double *objs, int nvar, int nobj)
{
    double f1, f2, g, h;
    int i;
    f1 = vars[0];
    g = 0.0;
    for (i = 1; i < nvar; i++)
        g += vars[i]*vars[i] - 10.0*cos(4.0*PI*vars[i]);
    
    g += 1.0 + 10.0*(nvar - 1);
    h  = 1.0 - sqrt(f1/g);
    f2 = g*h;
    objs[0] = f1;
    objs[1] = f2;
}

void zdt6 (double *vars, double *objs,int nvar,int nobj)
{
    double f1, f2, g, h;
    int i;
    f1 = 1.0 - (exp(-4.0*vars[0]))*pow((sin(6.0*PI*vars[0])),6.0);
    g = 0.0;
    for (i = 1; i < nvar; i++)
        g += vars[i];
    
    g = g/9.0;
    g = pow(g,0.25);
    g = 1.0 + 9.0*g;
    h = 1.0 - pow((f1/g),2.0);
    f2 = g*h;
    objs[0] = f1;
    objs[1] = f2;
}

void dtlz1 (double *vars, double *objs, int nvar, int nobj)
{
    int i, j;
    int k = nvar - nobj + 1;
    double g = 0.0, f;
    for (i = nvar - k + 1; i <= nvar; i++)
        g += pow(vars[i-1]-0.5,2) - cos(20.0 * PI * (vars[i-1]-0.5));
    
    g = 100.0 * (k + g);
    for (i = 1; i <= nobj; i++)
    {
	f = 0.5 * (1 + g);
	for (j = nobj - i; j >= 1; j--)
	    f *= vars[j-1];
	
	if (i > 1)
	    f *= 1.0 - vars[(nobj - i + 1) - 1];
	
	objs[i-1] = f;
    }
}

void dtlz2 (double *vars, double *objs, int nvar, int nobj)
{
    int i, j;
    int k = nvar - nobj + 1;
    double g = 0.0, f;
    for (i = nvar - k + 1; i <= nvar; i++)
	g += pow(vars[i-1]-0.5,2);
    
    for (i = 1; i <= nobj; i++)
    {
	f = (1 + g);
	for (j = nobj - i; j >= 1; j--)
	    f *= cos(vars[j-1] * PI * 0.5);
        
	if (i > 1)
	    f *= sin(vars[(nobj - i + 1) - 1] * PI * 0.5);
        
	objs[i-1] = f;
    }
}

void dtlz3 (double *vars, double *objs, int nvar, int nobj)
{
    int i, j;
    int k = nvar - nobj + 1;
    double g = 0.0, f;
    for (i = nvar - k + 1; i <= nvar; i++)
	g += pow(vars[i-1]-0.5,2) - cos(20 * PI * (vars[i-1]-0.5));

    g = 100 * (k + g);
    for (i = 1; i <= nobj; i++)
    {
	f = (1 + g);
	for (j = nobj - i; j >= 1; j--)
	    f *= cos(vars[j-1] * PI * 0.5);
        
	if (i > 1)
            f *= sin(vars[(nobj - i + 1) - 1] * PI * 0.5);

	objs[i-1] = f;
    }
}

void dtlz4 (double *vars, double *objs, int nvar, int nobj)
{
    int i, j;
    int k = nvar - nobj + 1;
    double g = 0.0, f, alpha = 100;
    for (i = nvar - k + 1; i <= nvar; i++)
	g += pow(vars[i-1]-0.5,2);

    for (i = 1; i <= nobj; i++)
    {
	f = (1 + g);
	for (j = nobj - i; j >= 1; j--)
	    f *= cos(pow(vars[j-1], alpha) * PI * 0.5);
        
	if (i > 1)
	    f *= sin(pow(vars[(nobj - i + 1) - 1], alpha) * PI * 0.5);
        
	objs[i-1] = f;
    }
}

void dtlz5 (double *vars, double *objs, int nvar, int nobj) 
{
    double g, f, p;
    int i, j, t, k;
    double *theta = (double*)malloc((nobj-1)*sizeof(double));
    g = 0.0;
    k = nvar - nobj + 1;
    
    for (i = nvar - k; i < nvar; i++)
        g += pow(vars[i] - 0.5, 2);
    
    theta[0] = vars[0] * PI * 0.5;
    p = PI / (4 * (1 + g));
    
    for(i = 1; i < nobj-1; i++)
        theta[i] = p * (1 + 2 * g * vars[i]);
    
    for(i = 0; i < nobj; i++) 
    {
        f = (1 + g);
        t = nobj - i - 1;
        for(j = 0; j < t; j++) 
            f *= cos(theta[j]);
            
        if(t < nobj - 1) 
            f *= sin(theta[t]);
        
        objs[i] = f;
    }
    free(theta);
}

void dtlz6 (double *vars, double *objs, int nvar, int nobj) 
{
    double g, f, *theta, p;
    int i, j, t, k;
    theta = (double*)malloc((nobj-1)*sizeof(double));
    g = 0.0;
    k = nvar - nobj + 1;
    
    for (i = nvar - k; i < nvar; i++)
        g += pow(vars[i], 0.1);
    
    theta[0] = vars[0] * PI * 0.5;
    p = PI / (4 * (1 + g));
    
    for(i = 1; i < nobj-1; i++)
        theta[i] = p * (1 + 2 * g * vars[i]);
    
    for(i = 0; i < nobj; i++) 
    {
        f = (1 + g);
        t = nobj - i - 1;
        for(j = 0; j < t; j++) 
            f *= cos(theta[j]);
            
        if(t < nobj - 1) 
            f *= sin(theta[t]);
        
        objs[i] = f;
    }
    free(theta);
}

void dtlz7 (double *vars, double *objs, int nvar, int nobj) 
{
    double g, h;
    int i, k;
    g = h = 0.0;
    k = nvar - nobj + 1;
    for(i = 0; i < nobj - 1; i++)
        objs[i] = vars[i];
    
    for(i = nvar - k; i < nvar; i++)
        g += vars[i];
    
    g = 1 + 9 * g / k;
    
    for(i = 0; i < nobj - 1; i++)
        h += vars[i] * (1 + sin(3 * PI * vars[i])) / (1 + g);
    
    h = nobj - h;
    objs[nobj-1] = (1 + g) * h;
}

void wfg (double *vars, double *objs, int nvar, int nobj)
{
    int i, l, k, l_factor = 10, k_factor = 2;
    int n = nvar;
    int m = nobj;
    vector<double> z;
    vector<double> f;
    
    l = l_factor * 2;
    if (m == 2)        
        k = 4;
    else
        k = k_factor * (m - 1);
    
    if (k + l != n)
    {
        printf("ERROR:Number of variables not expected (expected %d)\n",k+l);
        exit(0);
    }
    for (i = 0; i < n; i++)
        z.push_back(vars[i]);
    
    switch (nwfg) 
    {
        case 1:
            f = Problems::WFG1(z, k, m);
            break;
        case 2:
            f = Problems::WFG2(z, k, m);
            break;
        case 3:
            f = Problems::WFG3(z, k, m);
            break;
        case 4:
            f = Problems::WFG4(z, k, m);
            break;
        case 5:
            f = Problems::WFG5(z, k, m);
            break;
        case 6:
            f = Problems::WFG6(z, k, m);
            break;
        case 7:
            f = Problems::WFG7(z, k, m);
            break;
        case 8:
            f = Problems::WFG8(z, k, m);
            break;
        case 9:
            f = Problems::WFG9(z, k, m);
            break;
    }
    
    if (m != (int) f.size())
    {
        printf("ERROR: Return more objectives than expected\n");
        exit(0);
    }
    
    for (i = 0; i < m; i++)
        objs[i] = f[i];        
}

void (*getFunction(char *function))(double*, double*,int,int)
{
	if(strcmp(function,"zdt1") == 0)     { return zdt1; }
	else if(strcmp(function,"zdt2") == 0){ return zdt2; }
	else if(strcmp(function,"zdt3") == 0){ return zdt3; }
	else if(strcmp(function,"zdt4") == 0){ return zdt4; }
	else if(strcmp(function,"zdt6") == 0){ return zdt6; }	
	else if(strcmp(function,"dtlz1") == 0){return dtlz1;}
	else if(strcmp(function,"dtlz2") == 0){return dtlz2;}
	else if(strcmp(function,"dtlz3") == 0){return dtlz3;}
	else if(strcmp(function,"dtlz4") == 0){return dtlz4;}
        else if(strcmp(function,"dtlz5") == 0){return dtlz5;}
        else if(strcmp(function,"dtlz6") == 0){return dtlz6;}
        else if(strcmp(function,"dtlz7") == 0){return dtlz7;}
        else if(strcmp(function,"wfg1") == 0){ nwfg = 1; return wfg; }
	else if(strcmp(function,"wfg2") == 0){ nwfg = 2; return wfg; }
	else if(strcmp(function,"wfg3") == 0){ nwfg = 3; return wfg; }
	else if(strcmp(function,"wfg4") == 0){ nwfg = 4; return wfg; }
	else if(strcmp(function,"wfg5") == 0){ nwfg = 5; return wfg; }	
	else if(strcmp(function,"wfg6") == 0){ nwfg = 6; return wfg; }
	else if(strcmp(function,"wfg7") == 0){ nwfg = 7; return wfg; }
        else if(strcmp(function,"wfg8") == 0){ nwfg = 8; return wfg; }
        else if(strcmp(function,"wfg9") == 0){ nwfg = 9; return wfg; }
}