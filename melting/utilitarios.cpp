//UTILITARIOS

#include <cstdlib>
#include <iostream>
using std::cout;
using std::cin;
using std::endl;
#include <math.h>
#include <stdio.h>

#include "auxiliares.h"
#include "propriedades.h"
#include "sistema.h"
#include "solvers.h"
#include "utilitarios.h"

double norma(double **P, int nx, int ny, double h)
{
     int i, j;
	 double soma = 0;
	 double NORMA;
	
	 for(i=1; i<(nx-1); i++)
	      for(j=1; j<(ny-1); j++)
		       soma += P[i][j]*P[i][j];

	 soma = soma*h*h;

	 NORMA = sqrt(soma);
	 
	 return NORMA;
}

void normaliza(double **P, int nx, int ny, double h)
{
     int i, j;
	 double soma = 0;
	
	 for(i=1; i<(nx-1); i++)
	      for(j=1; j<(ny-1); j++)
		       soma += P[i][j];

	 soma = soma*h*h;

	 for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
		       P[i][j] = P[i][j] - soma;
}

double maximo(double **R, int nx, int ny)
{
     int i, j;
     double max = 0;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
          {
               if(fabs(R[i][j]) > max)
                    max = fabs(R[i][j]);
          }

     return max;
}

void zero(double ***e, int *nx, int *ny, int niveis)
{
    int t, i, j;

    for(t=0; t<niveis; t++)
         for(i=0; i<nx[t]; i++)
              for(j=0; j<ny[t]; j++)
                   e[t][i][j] = 0;
}

void igual(double **rhs, double **R, int nx, int ny)
{
     int i, j;

     for(i=0; i<nx; i++)
          for(j=0; j<ny; j++)
               rhs[i][j] = R[i][j];
}

void igual3D(double ***M1, double ***M2, int *nx, int *ny, int niveis)
{
     int t, i, j;

     for(t=0; t<(niveis-1); t++)
          for(i=0; i<nx[t]; i++)
               for(j=0; j<ny[t]; j++)
                    M1[t][i][j] = M2[t][i][j];
}


double escalar(double **M1, double **M2, int nx, int ny)
{
     int i, j;
	 double produto = 0;
    
   	 for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
		       produto += M1[i][j]*M2[i][j];
	
	 return produto;
}

void corrigeMAIS(double **P, double **e, int nx, int ny)
{
     int i, j;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               P[i][j] += e[i][j];
}

void corrigeMENOS(double **P, double **e, int nx, int ny)
{
     int i, j;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               P[i][j] -= e[i][j];
}


void interpola(double ***e, int *nx, int *ny, int t)   
{
     int i, j;

     int limite_x = nx[t]-2;
	 int limite_y = ny[t]-2;


     for(i=1; i<=limite_x; i++)
          for(j=1; j<=limite_y; j++)
          {
		       e[t+1][2*i-1][2*j-1] += 0.5*(e[t][i][j] + 0.5*(e[t][i-1][j] + e[t][i][j-1]));
		       
		       e[t+1][2*i-1][2*j] += 0.5*(e[t][i][j] + 0.5*(e[t][i-1][j] + e[t][i][j+1]));
		       
		       e[t+1][2*i][2*j-1] += 0.5*(e[t][i][j] + 0.5*(e[t][i][j-1] + e[t][i+1][j]));
		       
		       e[t+1][2*i][2*j] += 0.5*(e[t][i][j] + 0.5*(e[t][i+1][j] + e[t][i][j+1]));
          }    
}


void restringe(double **R, double ***rhs, int *nx, int *ny, int t)	
{
     int i, j;

     int limite_x = nx[t-1]-2;
	 int limite_y = ny[t-1]-2;

     for(i=1; i<=limite_x; i++)
          for(j=1; j<=limite_y; j++)
               rhs[t-1][i][j] = 0.25*(R[2*i-1][2*j-1] + R[2*i-1][2*j] + R[2*i][2*j-1] + R[2*i][2*j]);
}




