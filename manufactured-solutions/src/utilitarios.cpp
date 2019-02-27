//utilitarios.cpp

#include <cstdlib>
#include <iostream>

using std::cout;

#include <stdio.h>
#include <math.h>

#include "sistema.h"
#include "solvers.h"
#include "utilitarios.h"
#include "funcoes.h"



void erro(double ***Pex, double ***P, double ***E, int nx, int ny, int nz)
{
     int i, j, k;
        
     for(i=0; i<nx; i++)
          for(j=0; j<ny; j++)
               for(k=0; k<nz; k++)
         	        E[i][j][k] = Pex[i][j][k] - P[i][j][k];
}

void residuo(double ***R, double ***P, double ***rhs, double ***Ap, double ***Ae, double ***Aw, double ***An, double ***As, double ***Au, double ***Ad, int nx, int ny, int nz)
{
     int i, j, k;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
		       for(k=1; k<(nz-1); k++)
                    R[i][j][k] = rhs[i][j][k] - (Ae[i][j][k]*P[i+1][j][k] + Aw[i][j][k]*P[i-1][j][k] + An[i][j][k]*P[i][j+1][k] + As[i][j][k]*P[i][j-1][k] + Au[i][j][k]*P[i][j][k+1] + Ad[i][j][k]*P[i][j][k-1] + Ap[i][j][k]*P[i][j][k]);
}

double maximo(double ***R, int nx, int ny, int nz)
{
     int i, j, k;

     double max = 0;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
		       for(k=1; k<(nz-1); k++)
               {
                    if(fabs(R[i][j][k]) > max)
                         max = fabs(R[i][j][k]);
               }
               
     return max;
}

double L2(double ***M, int nx, int ny, int nz, double h)
{
     int i, j, k;
     
	 double l2 = 0;

	 for(i=1; i<(nx-1); i++)
	      for(j=1; j<(ny-1); j++)
	           for(k=1; k<(nz-1); k++)
			        l2 += pow(M[i][j][k],2.0)*pow(h,3.0);

     l2 = sqrt(l2);

	 return l2;
}

void zero(double ****e, int *nx, int *ny, int *nz, int niveis)
{
     int t, i, j, k;

     for(t=0; t<niveis; t++)
          for(i=0; i<nx[t]; i++)
               for(j=0; j<ny[t]; j++)
                    for(k=0; k<nz[t]; k++)
                         e[t][i][j][k] = 0;
}

void corrige(double ***P, double ***e, int nx, int ny, int nz)
{
     int i, j, k;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
		       for(k=1; k<(nz-1); k++)
                    P[i][j][k] += e[i][j][k];
}

void normaliza(double ***P, int nx, int ny, int nz, double h)
{
     int i, j, k;
	 double soma = P[1][1][1];	 

	 for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               for(k=1; k<(nz-1); k++)
			        P[i][j][k] = P[i][j][k] - soma;
}

void igual(double ***rhs, double ***R, int nx, int ny, int nz)
{
     int i, j, k;

     for(i=0; i<nx; i++)
          for(j=0; j<ny; j++)
               for(k=0; k<nz; k++)
                    rhs[i][j][k] = R[i][j][k];
}

double escalar(double ***M1, double ***M2, int nx, int ny, int nz)
{
     int i, j, k;
	 double produto = 0;

	 for(i=1; i<(nx-1); i++)
	      for(j=1; j<(ny-1); j++)
	           for(k=1; k<(nz-1); k++)
			        produto += M1[i][j][k]*M2[i][j][k];

	 return produto;
}


void interpola3d(double ****e, int *nx, int *ny, int *nz, int t)
{
     
     int i, j, k;

     int limite_x = nx[t]-2;
	 int limite_y = ny[t]-2;
	 int limite_z = nz[t]-2;


     for(i=1; i<=limite_x; i++)
          for(j=1; j<=limite_y; j++)
               for(k=1; k<=limite_z; k++)
               {
                    
                    e[t+1][2*i-1][2*j-1][2*k-1] = (1.0/64.0)*e[t][i-1][j-1][k-1] + (3.0/64.0)*e[t][i][j-1][k-1] + (9.0/64.0)*e[t][i][j][k-1] + (3.0/64.0)*e[t][i-1][j][k-1] + (3.0/64.0)*e[t][i-1][j-1][k] + (9.0/64.0)*e[t][i][j-1][k] + (27.0/64.0)*e[t][i][j][k] + (9.0/64.0)*e[t][i-1][j][k];    
                    
                    e[t+1][2*i-1][2*j-1][2*k] = (3.0/64.0)*e[t][i-1][j-1][k] + (9.0/64.0)*e[t][i][j-1][k] + (27.0/64.0)*e[t][i][j][k] + (9.0/64.0)*e[t][i-1][j][k] + (1.0/64.0)*e[t][i-1][j-1][k+1] + (3.0/64.0)*e[t][i][j-1][k+1] + (9.0/64.0)*e[t][i][j][k+1] + (3.0/64.0)*e[t][i-1][j][k+1];
                        
                    e[t+1][2*i-1][2*j][2*k-1] = (3.0/64.0)*e[t][i-1][j][k-1] + (9.0/64.0)*e[t][i][j][k-1] + (3.0/64.0)*e[t][i][j+1][k-1] + (1.0/64.0)*e[t][i-1][j+1][k-1] + (9.0/64.0)*e[t][i-1][j][k] + (27.0/64.0)*e[t][i][j][k] + (9.0/64.0)*e[t][i][j+1][k] + (3.0/64.0)*e[t][i-1][j+1][k];
                    
                    e[t+1][2*i-1][2*j][2*k] = (9.0/64.0)*e[t][i-1][j][k] + (27.0/64.0)*e[t][i][j][k] + (9.0/64.0)*e[t][i][j+1][k] + (3.0/64.0)*e[t][i-1][j+1][k] + (3.0/64.0)*e[t][i-1][j][k+1] + (9.0/64.0)*e[t][i][j][k+1] + (3.0/64.0)*e[t][i][j+1][k+1] + (1.0/64.0)*e[t][i-1][j+1][k+1];
		      
		            e[t+1][2*i][2*j-1][2*k-1] = (3.0/64.0)*e[t][i][j-1][k-1] + (1.0/64.0)*e[t][i+1][j-1][k-1] + (3.0/64.0)*e[t][i+1][j][k-1] + (9.0/64.0)*e[t][i][j][k-1] + (9.0/64.0)*e[t][i][j-1][k] + (3.0/64.0)*e[t][i+1][j-1][k] + (9.0/64.0)*e[t][i+1][j][k] + (27.0/64.0)*e[t][i][j][k];
		            
		            e[t+1][2*i][2*j-1][2*k] = (9.0/64.0)*e[t][i][j-1][k] + (3.0/64.0)*e[t][i+1][j-1][k] + (9.0/64.0)*e[t][i+1][j][k] + (27.0/64.0)*e[t][i][j][k] + (3.0/64.0)*e[t][i][j-1][k+1] + (1.0/64.0)*e[t][i+1][j-1][k+1] + (3.0/64.0)*e[t][i+1][j][k+1] + (9.0/64.0)*e[t][i][j][k+1];
		            
		            e[t+1][2*i][2*j][2*k-1] = (9.0/64.0)*e[t][i][j][k-1] + (3.0/64.0)*e[t][i+1][j][k-1] + (1.0/64.0)*e[t][i+1][j+1][k-1] + (3.0/64.0)*e[t][i][j+1][k-1] + (27.0/64.0)*e[t][i][j][k] + (9.0/64.0)*e[t][i+1][j][k] + (3.0/64.0)*e[t][i+1][j+1][k] + (9.0/64.0)*e[t][i][j+1][k];
		            
		            e[t+1][2*i][2*j][2*k] = (27.0/64.0)*e[t][i][j][k] + (9.0/64.0)*e[t][i+1][j][k] + (3.0/64.0)*e[t][i+1][j+1][k] + (9.0/64.0)*e[t][i][j+1][k] + (9.0/64.0)*e[t][i][j][k+1] + (3.0/64.0)*e[t][i+1][j][k+1] + (1.0/64.0)*e[t][i+1][j+1][k+1] + (3.0/64.0)*e[t][i][j+1][k+1];
		            		            
               } 
}



void restringe3d(double ***R, double ****rhs, int *nx, int *ny, int *nz, int t)
{

	 int i, j, k;
	
	 int limite_x = nx[t-1]-2;
     int limite_y = ny[t-1]-2;
     int limite_z = nz[t-1]-2;

     for(i=1; i<=limite_x; i++)
          for(j=1; j<=limite_y; j++)
               for(k=1; k<=limite_z; k++)
                    rhs[t-1][i][j][k] = 0.125*(R[2*i-1][2*j-1][2*k-1] + R[2*i-1][2*j-1][2*k] + R[2*i-1][2*j][2*k-1] + R[2*i-1][2*j][2*k] + R[2*i][2*j-1][2*k-1] + R[2*i][2*j-1][2*k] + R[2*i][2*j][2*k-1] + R[2*i][2*j][2*k]);
	
}





