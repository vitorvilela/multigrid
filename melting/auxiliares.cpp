//AUXILIARES

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


void Set_Flag(char **flag, int nx, int ny)
{
     int i, j;
     
     flag[1][ny-2] = 'A';
     flag[nx-2][ny-2] = 'B';
     flag[1][1] = 'C';
     flag[nx-2][1] = 'D';
     
     for(j=2; j<ny-2; j++)
     {
          flag[1][j] = 'W';
          flag[nx-2][j] = 'E';
     }
     
     for(i=2; i<nx-2; i++)
     {
          flag[i][1] = 'S';
          flag[i][ny-2] = 'N';
     }
     
     for(i=2; i<nx-2; i++)
          for(j=2; j<ny-2; j++)
               flag[i][j] = 'P';
}


void Set_Initial_Temperature(double **T, double **Told, int nx, int ny, double Ti)
{
     int i, j;
              
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
          {
               T[i][j] = Ti;                    
               Told[i][j] = Ti;                                        
          }
}

void Homogeneous_Neumann(double **M, int nx, int ny)
{
     int i, j;
     
     for(j=0; j<ny; j++)
     {          
          M[0][j] = M[1][j];
          M[nx-1][j] = M[nx-2][j];          
     }
     
     for(i=0; i<nx; i++)
     {          
          M[i][0] = M[i][1];
          M[i][ny-1] = M[i][ny-2];        
     } 
}

