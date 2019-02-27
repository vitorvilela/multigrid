//FRONTEIRAS

#include <cstdlib>
#include <iostream>

using std::cout;
using std::cin;
using std::endl;

#include <stdio.h>
#include <math.h>

#include "sistema.h"
#include "solvers.h"
#include "utilitarios.h"
#include "exatas.h"
#include "fronteiras.h"



void Fronteira_P(double ***P, double h, int nx, int ny, int nz, char fron, int sit)
{
     if(fron == 'p')
          Periodica_P(P, h, nx, ny, nz);
}


void Fronteira_U(double ***U, double h, int nx, int ny, int nz, char fron)
{
     if(fron == 'p')
          Periodica_U(U, h, nx, ny, nz);
}


void Fronteira_V(double ***V, double h, int nx, int ny, int nz, char fron)
{
     if(fron == 'p')
          Periodica_V(V, h, nx, ny, nz);
}


void Fronteira_W(double ***W, double h, int nx, int ny, int nz, char fron)
{
     if(fron == 'p')
          Periodica_W(W, h, nx, ny, nz);
}



void Periodica_P(double ***P, double h, int nx, int ny, int nz)
{
     int i, j, k;
     double x, y, z;   
     
     //normaliza(P, nx, ny, nz, h);  
         
     for(j=0; j<ny; j++)                            
	      for(k=0; k<nz; k++)
	      {               
               P[0][j][k] = P[nx-2][j][k];
               P[nx-1][j][k] = P[1][j][k];               		
          }     

     for(i=0; i<nx; i++)               
	      for(k=0; k<nz; k++)
	      {		      
               P[i][0][k] = P[i][ny-2][k];
               P[i][ny-1][k] = P[i][1][k];           	   
	      }      
     
     for(i=0; i<nx; i++)     
	      for(j=0; j<ny; j++)
	      {        
     		   P[i][j][0] = P[i][j][nz-2];
     		   P[i][j][nz-1] = P[i][j][1];               
	      }      
}


void Periodica_U(double ***U, double h, int nx, int ny, int nz)
{
     int i, j, k;
     double x, y, z;     
         
     for(j=1; j<(ny-1); j++)                            
	      for(k=1; k<(nz-1); k++)
	      {               
               U[0][j][k] = U[nx-2][j][k];
               U[nx][j][k] = U[2][j][k];               		
          }
     

     for(i=1; i<(nx-1); i++)               
	      for(k=1; k<(nz-1); k++)
	      {		      
               U[i][0][k] = U[i][ny-2][k];
               U[i][ny-1][k] = U[i][1][k];           	   
	      }     
     
     for(i=1; i<(nx-1); i++)     
	      for(j=1; j<(ny-1); j++)
	      {        
     		   U[i][j][0] = U[i][j][nz-2];
     		   U[i][j][nz-1] = U[i][j][1];               
	      }      
}


void Periodica_V(double ***V, double h, int nx, int ny, int nz)
{
     int i, j, k;
     double x, y, z;     
         
     for(j=1; j<(ny-1); j++)                            
	      for(k=1; k<(nz-1); k++)
	      {               
               V[0][j][k] = V[nx-2][j][k];
               V[nx-1][j][k] = V[1][j][k];               		
          }     

     for(i=1; i<(nx-1); i++)               
	      for(k=1; k<(nz-1); k++)
	      {		      
               V[i][0][k] = V[i][ny-2][k];
               V[i][ny][k] = V[i][2][k];           	   
	      }      
     
     for(i=1; i<(nx-1); i++)     
	      for(j=1; j<(ny-1); j++)
	      {        
     		   V[i][j][0] = V[i][j][nz-2];
     		   V[i][j][nz-1] = V[i][j][1];               
	      }      
}


void Periodica_W(double ***W, double h, int nx, int ny, int nz)
{
     int i, j, k;
     double x, y, z;
              
     for(j=1; j<(ny-1); j++)                            
	      for(k=1; k<(nz-1); k++)
	      {               
               W[0][j][k] = W[nx-2][j][k];
               W[nx-1][j][k] = W[1][j][k];               		
          }     

     for(i=1; i<(nx-1); i++)               
	      for(k=1; k<(nz-1); k++)
	      {		      
               W[i][0][k] = W[i][ny-2][k];
               W[i][ny-1][k] = W[i][1][k];           	   
	      }
         
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
	      {        
     		   W[i][j][0] = W[i][j][nz-2];
     		   W[i][j][nz] = W[i][j][2];               
	      }     
}


void Dirichlet(int nx, int ny, int nz, double ***U, double ***V, double ***W, double ***Uex, double ***Vex, double ***Wex)
{	
     int i, j, k;
         
     i = 0;
     for(j=1; j<(ny-1); j++)
          for(k=1; k<(nz-1); k++)
               U[i][j][k] = Uex[i][j][k];
    
     i = nx;
     for(j=1; j<(ny-1); j++)
          for(k=1; k<(nz-1); k++)
               U[i][j][k] = Uex[i][j][k];

     j = 0;
     for(i=1; i<(nx-1); i++)
          for(k=1; k<(nz-1); k++)
               V[i][j][k] = Vex[i][j][k];

     j = ny;
     for(i=1; i<(nx-1); i++)
          for(k=1; k<(nz-1); k++)
               V[i][j][k] = Vex[i][j][k];
               
     k = 0;
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               W[i][j][k] = Wex[i][j][k];

     k = nz;
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               W[i][j][k] = Wex[i][j][k];
}


void Update_vel(int nx, int ny, int nz, double ***Ue, double ***Ve, double ***We, double ***Uex, double ***Vex, double ***Wex)
{	
     int i, j, k;
         
     i = 1;
     for(j=1; j<(ny-1); j++)
          for(k=1; k<(nz-1); k++)
               Ue[i][j][k] = Uex[i][j][k];
    
     i = (nx-1);
     for(j=1; j<(ny-1); j++)
          for(k=1; k<(nz-1); k++)
               Ue[i][j][k] = Uex[i][j][k];

     j = 1;
     for(i=1; i<(nx-1); i++)
          for(k=1; k<(nz-1); k++)
               Ve[i][j][k] = Vex[i][j][k];

     j = (ny-1);
     for(i=1; i<(nx-1); i++)
          for(k=1; k<(nz-1); k++)
               Ve[i][j][k] = Vex[i][j][k];
               
     k = 1;
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               We[i][j][k] = Wex[i][j][k];

     k = (nz-1);
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               We[i][j][k] = Wex[i][j][k];
}



void Cavidade_vel(int nx, int ny, int nz, double ***U, double ***V, double ***W, double Ut)
{	
    int i, j, k;
    
    i = 1;
    for(j=1; j<(ny-1); j++)    
         for(k=1; k<(nz-1); k++)
         {
              U[i-1][j][k] = -U[i+1][j][k];
              V[i-1][j][k] = -V[i][j][k];
              W[i-1][j][k] = -W[i][j][k];
         }    
            
    i = (nx-2);
    for(j=1; j<(ny-1); j++)
         for(k=1; k<(nz-1); k++)
         {
              U[i+2][j][k] = -U[i][j][k];
              V[i+1][j][k] = -V[i][j][k];
              W[i+1][j][k] = -W[i][j][k];
         }    
             
    j = 1;
    for(i=1; i<(nx-1); i++)
         for(k=1; k<(nz-1); k++)  
         {       
              U[i][j-1][k] = -U[i][j][k];
              V[i][j-1][k] = -V[i][j+1][k];
              W[i][j-1][k] = -U[i][j][k];
         }    
             
    j = (ny-2);
    for(i=1; i<(nx-1); i++)
         for(k=1; k<(nz-1); k++)
         {
              U[i][j+1][k] = 2*Ut - U[i][j][k]; 
              V[i][j+2][k] = -V[i][j][k];  
              W[i][j+1][k] = -W[i][j][k]; 
         }
         
    k = 1;
    for(j=1; j<(ny-1); j++)    
         for(i=1; i<(nx-1); i++)
         {
              U[i][j][k-1] = -U[i][j][k];
              V[i][j][k-1] = -V[i][j][k];
              W[i][j][k-1] = -W[i][j][k+1];
         }    
            
    k = (nz-2);
    for(j=1; j<(ny-1); j++)
         for(i=1; i<(nx-1); i++)
         {
              U[i][j][k+1] = -U[i][j][k];
              V[i][j][k+1] = -V[i][j][k];
              W[i][j][k+2] = -W[i][j][k];
         }      
}


void Update_vel_cavidade(int nx, int ny, int nz, double ***U, double ***V, double ***W)
{	
     int i, j, k;
         
     i = 1;
     for(j=1; j<(ny-1); j++)
          for(k=1; k<(nz-1); k++)
               U[i][j][k] = 0;
    
     i = (nx-1);
     for(j=1; j<(ny-1); j++)
          for(k=1; k<(nz-1); k++)
               U[i][j][k] = 0;

     j = 1;
     for(i=1; i<(nx-1); i++)
          for(k=1; k<(nz-1); k++)
               V[i][j][k] = 0;

     j = (ny-1);
     for(i=1; i<(nx-1); i++)
          for(k=1; k<(nz-1); k++)
               V[i][j][k] = 0;
               
     k = 1;
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               W[i][j][k] = 0;

     k = (nz-1);
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               W[i][j][k] = 0;
}


void Cavidade_pressao(int nx, int ny, int nz, double ***Pe, double h, int sit)
{
     int i, j, k;
     
     if(sit == 1)
	      normaliza(Pe, nx, ny, nz, h);

     i = 1;
     for(j=0; j<=(ny-1); j++)
          for(k=0; k<=(nz-1); k++)
               Pe[i-1][j][k] = Pe[i][j][k];

     i = (nx-2);
     for(j=0; j<=(ny-1); j++)
          for(k=0; k<=(nz-1); k++)
               Pe[i+1][j][k] = Pe[i][j][k];

     j = 1;
     for(i=0; i<=(nx-1); i++)
          for(k=0; k<=(nz-1); k++)
               Pe[i][j-1][k] = Pe[i][j][k];

     j = (ny-2);
     for(i=0; i<=(nx-1); i++)
          for(k=0; k<=(nz-1); k++)
               Pe[i][j+1][k] = Pe[i][j][k];
               
     k = 1;
     for(i=0; i<=(nx-1); i++)
          for(j=0; j<=(ny-1); j++)
               Pe[i][j][k-1] = Pe[i][j][k];

     k = (nz-2);
     for(i=0; i<=(nx-1); i++)
          for(j=0; j<=(ny-1); j++)
               Pe[i][j][k+1] = Pe[i][j][k];
}


