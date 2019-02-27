//EULER

#include <math.h>
#include <iostream>
using std::cout;
using std::endl;
using std::cin;

#include "euler.h"


void ue_laminar(int nx, int ny, int nz, double dt, double ***Ro, double mi, double ***Ue, 
                double ***U, double ***V, double ***W, double ***P, double ***Fx, double h)    
{
     int i, j, k;
     
     double dpdx;
     double duudx, dvudy, dwudz;
     double d2udx2, d2udy2, d2udz2;
     
     double ni;
        
     for(i=1; i<nx; i++)
          for(j=1; j<(ny-1); j++)
               for(k=1; k<(nz-1); k++)
               {
                    ni = mi/Ro[i][j][k];      
                          
                    dpdx = (P[i][j][k]-P[i-1][j][k])/h;
                    duudx = (1/(2*h))*( (U[i+1][j][k]+U[i][j][k])*(U[i+1][j][k]+U[i][j][k]) - (U[i][j][k]+U[i-1][j][k])*(U[i][j][k]+U[i-1][j][k]) );
                    dvudy = (1/(2*h))*( (V[i][j+1][k]+V[i-1][j+1][k])*(U[i][j+1][k]+U[i][j][k]) - (V[i][j][k]+V[i-1][j][k])*(U[i][j][k]+U[i][j-1][k]) );
                    dwudz = (1/(2*h))*( (W[i][j][k+1]+W[i-1][j][k+1])*(U[i][j][k+1]+U[i][j][k]) - (W[i][j][k]+W[i-1][j][k])*(U[i][j][k]+U[i][j][k-1]) );
                    d2udx2 = (1/(h*h))*(U[i+1][j][k]-2*U[i][j][k]+U[i-1][j][k]);
                    d2udy2 = (1/(h*h))*(U[i][j+1][k]-2*U[i][j][k]+U[i][j-1][k]); 
                    d2udz2 = (1/(h*h))*(U[i][j][k+1]-2*U[i][j][k]+U[i][j][k-1]);                     
                          
                    Ue[i][j][k] = U[i][j][k] + dt*( (1/Ro[i][j][k])*Fx[i][j][k] + (-1/Ro[i][j][k])*dpdx - (duudx + dvudy + dwudz) + ni*(d2udx2 + d2udy2 + d2udz2) );  
               } 
}


void ve_laminar(int nx, int ny, int nz, double dt, double ***Ro, double mi, double ***Ve, 
                double ***U, double ***V, double ***W, double ***P, double ***Fy, double h) 
{
     int i, j, k;
     
     double dpdy;
     double dvvdy, duvdx, dwvdz;
     double d2vdx2, d2vdy2, d2vdz2;
     
     double ni;
        
     for(i=1; i<(nx-1); i++)
          for(j=1; j<ny; j++)
               for(k=1; k<(nz-1); k++)
               {
                    ni = mi/Ro[i][j][k];       
                          
                    dpdy = (P[i][j][k]-P[i][j-1][k])/h;
                    dvvdy = (1/(2*h))*( (V[i][j+1][k]+V[i][j][k])*(V[i][j+1][k]+V[i][j][k]) - (V[i][j][k]+V[i][j-1][k])*(V[i][j][k]+V[i][j-1][k]) );
                    duvdx = (1/(2*h))*( (U[i+1][j][k]+U[i+1][j-1][k])*(V[i+1][j][k]+V[i][j][k]) - (U[i][j][k]+U[i][j-1][k])*(V[i][j][k]+V[i-1][j][k]) );
                    dwvdz = (1/(2*h))*( (W[i][j][k+1]+W[i][j-1][k+1])*(V[i][j][k+1]+V[i][j][k]) - (W[i][j][k]+W[i][j-1][k])*(V[i][j][k]+V[i][j][k-1]) );
                    d2vdx2 = (1/(h*h))*(V[i+1][j][k]-2*V[i][j][k]+V[i-1][j][k]);
                    d2vdy2 = (1/(h*h))*(V[i][j+1][k]-2*V[i][j][k]+V[i][j-1][k]); 
                    d2vdz2 = (1/(h*h))*(V[i][j][k+1]-2*V[i][j][k]+V[i][j][k-1]);                     
                          
                    Ve[i][j][k] = V[i][j][k] + dt*( (1/Ro[i][j][k])*Fy[i][j][k] + (-1/Ro[i][j][k])*dpdy - (duvdx + dvvdy + dwvdz) + ni*(d2vdx2 + d2vdy2 + d2vdz2) );  
               }	                 
}


void we_laminar(int nx, int ny, int nz, double dt, double ***Ro, double mi, double ***We, 
                double ***U, double ***V, double ***W, double ***P, double ***Fz, double h)    
{
     int i, j, k;
     
     double dpdz;
     double dwwdz, dvwdy, duwdx;
     double d2wdx2, d2wdy2, d2wdz2;
     
     double ni;
        
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               for(k=1; k<nz; k++)
               {
                    ni = mi/Ro[i][j][k];      
                          
                    dpdz = (P[i][j][k]-P[i][j][k-1])/h;
                    dwwdz = (1/(2*h))*( (W[i][j][k+1]+W[i][j][k])*(W[i][j][k+1]+W[i][j][k]) - (W[i][j][k]+W[i][j][k-1])*(W[i][j][k]+W[i][j][k-1]) );
                    dvwdy = (1/(2*h))*( (V[i][j+1][k]+V[i][j+1][k-1])*(W[i][j+1][k]+W[i][j][k]) - (V[i][j][k]+V[i][j][k-1])*(W[i][j][k]+W[i][j-1][k]) );
                    duwdx = (1/(2*h))*( (U[i+1][j][k]+U[i+1][j][k-1])*(W[i+1][j][k]+W[i][j][k]) - (U[i][j][k]+U[i][j][k-1])*(W[i][j][k]+W[i-1][j][k]) );
                    d2wdx2 = (1/(h*h))*(W[i+1][j][k]-2*W[i][j][k]+W[i-1][j][k]);
                    d2wdy2 = (1/(h*h))*(W[i][j+1][k]-2*W[i][j][k]+W[i][j-1][k]); 
                    d2wdz2 = (1/(h*h))*(W[i][j][k+1]-2*W[i][j][k]+W[i][j][k-1]);                     
                          
                    We[i][j][k] = W[i][j][k] + dt*( (1/Ro[i][j][k])*Fz[i][j][k] + (-1/Ro[i][j][k])*dpdz - (duwdx + dvwdy + dwwdz) + ni*(d2wdx2 + d2wdy2 + d2wdz2) );  
               } 
}



void f(int nx, int ny, int nz, double ***B, double ***Ue, double ***Ve, double ***We, double dt, double h)
{       
     int i, j, k;
       
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               for(k=1; k<(nz-1); k++)
                    B[i][j][k] = -(1/(dt*h))*( Ue[i+1][j][k]-Ue[i][j][k] + Ve[i][j+1][k]-Ve[i][j][k] + We[i][j][k+1]-We[i][j][k] );
}


void p(int nx, int ny, int nz, double ***P, double ***Pe)
{       
     int i, j, k;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
                for(k=1; k<(nz-1); k++)
                     P[i][j][k] += Pe[i][j][k];                             
}


void u(int nx, int ny, int nz, double dt, double ***Ro, double ***U, double ***Ue, double ***Pe, double h)    
{
     int i, j, k;
        
     for(i=1; i<nx; i++)
          for(j=1; j<(ny-1); j++)
                for(k=1; k<(nz-1); k++)
                    U[i][j][k] = Ue[i][j][k] - dt/(Ro[i][j][k]*h)*(Pe[i][j][k]-Pe[i-1][j][k]);                              
}


void v(int nx, int ny, int nz, double dt, double ***Ro, double ***V, double ***Ve, double ***Pe, double h) 
{
     int i, j, k;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<ny; j++)
                for(k=1; k<(nz-1); k++)
                    V[i][j][k] = Ve[i][j][k] - dt/(Ro[i][j][k]*h)*(Pe[i][j][k]-Pe[i][j-1][k]);                              
}


void w(int nx, int ny, int nz, double dt, double ***Ro, double ***W, double ***We, double ***Pe, double h) 
{
     int i, j, k;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
                for(k=1; k<nz; k++)
                    W[i][j][k] = We[i][j][k] - dt/(Ro[i][j][k]*h)*(Pe[i][j][k]-Pe[i][j][k-1]);                              
}


void continuidade(double ***U, double ***V, double ***W, int nx, int ny, int nz, double h, int tempo)
{
     int i, j, k;
     
     double CONS = 0;

     double maximo = 0;
	   
     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
               for(k=1; k<(nz-1); k++)
               {
                    CONS = (1/h)*(U[i+1][j][k]-U[i][j][k] + V[i][j+1][k]-V[i][j][k] + W[i][j][k+1]-W[i][j][k]);

                    if(fabs(CONS) > maximo)
                         maximo = fabs(CONS);
               }

     cout << "\n\nTEMPO: " << tempo << " ---> CONS: " << maximo << "\n\n";                                 
}

