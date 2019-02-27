//EXATAS

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



void P_EX(double ***Pex, double h, int nx, int ny, int nz, char imp, double t)
{        
     int i, j, k;
     double x, y, z;
     const double Pi = 4*atan(1);
     
     for(i=0; i<nx; i++)
     {
          x = 0.5*h + (i-1)*h;
          for(j=0; j<ny; j++)
          {
               y = 0.5*h + (j-1)*h;			
		       for(k=0; k<nz; k++)
		       {
                    z = 0.5*h + (k-1)*h;
				    Pex[i][j][k] = cos(2*Pi*x+2*Pi*y+2*Pi*z+t);
		       }
          }
     }
     
     if(imp == 's')
     {
          FILE *sai_e;
          char titulo_e[10];
          sprintf(titulo_e, "Pex.dat");
       	  sai_e = fopen(titulo_e, "w");
          fprintf(sai_e,"variables = \"x\" , \"y\" , \"z\" , \"Pex\" \n");
          fprintf(sai_e,"zone i=%04d j=%04d k=%04d f=point \n", nx-2, ny-2, nz-2);
        	
		  double vadx, vady, vadz;
          
          vady = 0.0;
          for(j=1; j<(ny-1); j++)
          {
               vadx = 0.0;
               for(i=1; i<(nx-1); i++)
               {
                    vadz = 0.0;
				    for(k=1; k<(nz-1); k++)
				    {
                         fprintf(sai_e, "%e   %e   %e   %e\n", vadx, vady, vadz, Pex[i][j][k]);
                         vadz += h;
				    }
			        vadx += h;
               }
               vady += h;
          }
          fclose(sai_e);
     }     
}


void U_EX(double ***Uex, double h, int nx, int ny, int nz, char imp, double t)
{        
     int i, j, k;
     double x, y, z;
     const double Pi = 4*atan(1);
     
     for(i=0; i<(nx+1); i++)
     {
          x = (i-1)*h;
          for(j=0; j<ny; j++)
          {
               y = 0.5*h + (j-1)*h;			
		       for(k=0; k<nz; k++)
		       {
                    z = 0.5*h + (k-1)*h;
				    Uex[i][j][k] = pow(sin(2*Pi*x+2*Pi*y+2*Pi*z+t),2);
		       }
          }
     }
     
     if(imp == 's')
     {
          FILE *sai_e;
          char titulo_e[10];
          sprintf(titulo_e, "Uex.dat");
       	  sai_e = fopen(titulo_e, "w");
          fprintf(sai_e,"variables = \"x\" , \"y\" , \"z\" , \"Uex\" \n");
          fprintf(sai_e,"zone i=%04d j=%04d k=%04d f=point \n", nx-2, ny-2, nz-2);
        	
		  double vadx, vady, vadz;
          
          vady = 0.0;
          for(j=1; j<(ny-1); j++)
          {
               vadx = 0.0;
               for(i=1; i<(nx-1); i++)
               {
                    vadz = 0.0;
				    for(k=1; k<(nz-1); k++)
				    {
                         fprintf(sai_e, "%e   %e   %e   %e\n", vadx, vady, vadz, Uex[i][j][k]);
                         vadz += h;
				    }
			        vadx += h;
               }
               vady += h;
          }
          fclose(sai_e);
     }     
}


void V_EX(double ***Vex, double h, int nx, int ny, int nz, char imp, double t)
{        
     int i, j, k;
     double x, y, z;
     const double Pi = 4*atan(1);
     
     for(i=0; i<nx; i++)
     {
          x = 0.5*h + (i-1)*h;
          for(j=0; j<(ny+1); j++)
          {
               y = (j-1)*h;			
		       for(k=0; k<nz; k++)
		       {
                    z = 0.5*h + (k-1)*h;
				    Vex[i][j][k] = pow(cos(2*Pi*x+2*Pi*y+2*Pi*z+t),2);
		       }
          }
     }
     
     if(imp == 's')
     {
          FILE *sai_e;
          char titulo_e[10];
          sprintf(titulo_e, "Vex.dat");
       	  sai_e = fopen(titulo_e, "w");
          fprintf(sai_e,"variables = \"x\" , \"y\" , \"z\" , \"Vex\" \n");
          fprintf(sai_e,"zone i=%04d j=%04d k=%04d f=point \n", nx-2, ny-2, nz-2);
        	
		  double vadx, vady, vadz;
          
          vady = 0.0;
          for(j=1; j<(ny-1); j++)
          {
               vadx = 0.0;
               for(i=1; i<(nx-1); i++)
               {
                    vadz = 0.0;
				    for(k=1; k<(nz-1); k++)
				    {
                         fprintf(sai_e, "%e   %e   %e   %e\n", vadx, vady, vadz, Vex[i][j][k]);
                         vadz += h;
				    }
			        vadx += h;
               }
               vady += h;
          }
          fclose(sai_e);
     }     
}


void W_EX(double ***Wex, double h, int nx, int ny, int nz, char imp, double t)
{        
     int i, j, k;
     double x, y, z;
     const double Pi = 4*atan(1);
     
     for(i=0; i<nx; i++)
     {
          x = 0.5*h + (i-1)*h;
          for(j=0; j<ny; j++)
          {
               y = 0.5*h + (j-1)*h;			
		       for(k=0; k<(nz+1); k++)
		       {
                    z = (k-1)*h;
				    Wex[i][j][k] = 1;
		       }
          }
     }
     
     if(imp == 's')
     {
          FILE *sai_e;
          char titulo_e[10];
          sprintf(titulo_e, "Wex.dat");
       	  sai_e = fopen(titulo_e, "w");
          fprintf(sai_e,"variables = \"x\" , \"y\" , \"z\" , \"Wex\" \n");
          fprintf(sai_e,"zone i=%04d j=%04d k=%04d f=point \n", nx-2, ny-2, nz-2);
        	
		  double vadx, vady, vadz;
          
          vady = 0.0;
          for(j=1; j<(ny-1); j++)
          {
               vadx = 0.0;
               for(i=1; i<(nx-1); i++)
               {
                    vadz = 0.0;
				    for(k=1; k<(nz-1); k++)
				    {
                         fprintf(sai_e, "%e   %e   %e   %e\n", vadx, vady, vadz, Wex[i][j][k]);
                         vadz += h;
				    }
			        vadx += h;
               }
               vady += h;
          }
          fclose(sai_e);
     }     
}


void FX(double ***Fx, double ***Ro, double mi, double h, int nx, int ny, int nz, char imp, double t)
{        
     int i, j, k;
     double x, y, z;
     const double Pi = 4*atan(1);
     
     double dudt, ududx, vdudy, wdudz, dpdx, d2udx2, d2udy2, d2udz2;     
     
     for(i=0; i<nx; i++)
     {
          x = 0.5*h + (i-1)*h;
          for(j=0; j<ny; j++)
          {
               y = 0.5*h + (j-1)*h;			
		       for(k=0; k<nz; k++)
		       {
                    z = 0.5*h + (k-1)*h;
                    
                    dudt = 2*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t);
                    ududx = sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*(4*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi);
                    vdudy = cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*(4*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi);
                    wdudz = 4*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi;
                    dpdx = -2*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi;
                    d2udx2 = 8*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi*Pi-8*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi*Pi;
                    d2udy2 = 8*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi*Pi-8*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi*Pi;
                    d2udz2 = 8*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi*Pi-8*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi*Pi;
                                        
				    Fx[i][j][k] = (dudt+ududx+vdudy+wdudz)+dpdx-mi*(d2udx2+d2udy2+d2udz2);
		       }
          }
     }
     
     if(imp == 's')
     {
          FILE *sai_e;
          char titulo_e[10];
          sprintf(titulo_e, "Fx.dat");
       	  sai_e = fopen(titulo_e, "w");
          fprintf(sai_e,"variables = \"x\" , \"y\" , \"z\" , \"Fx\" \n");
          fprintf(sai_e,"zone i=%04d j=%04d k=%04d f=point \n", nx-2, ny-2, nz-2);
        	
		  double vadx, vady, vadz;
          
          vady = 0.0;
          for(j=1; j<(ny-1); j++)
          {
               vadx = 0.0;
               for(i=1; i<(nx-1); i++)
               {
                    vadz = 0.0;
				    for(k=1; k<(nz-1); k++)
				    {
                         fprintf(sai_e, "%e   %e   %e   %e\n", vadx, vady, vadz, Fx[i][j][k]);
                         vadz += h;
				    }
			        vadx += h;
               }
               vady += h;
          }
          fclose(sai_e);
     }     
}


void FY(double ***Fy, double ***Ro, double mi, double h, int nx, int ny, int nz, char imp, double t)
{        
     int i, j, k;
     double x, y, z;
     const double Pi = 4*atan(1);
     
     double dvdt, udvdx, vdvdy, wdvdz, dpdy, d2vdx2, d2vdy2, d2vdz2;     
     
     for(i=0; i<nx; i++)
     {
          x = 0.5*h + (i-1)*h;
          for(j=0; j<ny; j++)
          {
               y = 0.5*h + (j-1)*h;			
		       for(k=0; k<nz; k++)
		       {
                    z = 0.5*h + (k-1)*h;
                    
                    dvdt = -2*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*sin(2*Pi*x+2*Pi*y+2*Pi*z+t);
                    udvdx = sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*(-4*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi);
                    vdvdy = cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*(-4*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi);
                    wdvdz = -4*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi;
                    dpdy = -2*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi;
                    d2vdx2 = -8*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi*Pi+8*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi*Pi;
                    d2vdy2 = -8*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi*Pi+8*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi*Pi;
                    d2vdz2 = -8*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*cos(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi*Pi+8*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi*Pi;
                                        
				    Fy[i][j][k] = (dvdt+udvdx+vdvdy+wdvdz)+dpdy-mi*(d2vdx2+d2vdy2+d2vdz2);
		       }
          }
     }
     
     if(imp == 's')
     {
          FILE *sai_e;
          char titulo_e[10];
          sprintf(titulo_e, "Fy.dat");
       	  sai_e = fopen(titulo_e, "w");
          fprintf(sai_e,"variables = \"x\" , \"y\" , \"z\" , \"Fy\" \n");
          fprintf(sai_e,"zone i=%04d j=%04d k=%04d f=point \n", nx-2, ny-2, nz-2);
        	
		  double vadx, vady, vadz;
          
          vady = 0.0;
          for(j=1; j<(ny-1); j++)
          {
               vadx = 0.0;
               for(i=1; i<(nx-1); i++)
               {
                    vadz = 0.0;
				    for(k=1; k<(nz-1); k++)
				    {
                         fprintf(sai_e, "%e   %e   %e   %e\n", vadx, vady, vadz, Fy[i][j][k]);
                         vadz += h;
				    }
			        vadx += h;
               }
               vady += h;
          }
          fclose(sai_e);
     }     
}



void FZ(double ***Fz, double ***Ro, double mi, double h, int nx, int ny, int nz, char imp, double t)
{        
     int i, j, k;
     double x, y, z;
     const double Pi = 4*atan(1);
     
     double dwdt, udwdx, vdwdy, wdwdz, dpdz, d2wdx2, d2wdy2, d2wdz2;     
     
     for(i=0; i<nx; i++)
     {
          x = 0.5*h + (i-1)*h;
          for(j=0; j<ny; j++)
          {
               y = 0.5*h + (j-1)*h;			
		       for(k=0; k<nz; k++)
		       {
                    z = 0.5*h + (k-1)*h;
                    
                    dwdt = 0;
                    udwdx = 0;
                    vdwdy = 0;
                    wdwdz = 0;
                    dpdz = -2*sin(2*Pi*x+2*Pi*y+2*Pi*z+t)*Pi;
                    d2wdx2 = 0;
                    d2wdy2 = 0;
                    d2wdz2 = 0;
                                        
				    Fz[i][j][k] = (dwdt+udwdx+vdwdy+wdwdz)+dpdz-mi*(d2wdx2+d2wdy2+d2wdz2);
		       }
          }
     }
     
     if(imp == 's')
     {
          FILE *sai_e;
          char titulo_e[10];
          sprintf(titulo_e, "Fz.dat");
       	  sai_e = fopen(titulo_e, "w");
          fprintf(sai_e,"variables = \"x\" , \"y\" , \"z\" , \"Fz\" \n");
          fprintf(sai_e,"zone i=%04d j=%04d k=%04d f=point \n", nx-2, ny-2, nz-2);
        	
		  double vadx, vady, vadz;
          
          vady = 0.0;
          for(j=1; j<(ny-1); j++)
          {
               vadx = 0.0;
               for(i=1; i<(nx-1); i++)
               {
                    vadz = 0.0;
				    for(k=1; k<(nz-1); k++)
				    {
                         fprintf(sai_e, "%e   %e   %e   %e\n", vadx, vady, vadz, Fz[i][j][k]);
                         vadz += h;
				    }
			        vadx += h;
               }
               vady += h;
          }
          fclose(sai_e);
     }     
}



