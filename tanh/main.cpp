/*

"Mas graças a Deus que nos dá a vitória por nosso Senhor Jesus Cristo.
Portanto, meus amados irmãos, mantenham-se firmes, e que nada os abale. 
Sejam sempre dedicados à obra do Senhor, pois vocês sabem que, no Senhor,
o trabalho de vocês não será inútil." 1Co 15.57,58

BY: VITOR MACIEL VILELA FERREIRA

*/

#include <cstdlib>
#include <iostream>
using std::cout;
using std::cin;
using std::endl;

#include <stdio.h>
#include <math.h>

#include "system.h"
#include "solvers.h"
#include "tools.h"
#include "exact.h"

using namespace std;

int main(int argc, char *argv[])
{
     int pause;       
    
     int t, i, j, k; 
     
     const char boundaryType = 'n';
     const char preconditionType = 's';
     
     
     //COARSER GRID - Y DIRECTION = 2 - 4 - 8 - 16 - 32 - 64 - 128 - 256 (+2)
	 const double Y = 1.0;	
	 
	 int cgy = 4;
	 
	 const double H = Y/(cgy-2);	 
	 const int multX = 1;
	 const int multZ = 1;	 
     int cgx = cgy - 2;
     int cgz = cgy - 2;     
	 for(int m=1; m<multX; m++)
          cgx = 2*cgx;
     cgx = cgx + 2;
     for(int m=1; m<multZ; m++)
          cgz = 2*cgz;
     cgz = cgz + 2;          
	 const double X = pow(2,(multX-1))*Y;	
	 const double Z = pow(2,(multZ-1))*Y;
     
     
     //DENSITY CONSTANTS
     const double Xo = 0.5;
     const double Yo = 0.5;
     const double Zo = 0.5;
     const double ro = 0.125*X;
     const double d = 80.0;
     const double W0 = 1.0;
     const double W1 = 1000.0;
         
         		 
	 //MULTIGRID CONSTANTS    	
	 const int levels = 5;
	 
	 double *w = new double [levels];
	 int *upTimes = new int [levels];
     int *downTimes = new int [levels];
     
     w[0] = 1.0;               downTimes[0] = 6;               upTimes[0] = 6;
     w[1] = 1.0;               downTimes[1] = 2;               upTimes[1] = 2;
     w[2] = 1.0;               downTimes[2] = 2;               upTimes[2] = 2;
     w[3] = 1.0;               downTimes[3] = 2;               upTimes[3] = 2; 
     w[4] = 1.0;               downTimes[4] = 2;               upTimes[4] = 2;                      
    	 
	 int tempX = cgx - 2;
	 int tempY = cgy - 2;
	 int tempZ = cgz - 2;

	 int *nx = new int [levels];
	 int *ny = new int [levels];
	 int *nz = new int [levels];
	 double *h = new double [levels];
	 for(t=0; t<levels; t++)
	 {
          h[t] = Y/tempY;
		  nx[t] = tempX + 2;
		  ny[t] = tempY + 2;
		  nz[t] = tempZ + 2;
		  tempX = 2*tempX;
		  tempY = 2*tempY;
          tempZ = 2*tempZ;		
	 } 
	 	 
	 //4D VECTORS
	 double ****Ap = new double ***[levels];
     double ****Ae = new double ***[levels];
     double ****Aw = new double ***[levels];
     double ****An = new double ***[levels];
     double ****As = new double ***[levels];
	 double ****Af = new double ***[levels];
	 double ****Ab = new double ***[levels];
	 double ****Ro = new double ***[levels];	 
	 
	 double ****Pex = new double ***[levels];
	 double ****P = new double ***[levels];
     double ****B = new double ***[levels];
	 double ****E = new double ***[levels];
     	
	 double ****R = new double ***[levels];	 
	 double ****e = new double ***[levels];     
     double ****rhs = new double ***[levels]; 
     
     double ****z1 = new double ***[levels];
     double ****z2 = new double ***[levels];
     double ****R1 = new double ***[levels];
     double ****R2 = new double ***[levels];         
     double ****p = new double ***[levels];
     double ****W = new double ***[levels];
     
     for(t=0; t<levels; t++)
     {
          Ap[t] = new double **[nx[t]];
          Ae[t] = new double **[nx[t]];
          Aw[t] = new double **[nx[t]];
          An[t] = new double **[nx[t]];
          As[t] = new double **[nx[t]];
		  Af[t] = new double **[nx[t]];
		  Ab[t] = new double **[nx[t]];
		  Ro[t] = new double **[nx[t]];		  
		  
		  Pex[t] = new double **[nx[t]];
          P[t] = new double **[nx[t]];
          B[t] = new double **[nx[t]];
          E[t] = new double **[nx[t]];
          
          R[t] = new double **[nx[t]];          
          e[t] = new double **[nx[t]];
          rhs[t] = new double **[nx[t]];
          
          z1[t] = new double **[nx[t]];
          z2[t] = new double **[nx[t]];
          R1[t] = new double **[nx[t]];
          R2[t] = new double **[nx[t]];
          p[t] = new double **[nx[t]];          
          W[t] = new double **[nx[t]];	 
    		  
          for(i=0; i<nx[t]; i++)
          {
               Ap[t][i] = new double *[ny[t]];
               Ae[t][i] = new double *[ny[t]];
               Aw[t][i] = new double *[ny[t]];
               An[t][i] = new double *[ny[t]];
               As[t][i] = new double *[ny[t]];
			   Af[t][i] = new double *[ny[t]];
			   Ab[t][i] = new double *[ny[t]];
			   Ro[t][i] = new double *[ny[t]];			   
			   
			   Pex[t][i] = new double *[ny[t]];
			   P[t][i] = new double *[ny[t]];
			   B[t][i] = new double *[ny[t]];
			   E[t][i] = new double *[ny[t]];
			   
			   R[t][i] = new double *[ny[t]];			   
			   e[t][i] = new double *[ny[t]];
			   rhs[t][i] = new double *[ny[t]];
			   
			   z1[t][i] = new double *[ny[t]];
			   z2[t][i] = new double *[ny[t]];
			   R1[t][i] = new double *[ny[t]];
			   R2[t][i] = new double *[ny[t]];
			   p[t][i] = new double *[ny[t]];
			   W[t][i] = new double *[ny[t]];
			  
               for(j=0; j<ny[t]; j++)
               {
				    Ap[t][i][j] = new double [nz[t]];
				    Ae[t][i][j] = new double [nz[t]];
				    Aw[t][i][j] = new double [nz[t]];
				    An[t][i][j] = new double [nz[t]];
				    As[t][i][j] = new double [nz[t]];
				    Af[t][i][j] = new double [nz[t]];
				    Ab[t][i][j] = new double [nz[t]];
				    Ro[t][i][j] = new double [nz[t]];
				   				    
				    Pex[t][i][j] = new double [nz[t]];
				    P[t][i][j] = new double [nz[t]];
				    B[t][i][j] = new double [nz[t]];
				    E[t][i][j] = new double [nz[t]];
				    
				    R[t][i][j] = new double [nz[t]];				    
				    e[t][i][j] = new double [nz[t]];
				    rhs[t][i][j] = new double [nz[t]];
				    
				    z1[t][i][j] = new double [nz[t]];
				    z2[t][i][j] = new double [nz[t]];
				    R1[t][i][j] = new double [nz[t]];
				    R2[t][i][j] = new double [nz[t]];
				    p[t][i][j] = new double [nz[t]];
				    W[t][i][j] = new double [nz[t]];
				    
				    for(k=0; k<nz[t]; k++)
				    {
                         Ap[t][i][j][k] = 0;
                         Ae[t][i][j][k] = 0;
       	                 Aw[t][i][j][k] = 0;
                         An[t][i][j][k] = 0;
                         As[t][i][j][k] = 0;
					     Af[t][i][j][k] = 0;
					     Ab[t][i][j][k] = 0;
					     Ro[t][i][j][k] = 1;
										     
					     Pex[t][i][j][k] = 0;
					     P[t][i][j][k] = 0;
					     B[t][i][j][k] = 0;
					     E[t][i][j][k] = 0;
					     
					     R[t][i][j][k] = 0;					     
					     e[t][i][j][k] = 0;
					     rhs[t][i][j][k] = 0;
					     
                         z1[t][i][j][k] = 0;	
                         z2[t][i][j][k] = 0;
                         R1[t][i][j][k] = 0;
                         R2[t][i][j][k] = 0;
                         p[t][i][j][k] = 0;
                         W[t][i][j][k] = 0;	                         				     
				    }
               }
          }
     } 
            
     EXACT(Pex, h, nx, ny, nz, levels, 's');  
      
	 DENSITY_FIELD(Ro, Xo, Yo, Zo, h, nx, ny, nz, d, ro, W0, W1, levels, 'n');
     
     COEFFICIENTS(Ap, Aw, Ae, An, As, Af, Ab, Ro, h, nx, ny, nz, levels);
     
     FORCE_TERM(B, Pex, Ap, Aw, Ae, An, As, Af, Ab, nx, ny, nz, levels);
     
      	 
     t = levels-1;  
          
                  
     cout << "START?" << "\n\n";
     cin >> pause; 
     
	 
     time_t start, end;
     double dif;
     time (&start);	 
	 
                   
     double eps = 1E-08;
     	 
	 //SOLVERS     	  
     
	 MG(Pex[t], E[t], e, R, rhs, P[t], B[t], Ap, Ae, Aw, An, As, Af, Ab, upTimes, downTimes, X, Y, Z, h, nx, ny, nz, levels, w, eps, boundaryType);

     //CONJUGATE_GRADIENT(Pex[t], P[t], B[t], E[t], z1[t], z2[t], R1[t], R2[t], p[t], W[t], Ap, Ae, Aw, An, As, Af, Ab, X, Y, Z, h, nx, ny, nz, w, eps, levels, boundaryType, preconditionType);
	

     //SOR
     /*
	 int iteration = 0;	 
	 double rsd = 1.0;
     double maximumError = 0.0;
	 
	 FILE *p1, *p2;
	 p1 = fopen("RSDxNC.dat", "w");	
     p2 = fopen("ERRORxNC.dat", "w");
     fprintf(p1, "'NC'       'RSD'\n\n"); 
     fprintf(p2, "'NC'       'EMAX'\n\n"); 
     
     do
     {
          iteration++;
          
          SOR(P[t], B[t], Ap[t], Ae[t], Aw[t], An[t], As[t], Af[t], Ab[t], X, Y, Z, h[t], nx[t], ny[t], nz[t], w[t], upTimes[t], 1, boundaryType);          
          
          RESIDUE(R[t], P[t], B[t], Ap[t], Ae[t], Aw[t], An[t], As[t], Af[t], Ab[t], nx[t], ny[t], nz[t]);
          
          rsd = MAXIMUM(R[t], nx[t], ny[t], nz[t]);
          
          fprintf(p1, "%i       %e\n", iteration, rsd);          
          cout << "\nNC: " << iteration << "      -      rsd: " << rsd << "\n";          
          
          maximumError = INFINITY_NORM(Pex[t], P[t], E[t], nx[t], ny[t], nz[t]);	             	      
	      fprintf(p2, "%i       %e\n", iteration, maximumError);               
     }
     while(rsd > eps);
      
     fclose(p1);
     fclose(p2);
     */
     //SOR


     time (&end);
     dif = difftime (end, start);
     cout << "\n\nTIME = " << dif << "\n";

   
     //ERROR     
     double MAXoo;
     double MAXl2;
     
     MAXoo = INFINITY_NORM(Pex[t], P[t], E[t], nx[t], ny[t], nz[t]);	 
	 cout << "\n\nINFINITY_ERROR = " << MAXoo << endl;
	 MAXl2 = L2_NORM(E[t], h[t], nx[t], ny[t], nz[t]);
	 cout << "\nL2_ERROR = " << MAXl2 << "\n\n";   
    
     system("PAUSE");
     return EXIT_SUCCESS;
     
}




