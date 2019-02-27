//system.cpp

#include <cstdlib>
#include <stdio.h>
#include <math.h>


void DENSITY_FIELD(double ****Ro, double Xo, double Yo, double Zo, double *h, int *nx, int *ny, int *nz, double d, double ro, double W0, double W1, int levels, char imp)
{
     int t, i, j, k;
     
     double x, y, z, r1, r2, r3, r4, r5;
     
     
     for(t=0; t<levels; t++)
          for(i=1; i<(nx[t]-1); i++)
          {
               x = 0.5*h[t] + (i-1)*h[t];
               
               for(j=1; j<(ny[t]-1); j++)
               {
                    y = 0.5*h[t] + (j-1)*h[t];    
                        
                    for(k=1; k<(nz[t]-1); k++)
	                { 
                         z = 0.5*h[t] + (k-1)*h[t];  
                     
                         r1 = sqrt(pow(x-Xo,2.0) + pow(y-Yo,2.0) + pow(z-Zo,2.0));
                         
                         r2 = sqrt(pow(x-0.2,2.0) + pow(y-0.2,2.0) + pow(z-0.2,2.0));
                         
                         r3 = sqrt(pow(x-0.8,2.0) + pow(y-0.8,2.0) + pow(z-0.8,2.0));
                         
                         r4 = sqrt(pow(x-0.2,2.0) + pow(y-0.8,2.0) + pow(z-0.8,2.0));
                         
                         r5 = sqrt(pow(x-0.8,2.0) + pow(y-0.8,2.0) + pow(z-0.2,2.0));
                                                					     
                         Ro[t][i][j][k] = W0 + W1*0.5*(1.0 - tanh(d*(r1-ro))) + W1*0.5*(1.0 - tanh(d*(r2-ro))) + W1*0.5*(1.0 - tanh(d*(r3-ro))) + W1*0.5*(1.0 - tanh(d*(r4-ro))) + W1*0.5*(1.0 - tanh(d*(r5-ro)));                         
		             }
               }
          }     
	
	 if(imp == 's')
	 { 
          t = levels-1;  
            
          FILE *out;
          char title[15];
 	      sprintf(title, "Ro %03i.dat", nx[t]-2);
       	  out = fopen(title, "w");
          fprintf(out, "variables = \"x\" , \"y\" , \"z\" , \"Ro\" \n");
          fprintf(out, "zone i=%04d j=%04d k=%04d f=point \n", nx[t]-2, ny[t]-2, nz[t]-2);
        		  
          y = 0.0;
          for(j=1; j<(ny[t]-1); j++)
          {
               x = 0.0;
               for(i=1; i<(nx[t]-1); i++)
               {
		            z = 0.0;
				    for(k=1; k<(nz[t]-1); k++)
				    {
                         fprintf(out, "%e   %e   %e   %e\n", x, y, z, Ro[t][i][j][k]);
                         z += h[t];
				    }
				    x += h[t];
                }
                y += h[t];
          }
          fclose(out);
     }   
}


void COEFFICIENTS(double ****Ap, double ****Aw, double ****Ae, double ****An, double ****As, double ****Af, double ****Ab, double ****Ro, double *h, int *nx, int *ny, int *nz, int levels)
{
     int t, i, j, k;
     
     double x, y, z;
     
     for(t=0; t<levels; t++)
          for(i=1; i<(nx[t]-1); i++)                       
                for(j=1; j<(ny[t]-1); j++)               
			        for(k=1; k<(nz[t]-1); k++)
                    {                          
                         An[t][i][j][k] = -(0.5*(1/(Ro[t][i][j+1][k])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         As[t][i][j][k] = -(0.5*(1/(Ro[t][i][j-1][k])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         Aw[t][i][j][k] = -(0.5*(1/(Ro[t][i-1][j][k])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         Ae[t][i][j][k] = -(0.5*(1/(Ro[t][i+1][j][k])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         Af[t][i][j][k] = -(0.5*(1/(Ro[t][i][j][k+1])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         Ab[t][i][j][k] = -(0.5*(1/(Ro[t][i][j][k-1])+1/(Ro[t][i][j][k])))*1/(h[t]*h[t]);
                         Ap[t][i][j][k] = -(Ae[t][i][j][k]+Aw[t][i][j][k]+An[t][i][j][k]+As[t][i][j][k]+Af[t][i][j][k]+Ab[t][i][j][k]);                    
                    }         
}


void FORCE_TERM(double ****B, double ****Pex, double ****Ap, double ****Aw, double ****Ae, double ****An, double ****As, double ****Au, double ****Ad, int *nx, int *ny, int *nz, int levels)
{
     int t, i, j, k;
     
     for(t=0; t<levels; t++)         
          for(i=1; i<(nx[t]-1); i++)
               for(j=1; j<(ny[t]-1); j++)
                    for(k=1; k<(nz[t]-1); k++)	
                         B[t][i][j][k] = Ap[t][i][j][k]*Pex[t][i][j][k]+Ae[t][i][j][k]*Pex[t][i+1][j][k]+
                                         Aw[t][i][j][k]*Pex[t][i-1][j][k]+An[t][i][j][k]*Pex[t][i][j+1][k]+
                                         As[t][i][j][k]*Pex[t][i][j-1][k]+Au[t][i][j][k]*Pex[t][i][j][k+1]+
					                     Ad[t][i][j][k]*Pex[t][i][j][k-1];
}


void RESIDUE(double ***R, double ***P, double ***rhs, double ***Ap, double ***Ae, double ***Aw, double ***An, double ***As, double ***Af, double ***Ab, int nx, int ny, int nz)
{
     int i, j, k;

     for(i=1; i<(nx-1); i++)
          for(j=1; j<(ny-1); j++)
		       for(k=1; k<(nz-1); k++)
                    R[i][j][k] = rhs[i][j][k] - (Ae[i][j][k]*P[i+1][j][k] + Aw[i][j][k]*P[i-1][j][k] + An[i][j][k]*P[i][j+1][k] + As[i][j][k]*P[i][j-1][k] + Af[i][j][k]*P[i][j][k+1] + Ab[i][j][k]*P[i][j][k-1] + Ap[i][j][k]*P[i][j][k]);
}
     

     
