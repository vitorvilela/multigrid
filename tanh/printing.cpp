//printing.cpp

#include <cstdlib>
#include <stdio.h>
#include <math.h>

void DENSITY_FIELD(double ***Ro, int nx, int ny, int nz, double h, double Xo, double Yo, double Zo, double d, double ro, double W0, double W1, char imp);
{
     int i, j, k;
     
     double x, y, z, r, f;
     
     for(i=1; i<(nx-1); i++)
     {
          x = 0.5*h + (i-1)*h;
               
          for(j=1; j<(ny-1); j++)
          {
               y = 0.5*h + (j-1)*h;    
                        
               for(k=1; k<(nz-1); k++)
	           {
                    z = 0.5*h + (k-1)*h;  
                     
                    r = sqrt(pow(x-Xo,2.0) + pow(y-Yo,2.0) + pow(z-Zo,2.0));
                    
                    f = 0.5*(1.0 - tanh(d*(r-ro))); 
                            					     
                    Ro[i][j][k] = W0 + W1*f;                         
		        }
          }
     }     
	
	 if(imp == 's')
	 { 
          FILE *out;
          char title[15];
 	      sprintf(title, "Ro %03i.dat", nx-2);
       	  out = fopen(title, "w");
          fprintf(out, "variables = \"x\" , \"y\" , \"z\" , \"Ro\" \n");
          fprintf(out, "zone i=%04d j=%04d k=%04d f=point \n", nx-2, ny-2, nz-2);
        		  
          y = 0.0;
          for(j=1; j<(ny-1); j++)
          {
               x = 0.0;
               for(i=1; i<(nx-1); i++)
               {
		            z = 0.0;
				    for(k=1; k<(nz-1); k++)
				    {
                         fprintf(out, "%e   %e   %e   %e\n", x, y, z, Ro[i][j][k]);
                         z += h;
				    }
				    x += h;
                }
                y += h;
          }
          fclose(out);
     }   
}
