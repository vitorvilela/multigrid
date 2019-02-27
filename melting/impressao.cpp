#include <iostream>

using std::cout;
using std::endl;

#include <math.h>
#include <stdio.h>


void perfil_cavidade(double **U, double **V, int nx, int ny, double h, FILE *U_y, FILE *V_x)
{
        int i, j;

        double y;
        int a = ceil((nx-2)/2);

        for(j=1; j<=(ny-2); j++)
        {
                y = 0.5*h + h*(j-1);
                fprintf(U_y, "%14.14lf        %14.14lf  \n", U[a][j], y);
        }
        
        double x;
        int b = ceil((ny-2)/2);

        for(i=1; i<=(nx-2); i++)
        {
                x = 0.5*h + h*(i-1);
                fprintf(V_x, "%14.14lf        %14.14lf  \n", x, V[i][b]);
        } 		
}


void campo(int nx, int ny, double h, double **Cp, double **K, double **T, double **f, int t)
{

   FILE *sai;

   char titulo[20];
   sprintf (titulo, "ENERGIA%08i.dat", t);
   sai = fopen(titulo, "w");

   fprintf(sai, "title= \"%s\" \n", titulo);
   fprintf(sai, "variables = \"x\" , \"y\" , \"Cp\", \"K\", \"T\", \"f\", \n");
   fprintf(sai, "zone i=%04d j=%04d f=point \n", nx-2, ny-2);

   double vadx,vady;

   vady = 0.0;

   for(int j=1; j<(ny-1); j++)
   {

        vadx = 0.0;

        for(int i=1; i<(nx-1); i++)
        {

                fprintf(sai, "%e  %e  %e  %e  %e  %e \n", vadx, vady, Cp[i][j], K[i][j], T[i][j], f[i][j]);
                vadx += h;

        }

        vady += h;

   }

   fclose(sai);

}

