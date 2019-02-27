//AUXILIARES

#ifndef auxiliaresH
#define auxiliaresH

void Set_Flag(char **flag, int nx, int ny);

void Set_Initial_Temperature(double **T, double **Told, int nx, int ny, double Ti);

void Homogeneous_Neumann(double **M, int nx, int ny);


#endif


