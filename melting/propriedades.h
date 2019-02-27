//PROPRIEDADES

#ifndef propriedadesH
#define propriedadesH

void set_K(double ***K, double **T, int *nx, int *ny, int niveis);

void set_Cp0(double **Cp0, double **T, int nx, int ny);

void set_Cp(double ***Cp, double **Cp0, double **T, double **Told, int *nx, int *ny, int niveis);

void set_f0(double **f0, double **T, int nx, int ny);

void set_f(double **f, double **f0, double **T, double **Told, int nx, int ny);




#endif


