//tools.h

#ifndef toolsH
#define toolsH


double INFINITY_NORM(double ***Pex, double ***P, double ***E, int nx, int ny, int nz);

double MAXIMUM(double ***R, int nx, int ny, int nz);

double L2_NORM(double ***E, double h, int nx, int ny, int nz);

void ZERO(double ****e, int *nx, int *ny, int *nz, int levels);

void CORRECT(double ***P, double ***e, int nx, int ny, int nz);

void INTEGRATION_REFERENCE(double ***P, double h, int nx, int ny, int nz);

void CELL_REFERENCE(double ***P, int nx, int ny, int nz);

void EQUAL(double ***rhs, double ***R, int nx, int ny, int nz);

double INTERNAL_PRODUCT(double ***M1, double ***M2, int nx, int ny, int nz);

void INTERPOLATION(double ****e, int *nx, int *ny, int *nz, int t);

void RESTRICTION(double ***R, double ****rhs, int *nx, int *ny, int *nz, int t);


#endif


