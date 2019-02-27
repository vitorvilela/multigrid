//funcoes.h

#ifndef funcoesH
#define funcoesH

void fronteira(double ***P, double h, int nx, int ny, int nz, double D, double L, double S, int fun, char fron, int sit);

void E1(double ***Pex, double h, int nx, int ny, int nz, char imp);

void D1(double ***P, double h, int nx, int ny, int nz, double D, double L, double S, int sit);

void N1(double ***P, double h, int nx, int ny, int nz, double D, double L, double S, int sit);

#endif



