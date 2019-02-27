/*
Cell.cpp
*/

#include "Cell.h"

Cell::Cell(double WF, double C, double EF, double V)
   : nextPtr(0), lastPtr(0)
{
   setWface(WF);
   setCenter(C);
   setEface(EF);
   setValue(V);
   setDx();
}


