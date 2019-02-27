/*
Variable.cpp
*/

#include <cstring>
#include <cassert>
#include "Variable.h"

Variable::Variable(char *VN, char *VT)
{
   setName(VN);
   setType(VT);
}

void Variable::setName(char *VN)
{
   NAME = new char( strlen(VN) + 1 );
   assert(NAME != 0);
   strcpy(NAME, VN);
}
   
void Variable::setType(char *VT)
{
   TYPE = new char( strlen(VT) + 1 );
   assert(TYPE != 0);
   strcpy(TYPE, VT);
}  
