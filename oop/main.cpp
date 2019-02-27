#include <cstdlib>
#include <iostream>
using  std::cout;
using std::endl;
using std::cin;

#include "Domain.h"
#include "Variable.h"
#include "Cell.h"

int main(int argc, char *argv[])
{
   //TESTS
   
   const Domain field;
   cout << "WEST = " << field.getWest() << endl;
   cout << "EAST = " << field.getEast() << endl;
   
   Variable temp("Temperature", "Centered");
   cout << "NAME = " << temp.getName() << endl;
   cout << "TYPE = " << temp.getType() << endl;
   
   Cell first(0.0, 0.25, 0.5, 0.0);
   cout << "WF = " << first.getWface() << endl;
   cout << "C = " << first.getCenter() << endl;
   cout << "EF = " << first.getEface() << endl;
   cout << "VALUE = " << first.getValue() << endl;
   
   
   int i;
   cin >> i; 
   
   
   
   
   
   system("PAUSE");
   return EXIT_SUCCESS;
}
