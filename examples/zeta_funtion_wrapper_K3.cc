// Copyright 2013-2018 Edgar Costa
// See LICENSE file for license details.


#include "wrapper.h"
#include <cstdlib>
#include <sstream>

using namespace std;
using namespace NTL;

int main()
{

  //format:
  // prime                  p = 23
  // monomials              x^4 x*y^3 z*w^3 x^3*y z^3*w
  // monomial coefficients  1 1 1 1
  // --> x^4 +  x*y^3 + z*w^3 + x^3*y + z^3*w = 0 over F_23
  char buffer[] = "23\n[[4 0 0 0] [1  3  0  0]  [0  0  1  3]  [3  1  0  0]  [0  0  3  1]]\n[1 1  1  1  1]";
  ZZX cp;
  Mat<ZZ> Frob;
  zeta_function(cp, Frob, buffer, true, 1);
  cout << cp << endl;

  return 0;
}
