// Copyright 2013-2018 Edgar Costa
// See LICENSE file for license details.


#include "wrapper.h"
#include <cstdlib>
#include <sstream>

using namespace std;
using namespace NTL;

int main()
{

  char buffer[] = "23\n[[4 0 0 0] [1  3  0  0]  [0  0  1  3]  [3  1  0  0]  [0  0  3  1]]\n[1 1  1  1  1]";
  ZZX cp;
  zeta_function(cp, buffer, true, 1);
  cout << cp << endl;

  return 0;
}
