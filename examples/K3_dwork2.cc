// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
//
// tests that the Frob  function  the same under a change of variables


#include "wrapper.h"
#include "timing.h"
#include "tools.h"
#include <time.h>
#include <cstdlib>

#include <sstream>

using namespace std;
using namespace NTL;

int usage(char name[]){
    cout << "Computes the zeta function of: x^4 + y^4 + z^4 + w^4 + lambda * x*y*z*w == 0" << endl;
    cout << "Usage: " << name << " <p> <lambda>" << endl;
    return 1;
}



int main(int argc, char* argv[]) {
    if(argc != 3) 
        return usage(argv[0]);
    int64_t p = atoi(argv[1]);
    if( not ProbPrime(ZZ(p)) ) {
        cout << p << " is not a prime!" << endl << endl; 
        usage(argv[0]);
    }
    int64_t lambda = atoi(argv[2]);
    map< Vec<int64_t>, int64_t , vi64less> f;

    // a very lazy way to create f
    for(int64_t i = 0; i < 4; ++i) {
        Vec<int64_t> v;
        v.SetLength(4, 0);
        v[i] = 4;
        f[v] = 1;
    }
    Vec<int64_t> v;
    v.SetLength(4, 1);
    f[v] = lambda;

    ZZX cp;
    zeta_function(cp, f, p, true);
    cout << cp << endl;
    return 0;
}

