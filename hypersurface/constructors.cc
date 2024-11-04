// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

void hypersurface::init_after_dR() {
    p = dR->p;
    precision = dR->precision;
    n = dR->n;
    d = dR->d;
    verbose = dR->verbose;

    fpow.SetLength(1);
    int64_t i;
    Vec<int64_t> u;
    u.SetLength(n + 1);
    for( i = 0; i <= n ; i++)
        u[i] = 0;
    NTL::set( fpow[0][u] );
    fpow[0][u]=1;
}

hypersurface::hypersurface(int64_t p, int64_t precision, int64_t n, int64_t d, bool verbose) {
    dR = make_shared<de_Rham_local>(p, precision, n, d, verbose);
    init_after_dR();

}

hypersurface::hypersurface(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<ZZ_p> f_vector, bool verbose) {
    dR = make_shared<de_Rham_local>(p, precision, n, d, f_vector, verbose);
    init_after_dR();
}



hypersurface::hypersurface(int64_t p, int64_t precision, map< Vec<int64_t>, ZZ_p, vi64less> f, bool verbose) {
    dR = make_shared<de_Rham_local>(p, precision, f, verbose);
    init_after_dR();
}

hypersurface::hypersurface(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<zz_p> fbar_vector, bool verbose) {
    dR = make_shared<de_Rham_local>(p, precision, n, d, fbar_vector, verbose);
    init_after_dR();
}



hypersurface::hypersurface(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose) {
    dR = make_shared<de_Rham_local>(p, precision, fbar, verbose);
    init_after_dR();
}

hypersurface::hypersurface(const char * filename) {
    dR = make_shared<de_Rham_local>(filename);
    init_after_dR();
}
