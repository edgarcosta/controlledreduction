// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

hypersurface::hypersurface(int64_t p, int64_t precision, int64_t n, int64_t d, bool verbose)
{
    this->p = p;
    this->precision = precision;
    this->n = n;
    this->d = d;
    this->verbose = verbose;
    dR = make_shared<de_Rham_local>(p, precision, n, d, verbose);
    //dR = new de_Rham_local(p, precision, n, d, verbose);
    fpow.SetLength(1);

    Vec<int64_t> u;
    int64_t i;

    u.SetLength(n + 1);
    for( i = 0; i <= n ; i++)
        u[i] = 0;
    set( fpow[0][u] );
    fpow[0][u]=1;
}

hypersurface::hypersurface(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<zz_p> fbar_vector, bool verbose)
{
    this->p = p;
    this->precision = precision;
    this->n = n;
    this->d = d;
    this->verbose = verbose;
    dR = make_shared<de_Rham_local>(p, precision, n, d, fbar_vector, verbose);
    //dR = new de_Rham_local(p, precision, n, d, fbar_vector, verbose);
    fpow.SetLength(1);

    Vec<int64_t> u;
    int64_t i;

    u.SetLength(n + 1);
    for( i = 0; i <= n ; i++)
        u[i] = 0;
    set( fpow[0][u] );
    fpow[0][u]=1;
}



hypersurface::hypersurface(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose)
{
    this->p = p;
    this->precision = precision;
    dR = make_shared<de_Rham_local>(p, precision, fbar, verbose);
    //dR = new de_Rham_local(p, precision, fbar, verbose);
    n = dR->n;
    d = dR->d;
    this->verbose = verbose;
    fpow.SetLength(1);

    Vec<int64_t> u;
    int64_t i;

    u.SetLength(n + 1);
    for( i = 0; i <= n ; i++)
        u[i] = 0;
    set( fpow[0][u] );
    assert(fpow[0].size() == 1);

}

hypersurface::hypersurface(const char * filename)
{
    dR = make_shared<de_Rham_local>(filename);
    //dR = new de_Rham_local(filename);
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
    set( fpow[0][u] );
    fpow[0][u]=1;
}
