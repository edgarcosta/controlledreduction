// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// hypersurface_nd.h: header file for the class hypersurface_non_degenerate, which contains the routines to compute the Frobenius matrix

#ifndef HYPERSRUFACE_ND_H_
#define HYPERSURFACE_ND_H_

#include "hypersurface.h"
#include "dr_nd.h"
#include <omp.h>
#include <iostream>

using namespace std;
using namespace NTL;


class hypersurface_non_degenerate : public hypersurface{

    public:
        //shared_ptr<de_Rham_non_degenerate_local> dR_ND;
        de_Rham_non_degenerate_local* dR_ND; 
        hypersurface_non_degenerate(){};
        ~hypersurface_non_degenerate(){};
        hypersurface_non_degenerate(int64_t p, int64_t precision, int64_t n, int64_t d, bool verbose);
        hypersurface_non_degenerate(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<zz_p> fbar_vector, bool verbose);
        hypersurface_non_degenerate(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose);

        hypersurface_non_degenerate(const char * filename);
        virtual bool save(const char * filename);


        Vec<ZZ_p> frob_ND(const int64_t coordinate, const int64_t N);

        Vec<ZZ_p> frob_ND_ZZ(const int64_t coordinate, const int64_t N, const int64_t loop);
        Vec<ZZ_p> frob_ND_ZZ_p(const int64_t coordinate, const int64_t N, const int64_t loop);
        Vec<ZZ_p> frob_ND_flint(const int64_t coordinate, const int64_t N);
        Mat<ZZ_p> frob_matrix_ND(Vec<int64_t> N);
};

#endif  // HYPERSRUFACE_ND_H_

