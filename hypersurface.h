// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// hypersurface.h: header file for the class hypersurface, which contains the routines to compute the Frobenius matrix
  


#ifndef HYPERSURFACE_H
#define HYPERSURFACE_H


#include "dr.h"
#include <memory>
#include <iostream>

using namespace std;
using namespace NTL;

class hypersurface{

    public:
        int64_t p;
        int64_t precision;
        int64_t d;
        int64_t n;
        bool verbose;
        //we need it as pointer, so we are able to extend it to a de_Rham_non_degenerate_local
        //shared_ptr<de_Rham_local> dR;
        de_Rham_local* dR;
        /*
         * fpow[j] stores the coefficients of Frob(f^j)
         * st u-> c[u] is the coefficient of x^pu of Frob(f^j)
         */
        Vec< map< Vec<int64_t>, ZZ_p, vi64less> > fpow;

        hypersurface() {};
        virtual ~hypersurface(){};
        
        hypersurface(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<zz_p> fbar_vector, bool verbose);
        hypersurface(int64_t p, int64_t precision, int64_t n, int64_t d, bool verbose);
        hypersurface(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose);
        hypersurface(const char * filename);
        virtual bool save(const char * filename);

        /*
         * assures the array fpow 
         * (where we store the powers of f)
         * has enough powers of f
         */ 
        void compute_fpow(int64_t power);



        /*
         * returns the coordinates of the p-adic approximation of frob(e_i)
         * using N terms
         */
        Vec<ZZ_p> frob_J(const int64_t coordinate, const int64_t N)
        {
            // ZZ_p is broken! FIXME
            return frob_J_ZZ(coordinate, N, 0);
        };
        // working over ZZ_p
        Vec<ZZ_p> frob_J_ZZ_p(const int64_t coordinate, const int64_t N, const int64_t loop);
        // using reduce_vector_J_poly_ZZ
        Vec<ZZ_p> frob_J_ZZ(const int64_t coordinate, const int64_t N, const int64_t loop);

        /*
         * returns the matrix coordinates of the p-adic approximation of frob(e_i)
         * using N_k terms, where (k+1) is the pole order of e_i
         */
        Mat<ZZ_p> frob_matrix_J(Vec<int64_t> N);

};

#endif
