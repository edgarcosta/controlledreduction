// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// dr_nd.h: header file for de_Rham_non_degenerate_local, which implements the de Rham structure for non degenerate varieties to perform controlled_reductions
   

#ifndef DR_ND_H
#define DR_ND_H

#include "dr.h"


using namespace std;
using namespace NTL;


class de_Rham_non_degenerate_local : public de_Rham_local{
    /*
     * Some notation:
     * let W_u = {x^u ((m-1)! G0 + m! G1 / f + ... + (m+n-1)! Gn / f^n) Omega / (x0 ... xn f^m), Gi \in k[x0, ... , xn]/J_\empty and deg Gi = d*i}
     * and sum(u) = d * m
     */

    public:
        /*
         * the Hilbert polynomial for S = \emtpy
         */
        Vec<int64_t> Hilbert_ND;

        /*
         * the basis of the coKernels of J_\emtpyset
         * at tge degrees multiples of the degree of f
         * basis elements are of the x^t \Omega / f^m
         * with sum(t) = m*d
         */
        Vec< Vec<int64_t> > coKernels_ND_basis;

        /*
         * the dictionary for reverse lookup
         */
        map< Vec<int64_t>, int64_t, vi64less> coKernels_ND_basis_dict;

        /*
         * the equivalent of solve_J (as in dr.h) but for the ND case, ie, for the ideal J_\emptyset
         * solve_ND[i] is the matarix that will tell us how write a monomial of defree i in terms of 
         * xi dfi and the coKenerls elements
         */
        map< int64_t, pair< Vec<int64_t>, Mat<ZZ_p> > > solve_ND;

        /*
         * the dictionary where we will store the inlusion matrices from x^u (l-1)! Gl/ f^l (x0 ... xn) -> (l-1)! H / f^l
         * the factorial keeps being implicit
         * sum(u) = d, u[i] > 0 and l < n
         * deg Gl = l*d
         * and Gl \in coKernel of J_\emptyset at level l*d
         * deg H = (l+1) * d - (n+1)
         */
        map< Vec<int64_t>, Vec< Mat<ZZ_p> >, vi64less > inclusion_matrix_ND_dict;

        /*
         * the dictionary where we will store the matrices that map W_u to the basis for middle cohomology
         * u -> M(u):W_u -> <cohomlogy basis>
         * where sum(u) = d
         * and u[i] > 0
         */
        map< Vec<int64_t>, Mat<ZZ_p>, vi64less> coKernels_ND_to_basis_dict;

        /*
         * the dictionary where we will store the reduction matrices from W_{u+v} -> W_u
         * v -> M(u): W_{u+v} -> W_u
         * where M(u) is stored as a polynomial
         */
        map< Vec<int64_t>, map<Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less> reduction_matrix_ND_dict;
         /*
         * same map as above but stores matrices already over ZZ to avoid unecessary conversions
         */
        map< Vec<int64_t>, map<Vec<int64_t>, Mat<ZZ>, vi64less> , vi64less> reduction_matrix_ND_ZZ_dict;
        /*
         * the dictionary where we will store the reduction matrices from W_{u+v} -> W_u
         * v -> M(u,x): W_{u+(x+1)*v} -> W_{u+x*v}
         * where M(u,x) is a stored as polynomial
         */
        map< Vec<int64_t>, map<Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less> reduction_matrix_ND_poly_dict;
        /*
         * same but they are stored over ZZ
         */
        map< Vec<int64_t>, map<Vec<int64_t>, Mat<ZZ>, vi64less> , vi64less> reduction_matrix_ND_poly_ZZ_dict;



        /*
         * FUNCTIONS
         */

        /*
         * constructors
         */

        de_Rham_non_degenerate_local(){};
        de_Rham_non_degenerate_local(int64_t p, int64_t precision, int64_t n, int64_t d, bool verbose, bool save_memory = true);
        de_Rham_non_degenerate_local(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<ZZ_p> fbar_vector, bool verbose, bool save_memory = true);
        de_Rham_non_degenerate_local(int64_t p, int64_t precision, map< Vec<int64_t>, ZZ_p, vi64less> fbar, bool verbose, bool save_memory = true);
        de_Rham_non_degenerate_local(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<zz_p> fbar_vector, bool verbose, bool save_memory = true);
        de_Rham_non_degenerate_local(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose, bool save_memory = true);
        de_Rham_non_degenerate_local(const char* filename);
        bool save(const char * filename);

        ~de_Rham_non_degenerate_local(){};

        void init_ND(int64_t p, int64_t precision, map< Vec<int64_t>, ZZ_p, vi64less> f, bool verbose, bool save_memory);
        

        /*
         * computes all the reduction and inclusion matrices
         */
        void compute_everything_ND(bool J=true, bool ND_ZZ=true);

        /*
         * computes the pair solve_ND[l] and adds it to the dictionary (if not computed already)
         * returns the pointer for the pair
         */
        pair< Vec<int64_t>, Mat<ZZ_p> >* get_solve_ND(const int64_t level);


        /*
         * returns the reduction matrix W_{u+v} to W_u
         */
        Mat<ZZ_p> get_reduction_matrix_ND(const Vec<int64_t> u, const Vec<int64_t> v);
         /*
         * retuns the reduction matrix (over ZZ) from W_{u+v} to W_u
         * that will be correct in ZZ_p
         */
        Mat<ZZ> get_reduction_matrix_ND_ZZ(const Vec<int64_t> u, const Vec<int64_t> v);

        /*
         * reduces a vector from W_{u+v} to W_u 
         * using the matrices over ZZ or ZZ_p
         */
        void reduce_vector_ND(Vec<ZZ_p> &result, const Vec<int64_t> u, const Vec<int64_t> v, const Vec<ZZ_p> G);
        void reduce_vector_ND_ZZ(Vec<ZZ> &result, const Vec<int64_t> u, const Vec<int64_t> v, const Vec<ZZ> G);

        /*
         * reduces a vector from W_{u+k*v) to W_u
         */
        void reduce_vector_ND_poly(Vec<ZZ_p> &result, const Vec<int64_t> u, const Vec<int64_t> v, const int64_t k, const Vec<ZZ_p> G);
        void reduce_vector_ND_poly_ZZ(Vec<ZZ> &result, const Vec<int64_t> u, const Vec<int64_t> v, const int64_t k, const Vec<ZZ> G);

        /*
         * returns the polynomial with matrix coefficients that 
         * reduces a vector from W_{u + x*v} to W_u
         */
        void get_ND_poly_flint(fmpz_mat_struct * result, const Vec<int64_t> u, const Vec<int64_t> v);
        /*
         * reduces a vector from W_{u+k*v) to W_u
         * this function is threadsafe!
         */
        void reduce_vector_ND_poly_flint(fmpz * result, fmpz_mat_struct * poly, const int64_t k, const fmpz * G, fmpz_t modulus);



        /*
         *  computes the reduction matrix W_{u+v} to W_u as polynomial in u0, ... un
         *  and automatically adds it to the dictionary
         */
        map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less>::const_iterator compute_reduction_matrix_ND(const Vec<int64_t> v);
        map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ>, vi64less> , vi64less>::const_iterator compute_reduction_matrix_ND_ZZ(const Vec<int64_t> v);
        /*
         * computes the reduction matrix W_{u+(x+1)v} to W_{u + x*v} to a 
         * polynomial in u0, ..., un and x
         * and automatically adds it to the dictionary
         */
        map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less>::const_iterator compute_reduction_matrix_ND_poly(const Vec<int64_t> v);
        map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ>, vi64less> , vi64less>::const_iterator compute_reduction_matrix_ND_poly_ZZ(const Vec<int64_t> v);


        /*
         * returns the inclusion matrix
         * from x^u (l-1)! G_l/f^l (x0...xn) to -> (l-1)! H / f^l (the factorial is implicit)
         * deg H = (l + 1) * d - (n + 1) 
         */
        Mat<ZZ_p>* get_inclusion_matrix_ND(const Vec<int64_t> u, const int64_t l);

        /*
         * compute the inclusions matrices from
         * x^u G_l/f^l (x0...xn) to -> H / f^l
         * deg H = (l + 1) * d - (n + 1)
         * and adds it automatically to the dictionary
         */
        void compute_inclusion_matrix_ND(const Vec<int64_t> u);

        /*
         * returns the matrix that maps W_u to the cohomology basis
         */
        Mat<ZZ_p>* get_coKernels_ND_to_basis(const Vec<int64_t> u);

        /*
         * computes coKernels_ND_to_basis[u] and adds it to the dictionary
         */
        void compute_coKernels_ND_to_basis(const Vec<int64_t> u);


        /*
         * computes the matrix map
         * (H0,...,Hn) -> H0 * x0 * df0 + ... + Hn * xn * dfn
         * Where Hi have degree = l - d
         */
        void matrix_ND(Mat<ZZ_p>& result, int64_t l);


        /*
         * Computes the coordinates of (m-1)! x^u \Omega / x0 ... xn f^m
         * in the cohomology basis
         */
        Vec<ZZ_p> monomial_to_basis_ND(Vec<int64_t> u);

        /*
         * Test Functions
         */

        /*
         * checks that a basis element x^w \Omega / f^k ~ x^w f^l \Omega / f^(l+k)
         * for l <= N
         */
        bool test_monomial_to_basis_ND(int64_t N);

        /*
         * tests the path independence while performing reductions
         */
        bool test_paths_ND(int64_t trials, int64_t paths);

};




bool isND(const map< Vec<int64_t>, zz_p, vi64less> &f);


/* 
 * returns a change of variables that will make f a non degenerate hypersurface
 * if the matrix returned is the zero matrix then we weren't able to find such change of variables in N tries
 */
Mat<zz_p> find_change_of_variables( map< Vec<int64_t>, zz_p, vi64less> f, int64_t N = 200 );

#endif
