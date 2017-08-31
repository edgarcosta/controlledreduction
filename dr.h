// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
//
// dr.h: header file for de_Rham_local class, which implements the de Rham structure to perform controlled_reductions

#ifndef DR_H
#define DR_H

#include "conv.h"
#include "solve_system.h"
#include "tools.h"
#include "matrix.h"
#include "vec_int64.h"
#include <NTL/ZZX.h>
#include <stdio.h>
#include <map>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;
using namespace NTL;

class de_Rham_local{
    /*
     * Some notation:
     * V_u = { (m-1)! x^u G Omega / x0 ... xn f^m}
     * with G a dense polynomial of degree d * n - n and x0 ... xn | x^u G
     */
    public:
        int64_t p;
        int64_t precision;
        bool save_memory;
        bool verbose;

        map< Vec<int64_t>, zz_p, vi64less> fbar;
        map< Vec<int64_t>, ZZ_p, vi64less> f;
        
        int64_t d;
        int64_t n;
        
        /*
         * tuple_list[l] = the list (n+1)-tuples that are a partition of l
         */
        Vec< Vec< Vec<int64_t> > > tuple_list;

        /*
         * tuple_dict[l] is the dictionary for reverse lookup of (n+1)-tuples that are a partition of l
         */
        Vec< map< Vec<int64_t>, int64_t, vi64less> > tuple_dict;

        /*
         * the Hilbert polynomial as a vector
         */
        Vec<int64_t> Hilbert_J;
        
        /*
         * the basis for the midle cohomology
         */
        Vec< Vec<int64_t> > coKernels_J_basis;
        /*
         * the dictionary for reverse lookup
         */
        map< Vec<int64_t>, int64_t, vi64less> coKernels_J_basis_dict;
        
        /*
         * the equivalent of solve_Ji as a dictionary i->solve_Ji
         * solve_Ji is the matrix that will tell us how to write a monomial of degree i in terms of 
         * dfi and coKernel elements
         */
        map< int64_t, pair< Vec<int64_t>, Mat<ZZ_p> > > solve_J;

        /*
         * the dictionary where we will store the inclusion matrices 
         * u -> inclusion_matrix_J(u) (sum(u) = n)
         * where  inclusion_matrix_J(u) is the matrix that maps
         * V_u - > {(n-1)! H Omega/f^n : deg H = d*n - n - 1}
         */
        map< Vec<int64_t>, Mat<ZZ_p>, vi64less> inclusion_matrix_J_dict;

        /*
         * the dictionary where we will store the final reduction matrices
         * pole -> matrix that reduces (pole -1)! G \Omega / f^pole to the cohomology basis elements
         */
        Vec< Mat<ZZ_p> >  final_reduction_matrix_J_dict;
        
        /*
         * the dictionary where we will store the reduction matrices as polynomials
         * v -> reduction_matrix(v)(u) \in ZZ_p[u0, ..., un] represented as a vector M
         * consntant term = M[0]
         * ui term = M[i+1]
         * where  reduction_matrix(v)(u) is the map from V_{u+v} to V_{u}
         */
        map< Vec<int64_t> ,  Vec< Mat<ZZ_p> > , vi64less>  reduction_matrix_J_dict;
        /*
         * same map as above but stores matrices already over ZZ to avoid unecessary conversions
         */
        map< Vec<int64_t> ,  Vec< Mat<ZZ> > , vi64less>  reduction_matrix_J_ZZ_dict;
        
        /*
         * FUNCTIONS
         */

        /*
         * constructors
         */
        de_Rham_local(){};
        de_Rham_local(int64_t p, int64_t precision, int64_t n, int64_t d, bool verbose, bool save_memory = true);
        de_Rham_local(int64_t p, int64_t precision, int64_t n, int64_t d, Vec< zz_p> fbar_vector, bool verbose, bool save_memory = true);
        de_Rham_local(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose , bool save_memory = true);
        de_Rham_local(const char * filename);
        virtual bool save(const char * filename);
        virtual ~de_Rham_local(){};
        
        void init(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose, bool save_memory);

        /*
         * computes all the reduction and inclusion matrices
         */
        void compute_everything_J();

        /*
         * computes the pair solve_Jl and adds it to the dictionary (if not compute already)
         * returns the pointer for the pair
         */
        pair< Vec<int64_t>, Mat<ZZ_p> >* get_solve_J(const int64_t level);
        

        /*
         * returns the reduction matrix from V_{u+v} to V_u
         */
        Mat<ZZ_p> get_reduction_matrix_J(const Vec<int64_t> u, const Vec<int64_t> v);

        /*
         * retuns the reduction matrix (over ZZ) from V_{u+v} to V_u
         * that will be correct in ZZ_p
         */
        Mat<ZZ>  get_reduction_matrix_J_ZZ(const Vec<int64_t> u, const Vec<int64_t> v);
    
        Vec<ZZ_p> reduce_vector_J(const Vec<int64_t> u, const Vec<int64_t>, const Vec<ZZ_p> G);

        void reduce_vector_J_ZZ(Vec<ZZ> &result, const Vec<int64_t> u, const Vec<int64_t> v, const Vec<ZZ> G);
        /*
         * reduces a vector from V_{u+kv) to V_u
         */
        void reduce_vector_J_poly(Vec<ZZ_p> &result, const Vec<int64_t> u, const Vec<int64_t> v, const int64_t k, const Vec<ZZ_p> G);
        void reduce_vector_J_poly_ZZ(Vec<ZZ> &result, const Vec<int64_t> u, const Vec<int64_t> v, const int64_t k, const Vec<ZZ> G);


        /*
         * returns the iterator for the reduction matrix from V_{u+v} to V_u as a polynomial in u0,...,un
         * if not computed adds it to the map
         */
        map< Vec<int64_t>, Vec<Mat<ZZ_p> >, vi64less>::const_iterator compute_reduction_matrix_J(const Vec<int64_t> v);
        map< Vec<int64_t>, Vec<Mat<ZZ> >, vi64less>::const_iterator compute_reduction_matrix_J_ZZ(const Vec<int64_t> v);

        
        /*
         * retuns a pointer to the inclusion_matrix_J(u)
         */
        Mat<ZZ_p>* get_inclusion_matrix_J(Vec<int64_t> u);

        /*
         * computes the inclusion_matrix_J(u) and adds it to the dictionary 
         */
        void compute_inclusion_matrix_J(Vec<int64_t> u);


        /*
         * computes all the final reduction matrices for 0 <= pole <= k
         * that reduce (pole -1)! G \Omega / f^pole to the cohomology basis elements
         * and adds them to the list.
         */
        void compute_final_reduction_matrix_J(int64_t k);

        /*
         * returns the final reduction matrix for
         * (k -1)! G \Omega / f^k to the cohomology basis elements
         */
        Mat<ZZ_p>* get_final_reduction_matrix_J(int64_t k);


        
        /*
         * computes the matrix map
         *  (H0,...,Hn) -> H0 * df0 + ... + Hn * dfn
         *  Where Hi have degree = l - (d-1).
         */
        void matrix_J(Mat<ZZ_p> &result, int64_t l);



        /*
         * Computes the coordinates of (m-1)! x^u \Omega / x0 ... xn f^m
         * in the cohomology basis
         */
        Vec<ZZ_p> monomial_to_basis_J(Vec<int64_t> u);

        /*
         * Test Functions
         */

        /*
         * checks that a basis element x^w \Omega / f^k ~ x^w f^l \Omega / f^(l+k)
         * for l <= N
         */
        void test_monomial_to_basis_J(int64_t N);

        /*
         * tests the path independence while performing reductions
         */
        void test_paths_J(int64_t trials, int64_t paths);

};





bool isSmooth( const map< Vec<int64_t>, zz_p, vi64less> &f );
        

#endif
