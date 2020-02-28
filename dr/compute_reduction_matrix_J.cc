// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "dr.h"

/*
 * computes the reduction matrix from V_{u+v} to V_u as a polynomial in u0,...,un
 * and automatically adds it to the dictionary
 * Recall V_u = {  (m-1)! x^u G Omega / x0 ... xn  f^m}
 * G dense polynomial of degree d*n-n and x0...xn | x^u G
 * and we assume that u_i == 0 for every i such that v_i == 0
*/
map< Vec<int64_t>, Vec<Mat<ZZ_p> >, vi64less>::const_iterator de_Rham_local::compute_reduction_matrix_J(const Vec<int64_t> v)
{
    map< Vec<int64_t> ,  Vec<Mat<ZZ_p> >, vi64less>::const_iterator it;
    it = reduction_matrix_J_dict.find(v);
    if(it != reduction_matrix_J_dict.end())
    {
        return it;
    }
    else
    {
        if(verbose)
            cout<<"Computing the reduction matrix J for v = "<<v<<endl;

        Vec< Vec<int64_t> >* list_G;
        Vec< Vec<int64_t> >* list_F;
        map< Vec<int64_t>, int64_t, vi64less>* dict_G;
        map< Vec<int64_t>, int64_t, vi64less>* dict_w;


        Mat<ZZ_p>* solve_top;
        solve_top = &( ( get_solve_J( (d - 2)*(n + 1) + 1) )->second );


        list_G = &tuple_list[ d * n - n ];
        dict_G = &tuple_dict[ d * n - n ];
        dict_w = &tuple_dict[ (d * n - n) + d - (n+1) ]; //deg G + sum(v) - (n+1) = (d-2)*(n+1) + 1
        list_F = &tuple_list[ d * n - 2 * n ]; //deg G + d -(n+1) - (d-1) = (d-2)*(n+1) + 1 - ( d - 1)


        int64_t dim_Fi = list_F->length();

        Vec< Mat<ZZ_p> >* M;
        M = &(reduction_matrix_J_dict[v]);
        M->SetLength(n+2);

        for(int64_t i = 0; i < n + 2; i++)
        {
            (*M)[i].SetDims(list_G->length(),list_G->length());
        }


        #ifdef _OPENMP
        ZZ_pContext context;
        context.save();
        #endif

        #pragma omp parallel for
        for(int64_t coordinate_of_monomial = 0; coordinate_of_monomial < (int64_t) list_G->length() ; ++coordinate_of_monomial)
        {
            #ifdef _OPENMP
            context.restore();
            #endif
            Vec<int64_t> one_minus_ei;
            one_minus_ei.SetLength(n+1);
            for(int64_t i = 0; i<= n ; i++)
              one_minus_ei[i] = 1;
            Vec<int64_t> w;
            w = (*list_G)[coordinate_of_monomial] + v; // sum(w) = d * n - n + d
            bool boolean = true;
            for(int64_t i = 0; i <= n; i++)
            {
                if( w[i] > 0)
                    w[i]--;
                else
                {
                    boolean = false;
                    break;
                }
            }
            if(boolean) // sum(w) = d * n - 2n + d - 1
            {
                int64_t sum=0;
                for(int64_t i =0 ; i <= n ; ++i)
                    sum += w[i];
                assert(sum == d * n - 2*n + d - 1);
                int64_t coordinate_of_w = (*dict_w)[w];
                // F = solve_top.column(coordinate_of_w)
                Vec<ZZ_p> F;
                F.SetLength(solve_top->NumRows());
                for(int64_t i = 0; i < (int64_t) solve_top->NumRows(); ++i) {
                    F[i] = solve_top->get(i, coordinate_of_w);
                }
                for(int64_t i = 0; i <= n; i++) {
                    one_minus_ei[i]--;
                    for(int64_t coordinate_of_monomial_of_Fi = 0; coordinate_of_monomial_of_Fi < dim_Fi; ++coordinate_of_monomial_of_Fi) {
                        // x^(u+v) * G / x0 ... xn = \sum c_w x^u * x^w
                        // x^u * x^w = x^u \sum F_i di f ~ \sum di( x^u F_i )
                        // x^u * x^s * di f ~ (u[i]+s[i]) x^u * x^s / xi
                        // (u[i]+s[i]) * x^u * x^s x0 * ... * xi-1 * xi+1 * ... xn / x0 ... xn
                        int64_t row = (*dict_G)[ (*list_F)[coordinate_of_monomial_of_Fi] + one_minus_ei ];
                        (*M)[0][row][coordinate_of_monomial] += (*list_F)[coordinate_of_monomial_of_Fi][i] * F[i*dim_Fi + coordinate_of_monomial_of_Fi];
                        (*M)[i+1][row][coordinate_of_monomial] += F[i*dim_Fi + coordinate_of_monomial_of_Fi];
                    }
                    one_minus_ei[i]++;
                }
            }
        }
        return reduction_matrix_J_dict.find(v);
    }
}


