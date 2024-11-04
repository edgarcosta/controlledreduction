// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"


/*
 * computes the matrix that reduces
 * (k-1)! G \Omega / f^k to the cohomology basis elements
 */
void de_Rham_local::compute_final_reduction_matrix_J(int64_t k)
{
    assert( k > 0);
    assert( k <= n );
    int64_t degree_of_G, i, j;
    degree_of_G = k * d - (n + 1);
    Vec< Vec<int64_t> >* list_G;

    Mat<ZZ_p>* M;
    M = &final_reduction_matrix_J_dict[k-1];
    list_G = &tuple_list[ degree_of_G ];
    if( k == 1 )
    {
        M->SetDims( coKernels_J_basis.length() , list_G->length() );
        for( i = 0 ; i < (int64_t) coKernels_J_basis.length() ; i++)
        {
            for( j = 0; j < (int64_t) list_G->length(); j++)
            {
                if( (*list_G)[j] == coKernels_J_basis[i] )
                    NTL::set((*M)[i][j]);
            }
        }
    }
    else
    {
        int64_t fact, dim_Fi, coordinate_of_monomial, coordinate_of_monomial_of_Fi;
        Vec<ZZ_p> H;
        Vec<ZZ_p> F;
        Vec< Vec<int64_t> >* list_F;
        Vec< Vec<int64_t> >* list_H;
        map< Vec<int64_t>, int64_t, vi64less >* dict_H;
        Mat<ZZ_p>* M_low;

        Mat<ZZ_p>* U;
        Vec<int64_t>* B;
        pair< Vec<int64_t>, Mat<ZZ_p> >* pairBU;
        
        Mat<ZZ_p> Mtranspose;

        Mtranspose.SetDims(list_G->length(), coKernels_J_basis.length() );

        list_F = &tuple_list[degree_of_G - (d - 1)];
        list_H = &tuple_list[degree_of_G - d];
        dict_H = &tuple_dict[degree_of_G - d];

        H.SetLength(list_H->length());

        dim_Fi = list_F->length();
        
        pairBU = get_solve_J(degree_of_G);
        B = &(pairBU->first);
        U = &(pairBU->second);

        F.SetLength(U->NumRows());

        M_low = &final_reduction_matrix_J_dict[k-2];

        // (k - 1)!
        fact = 1;
        for( i = 1 ; i < k ; i++ )
        {
            fact *= i;
        }

        for( coordinate_of_monomial = 0; coordinate_of_monomial < (int64_t) list_G->length(); coordinate_of_monomial++)
        {
            for( i = 0 ; i < (int64_t) H.length(); i++)
                H[i] = 0;

            for( i = 0 ; i < (int64_t) F.length(); i++)
                F[i] = U->get(i, coordinate_of_monomial);

            for( i = 0 ; i <= n; i++)
            {
                for( coordinate_of_monomial_of_Fi = 0; coordinate_of_monomial_of_Fi < dim_Fi; coordinate_of_monomial_of_Fi++ )
                {
                    /*
                     * (k-1)! G = (k-1)! \sum Fi di f
                     * (k-1)! G \Omega / f^k ~ 
                     * (k-2)! \sum di Fi / f^(k - 1)= (k-2)! H / f^(k-1)
                     *  +
                     *  (k-1)! * basis_elements
                     */
                    if( (*list_F)[ coordinate_of_monomial_of_Fi ][i] > 0)
                        H[ (*dict_H)[ diff(  (*list_F)[ coordinate_of_monomial_of_Fi ], i) ] ] += (*list_F)[ coordinate_of_monomial_of_Fi ][i] * F[ i * dim_Fi + coordinate_of_monomial_of_Fi ];
                }
            }
            Mtranspose[ coordinate_of_monomial ] = (*M_low) * H;

            for( i = 0 ; i < (int64_t) B->length(); i++ )
            {
                Mtranspose[ coordinate_of_monomial ][ coKernels_J_basis_dict[ (*list_G)[ (*B)[i] ] ] ] += fact * F[ (n + 1) * dim_Fi + i];
            }
        }
        transpose(*M, Mtranspose);
    }
}

