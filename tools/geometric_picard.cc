// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "conv.h"
#include "tools.h"

#include <assert.h>

#include <flint/arith.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_mat.h>
#include <flint/fmpz_poly_factor.h>
#include <NTL/vector.h>
#include <NTL/ZZ.h>

#include <cstdint>


using namespace NTL;


int64_t number_of_roots_of_unity(fmpz_poly_mat_t cyclo, fmpz_poly_t f, int64_t p)
{
    int64_t degree = fmpz_poly_degree(f);
    int64_t result;
    bool resultQ = false;
    fmpz_t ppower;
    fmpz_t rem;
    fmpz* coeff;
    fmpz_init(ppower);
    fmpz_init(rem);
    fmpz_set_ui(ppower,p);
    //Performing change of variables T->T/p 
    for(int64_t i = 1; i <= degree; ++i)
    {
        //ppower = p^i
        coeff = fmpz_poly_get_coeff_ptr(f, i);
        fmpz_fdiv_qr(coeff, rem, coeff, ppower);
        if( !fmpz_is_zero(rem) )
        {
            //polynomial is not integral after the change of variables
            result = 0;
            resultQ = true;
            break;
        }
        fmpz_mul_ui(ppower, ppower, p);
    }


    int64_t k = 1;
    while( ((double)k)/(log(k)/log(2) + 1) <= degree && !resultQ )
    {
        if( fmpz_poly_equal(f, fmpz_poly_mat_entry(cyclo, 0, k)) )
        {
            result = (int64_t) degree;
            resultQ = true;
            break;
        }
        k++;
    }

    if ( !resultQ )
        result = 0;
    
    fmpz_clear(ppower);
    fmpz_clear(rem);
    return result;
}

int64_t geometric_picard(const Vec<ZZ> &charpoly_frob, const int64_t p, const int64_t dimension)
{
    fmpz_poly_t p_flint;
    fmpz_poly_init( p_flint ) ;
    fmpz_t tmp;
    fmpz_init(tmp);
    int64_t i;
    for(i = 0; i < (int64_t) charpoly_frob.length(); ++i)
    {
        conv(tmp, charpoly_frob[charpoly_frob.length() - 1 - i]);
        fmpz_poly_set_coeff_fmpz(p_flint , i , tmp) ;
    }
    
    fmpz_poly_factor_t fac;
    fmpz_poly_factor_init(fac);
    fmpz_poly_factor_zassenhaus(fac, p_flint);
    long max_degree = 0;
    for( i = 0; i < (int64_t) fac->num; ++i)
        max_degree = std::max(max_degree, fmpz_poly_degree(fac->p + i));
    fmpz_poly_mat_t cyclo;
    int64_t N = 1;
    while( ((double)N)/(log(N)/log(2) + 1) <= max_degree )
    {
        N++;
    }

    fmpz_poly_mat_init(cyclo, 1, N);
    for(i = 0; i < N; ++i)
       fmpz_poly_cyclotomic(fmpz_poly_mat_entry(cyclo,0,i),i);

    int64_t count = 0;
    for(i = 0; i < (int64_t) fac->num; ++i)
    {
        count += fac->exp[i] * number_of_roots_of_unity(cyclo, fac->p + i, p);
    }
    fmpz_poly_mat_clear(cyclo);
    fmpz_poly_factor_clear(fac);
    fmpz_clear(tmp);
    fmpz_poly_clear(p_flint);
    if (dimension%2 == 0)
        count++;
    assert(dimension != 2 || (int64_t) charpoly_frob.length()%2 == count %2);
    return count;
}

