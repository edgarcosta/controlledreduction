// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "tools.h"

#include <assert.h>
#include <stdint.h>

#include "vec_int64.h"


map< Vec<int64_t>, int64_t, vi64less> change_of_variables_monomial(const Vec<int64_t> u,const Vec<int64_t> v)
{
    int64_t n =  u.length();
    assert( n ==  v.length() );
    map< Vec<int64_t>, int64_t, vi64less > monomial_extended, H, Hnew;
    map< Vec<int64_t>, int64_t, vi64less >::const_iterator itH, it_monomial;
    map< Vec<int64_t>, int64_t, vi64less >::iterator itHnew;


    Vec<int64_t> zero, monomial, tmp_monomial;
    zero.SetLength(n+1);
    for(int64_t i = 0; i <= n; ++i)
        zero[i] = 0;

    // monomials_extended[i] = (xi - x v[i] )^u[i] = \sum xi^(u[i] - j) x^j (u[i] choose j) (-1)^j v[i]^j
    H[zero] = 1; 
    for(int64_t i = 0; i < n; ++i)
    {
        monomial_extended.clear();
        for(int64_t j = 0; j <= u[i] ; ++j)
        {
            tmp_monomial = zero;
            tmp_monomial[i] = u[i] - j;
            tmp_monomial[n] = j;
            monomial_extended[tmp_monomial] = binomial(u[i], j) * power_long( v[i],  j); 
        }
        
        Hnew.clear();
        for(itH = H.begin(); itH != H.end(); ++itH)
        {
            for( it_monomial = monomial_extended.begin(); it_monomial != monomial_extended.end(); ++it_monomial)
            {
                monomial = itH->first + it_monomial->first;
                itHnew = Hnew.find(monomial);

                if(itHnew == Hnew.end() )
                    Hnew[monomial] = itH->second * it_monomial->second;
                else
                    itHnew->second += itH->second *it_monomial->second; 
            }
            
        }

        H.swap(Hnew);
    }
    return H;
}
