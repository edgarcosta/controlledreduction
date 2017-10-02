// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "matrix.h"

#include <assert.h>
#include <stdint.h>

#include <NTL/matrix.h>
#include <NTL/ZZ.h>

using namespace NTL;

Vec<ZZ> charpoly_frob(const Mat<ZZ> M,  Vec<int64_t> prec, const int64_t p, const int64_t dimension)
{
    Vec<ZZ> cp = charpoly(M);
    Vec<ZZ> mod;
    ZZ rem;
    int signal = 1;
    int64_t degree = cp.length() - 1;

    assert(prec.length() == cp.length() );

    mod.SetLength( prec.length() );

    for (int64_t i = 0; i <= degree; ++i) {
        mod[i] =power_ZZ(p, prec[i]);
        cp[i] = cp[i] % mod[i];
    }
    cp[degree] = 1;
    if( dimension % 2 == 0 )
    {
        // figure out the sign if dimension is even, i.e., weight is even
        // for odd case the sign is 1
        for (int64_t i = 1; i <= degree/2; ++i) {
            ZZ p_power = power_ZZ(p, std::min(prec[i], prec[degree - i] + ((degree-2*i)*dimension) /2));
            if( cp[i] % p_power !=0 &&  cp[degree-i] % p_power != 0) {
                if (0 == (cp[i] + cp[degree - i] * power_ZZ(p, ((degree-2*i)*dimension) /2)) %  p_power )
                    signal = -1;
                else {
                    signal = 1;
                    assert(0 == (cp[i] - cp[degree - i] * power_ZZ(p, ((degree-2*i)*dimension) /2)) %  p_power );
                }
                break;
            }
        }
    }
    cp[0] = signal *  power_ZZ(p, (degree * dimension) /2));
    //apply the symmetry
    for (int64_t i = 0; i <= degree/2; ++i) {
        if( prec[i] >=  prec[degree - i] + ((degree-2*i)*dimension) /2 )
        {
            prec[degree - i] = prec[i] - ((degree-2*i)*dimension) /2 ;
            mod[degree - i] = power_ZZ(p, prec[degree - i]);
            DivRem(cp[degree - i], rem, signal * cp[i],  power_ZZ(p, ((degree-2*i)*dimension) /2) );
            assert(IsZero(rem));
            cp[degree - i] = cp[degree - i] % mod[degree -i];
        }
        else
        {
            prec[i] = prec[degree - i] + ((degree-2*i)*dimension) /2;
            mod[i] =  power_ZZ(p, prec[i]);
            cp[i] = signal * cp[degree - i] *  power_ZZ(p, ((degree-2*i)*dimension) /2);
            cp[i] = cp[i] % mod[i];
        }
    }
    //calculate the i-th power sum of the roots and correct cp allong the way
    Vec<ZZ> s,e;
    e.SetLength(degree + 1);
    for(int64_t k = 0; k <= degree; k++) {
        e[k] = (k%2 == 0)? cp[degree - k] : -cp[degree - k];
        if(k > 0)
            assert( (log(k)/log(p) + prec[degree-k] > log(2*degree)/log(p) + 0.5*dimension*k) );
    }
    ZZ sum;
    ZZ pN;
    s.SetLength(degree+1);
    for ( int64_t k = 1; k <= degree; k++) {
        sum = 0;
        for ( int64_t i = 1; i < k ; ++i)
            sum += ((i%2 == 1)? 1 : -1) * e[k-i] * s[i];
        s[k] = ( (k%2==0)? 1 : -1 ) * (sum - k * e[k] );
        pN = k * mod[degree - k];
        s[k] = s[k] % pN;
        if( sqr(s[k]) > degree*degree*power_ZZ(p, dimension*k) )
            s[k] = - ( (-s[k]) % pN );

        DivRem(e[k], rem, sum + ( (k%2 == 0)? -s[k] : s[k] ) , ZZ( k ));
        assert(IsZero(rem));
        cp[degree - k] = (k%2 == 0)? e[k]: -e[k];
    }
    return cp;
}
