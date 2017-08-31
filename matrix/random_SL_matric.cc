// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include <assert.h>
#include <NTL/matrix.h>
#include <NTL/mat_lzz_p.h>

#include <cstdint>

using namespace NTL;

Mat<zz_p> random_SL_matrix(int64_t n)
{
    int64_t i, j, randomRow;
    Mat<zz_p> result, copy_result;
    Vec<zz_p> vector;
    zz_p det;
    copy_result.SetDims(n,n);
    
    bool invertible = false;

    while(!invertible)
    {
        result.SetDims(n, n);
        for( i = 0; i < n ; ++i)
        {
            for( j = 0; j < n ; j++ )
            {
                result[i][j] = random_zz_p();
                copy_result[i][j] = result[i][j];
            }
        }
        invertible = (n == (int64_t) gauss(copy_result) );
    }

    determinant(det, result);
    
    randomRow = RandomBnd(n);
    mul(result[ randomRow ], inv(det), result[ randomRow ] );
    determinant(det, result);
    assert( det == 1 );
    return result;
}

