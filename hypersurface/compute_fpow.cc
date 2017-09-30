// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"
/*
 * assures the array fpow 
 * (where we store the powers of f)
 * has enough powers of f
 */

void hypersurface::compute_fpow(int64_t power)
{
    if( power >= (int64_t) fpow.length() )
    {
        if(verbose)
            cout<< "compute_fpow("<<power<<")"<<endl;
        int64_t k;
        //avoiding reallocation, as maps are not  "reallocatable"
        Vec< map< Vec<int64_t>, ZZ_p, vi64less> > new_fpow;
        new_fpow.SetMaxLength(power + 1);
        new_fpow.SetLength( fpow.length());
        for( k = 0; k < (int64_t) fpow.length(); k++)
            new_fpow[k].insert(fpow[k].begin(), fpow[k].end()); //= fpow[k];

        //and now fpow has enough space to grow
        swap(fpow, new_fpow);

        map< Vec<int64_t>, ZZ_p, vi64less>::const_iterator itfkminus1, itf;
        map< Vec<int64_t>, ZZ_p, vi64less>::iterator itY;
        Vec<int64_t> u;

        for( k =  (int64_t) fpow.length() ; k <= power ; k++ )
        {
            map< Vec<int64_t>, ZZ_p, vi64less> Y;

            for( itfkminus1 = fpow.get(k - 1).begin() ; itfkminus1 != fpow.get(k - 1).end() ; itfkminus1++ )
            {
                for( itf = dR->f.begin() ; itf != dR->f.end() ; itf++ )
                {
                    u = itfkminus1->first + itf->first;
                    itY = Y.find(u);
                    if( itY == Y.end() )
                    {
                        Y[u] = itfkminus1->second * itf->second;
                    }
                    else
                    {
                        itY->second += itfkminus1->second * itf->second;
                    }
                }
            }

            append(fpow, Y);
        }
        assert((int64_t) fpow.length() == power + 1);
    }
}
