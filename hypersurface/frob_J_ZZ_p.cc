// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

/*
 * returns the coordinates of the p-adic approximation of frob(e_i)
 * using N terms, using reduce_vector_J_poly_ZZ
 */

Vec<ZZ_p> hypersurface::frob_J_ZZ_p(const int64_t coordinate, const int64_t N, const int64_t loop)
{
    assert(N > 0);
    int64_t bin, e, i, j, m, sum, end, val;
    ZZ_p fact;
    //ZZ factZZ;
    ZZ tmp, remainder;
    int64_t Djm;
    Vec<int64_t> ei, monomial, ones, src, dest, v, y;
    Vec<ZZ_p> F, G, G0, H_last, *H_temp;
    Mat<ZZ_p> Mred;
    map< Vec<int64_t>, int64_t, vi64less > *dict_V, *dict_G0;
    map< Vec<int64_t>, Vec<ZZ_p>, vi64less > H, Hnew;
    map< Vec<int64_t>, Vec<ZZ_p>, vi64less >::const_iterator Hit;
    map< Vec<int64_t>, Vec<ZZ_p>, vi64less >::iterator Hnewit;
    map< Vec<int64_t>, ZZ_p, vi64less>::const_iterator it;
    
    NTL::set(fact);

    v.SetLength(n+1);

    ones.SetLength(n+1);
    for( i = 0; i <= n; i++)
        ones[i] = 1;

    ei = dR->coKernels_J_basis[coordinate];

    compute_fpow(N-1);

    sum = 0;
    for(i = 0; i <= n; i++)
        sum += ei[i];

    m = (sum + n + 1)/d;

    dict_V = &dR->tuple_dict[ d * n - n ];
    // H = map of differentials stored in V_{ tau( p * u ) }
    // Stored as a u -> coordinates in V_{ tau( p * u ) }

    H_last.SetLength(binomial(d*n-1,n)); //degree = d*n-n+n-(n+1) = d*n - (n+1)

    for( e = m + N - 1; e > n/p; e--)
    {
        Hnew.clear();

        if(verbose)
            cout<<"e = "<<e<<endl<<"n\tadding new terms\n";

        // add new terms with pole order p*e = p*(m+j)

        if( e >= m )
        {
            j = e - m;
            Djm = 0;
            // Djm = sum( [ binomial(m + i -1, m -1) * binomial(i , j) * (-1) ** j for i in range(j,N) ] )
            for( i = j ; i < N ; i++)
            {
                bin = binomial(m + i - 1, m - 1) * binomial(i,j);
                if(j%2 == 0)
                    Djm += bin;
                else
                    Djm -= bin;
            }
            for( it = fpow.get(j).begin(); it != fpow.get(j).end(); ++it)
            {
                monomial = it->first + ei + ones;
                H_temp = &H[monomial];
                H_temp->SetLength( dict_V->size() );
                (*H_temp)[ (*dict_V)[ p * monomial - tweak( p * monomial, d * n - n) ] ] += Djm * it->second * fact;
                /*
                if( H.find(monomial) == H.end() )
                {
                    H_temp = &H[monomial];
                    H_temp->SetLength( dict_V->size() );
                    for(i = 0; i < dict_V->size(); i++)
                        clear( (*H_temp)[i] );
                }
                H[monomial][ (*dict_V)[ p * monomial - tweak( p * monomial, d * n - n) ] ] += Djm * it->second * fact;
                */
            }
        }

        /*
         * before
         * fact = (p * (m + N - 1) - 1)!/(p * e - 1)!
         * after
         * fact = (p * (m + N - 1) - 1)! / max(p * (e-1) - 1,0)!
         */
        //assert( fact == factorial<ZZ_p>( p * (m + N - 1) - 1, p*e ) );
        end = (e > 1) ? p : p - 1;

        for( i = 0; i < end ; i++)
            fact *= p*e - i - 1;
        /*
        if (e > 1)
            assert( fact == factorial<ZZ_p>(p * (m + N - 1) - 1, p*(e-1) ) );
        else
            assert( fact == factorial<ZZ_p>(p * (m + N - 1) - 1) );
        */
        // reduce all terms from pole order d * p * e   to  d * p * (e-1)

        for( Hit = H.begin(); Hit != H.end(); Hit++)
        {
            if( verbose )
                cout<<"\treducing "<< Hit->first <<endl;

            G = Hit->second;
            // choosing the main reduction direction
            sum = 0;

            for( i = 0; i <= n; i++)
            {
                if(Hit->first[i] > 0)
                {
                    v[i] = 1;
                    sum++;
                }
                else
                    v[i] = 0;
            }

            while( sum < d )
            {
                for(i = 0; i <= n; i++)
                {
                    if(Hit->first[i] > v[i])
                    {
                        v[i]++;
                        sum++;
                        break;
                    }
                }
            }

            end = (e > n/p + 1) ? p : p - n%p;

            switch(loop)
            {
                case 0:
                    i = 0;
                    if( end > d * n - n )
                    {
                        // tweak( p
                        dest = p * (Hit->first - v) + (p - end) * v + tweak( (d*n - n) * v, (d*n - n) );
                        assert( dest ==  tweak(p * (Hit->first - v) + (p - end)*v + (d*n - n) * v , d*n - n ) );
                        // dest == tweak( p * Hit->first - (end - (d*n - n))*v, d*n -n)
                        dR->reduce_vector_J_poly(G, dest, v, end - (d * n - n), G);
                        i = end - (d*n - n);
                    }
                    for( ; i < end; i++)
                    {
                        src = tweak(p * Hit->first - i * v, d * n - n);
                        dest = tweak(p * Hit->first - (i + 1) * v, d * n - n);
                        y = src - dest;
                        G = dR->reduce_vector_J( dest,  y, G);
                    }
                    break;
                case 1:
                    for(i = 0; i < end; i++)
                    {
                        src = tweak(p * Hit->first - i * v, d * n - n);
                        dest = tweak(p * Hit->first - (i + 1) * v, d * n - n);
                        y = src - dest;
                        G = dR->reduce_vector_J(dest, y, G);
                    }
                    break;
                case 2:
                    for(i = 0; i < end; i++)
                    {
                        src = tweak(p * Hit->first - i * v, d * n - n);
                        dest = tweak(p * Hit->first - (i + 1) * v, d * n - n);
                        y = src - dest;
                        Mred = dR->get_reduction_matrix_J(dest, y);
                        mul(G,Mred,G);
                    }
                    break;
            }

            if(e > n/p + 1)
            {
                y = Hit->first - v;
                Hnewit =  Hnew.find(y);
                if( Hnewit == Hnew.end() )
                    Hnew[y] = G;
                else
                    Hnewit->second += G;
            }
            else
            {
                sum = 0;
                for( i = 0; i <= n; i++)
                    sum += dest[i];

                mul(G, *(dR->get_inclusion_matrix_J(dest)), G);

                H_last += G;
            }
        }
        //reduced all terms from p * e to p * (e - 1)
        
        H.swap(Hnew);
    }

    // H_last has the poles with order n
    // final reduction step for H_last
    mul(F, *(dR->get_final_reduction_matrix_J(n)), H_last);
    for(e = n/p; e != 0 && e >= m; e-- )
    {
        if(verbose)
            cout<<"e = "<<e<<endl<<"\tadding new terms\n";

        dict_G0 = &dR->tuple_dict[ d * p * e - (n + 1) ];
        clear(G0);
        G0.SetLength( dict_G0->size() );

        j = e - m;
        Djm = 0;
        for( i = j ; i < N ; i++)
        {
            bin = binomial(m + i - 1, m - 1) * binomial(i,j);
            if(j%2 == 0)
                Djm += bin;
            else
                Djm -= bin;
        }

        for( it = fpow[j].begin(); it != fpow[j].end(); it++)
        {
            monomial = p * ( it->first + ei + ones ) - ones;
            G0[ (*dict_G0)[ monomial ] ] += fact * Djm * it->second;
        }
        F += (*(dR->get_final_reduction_matrix_J(p * e))) * G0;

        /*
         * before
         * fact = (p * (m + N - 1))!/(p * e)!
         * after
         * fact = (p * (m + N - 1) - 1)! / max(p * (e-1) - 1,0)!
         */
        
        end = (e > 1) ? p : p - 1;

        for( i = 0; i < end ; i++)
            fact *= p*e - i - 1;
    }

    assert( fact == factorial<ZZ_p>( p*(m + N - 1) - 1 ) );
    val = valuation_of_factorial(p*(m + N - 1) - 1, p);
    DivRem(tmp, remainder,  factorial<ZZ>( p*(m + N - 1) - 1 ) , power_ZZ(p,val));
    assert(IsZero(remainder));
    inv(fact, to_ZZ_p(tmp) );
    F *= fact;
    
    if( n - 1 >= val)
        F *= to_ZZ_p(power_ZZ(p, n - 1 - val));
    else
    {
        ZZ ppower = power_ZZ(p,val - ( n - 1) );
        for(i = 0; i < (int64_t) F.length(); i++)
        {
            //Not sure if DivRem(F[i], remainder, rep(F[i]), ppower) is safe
            DivRem(tmp, remainder, rep(F[i]), ppower);
            assert( IsZero(remainder) );
            F[i] = to_ZZ_p(tmp);
        }
    }

    return F;
}
