// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
//
// hypersurface_bs.cpp: routines for the class hypersurface hypersurface_non_degenerate

#include "hypersurface_nd.h"
#include "timing.h"

#ifdef _OPENMP
# include <omp.h>
#endif

// non degenerate
//
hypersurface_non_degenerate::hypersurface_non_degenerate(int64_t p, int64_t precision, int64_t n, int64_t d, bool verbose) {
    dR = dR_ND = make_shared<de_Rham_non_degenerate_local>(p, precision, n, d, verbose);
    init_after_dR();
}


hypersurface_non_degenerate::hypersurface_non_degenerate(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<ZZ_p> f_vector, bool verbose) {
    dR = dR_ND = make_shared<de_Rham_non_degenerate_local>(p, precision, n, d, f_vector, verbose);
    init_after_dR();
}


hypersurface_non_degenerate::hypersurface_non_degenerate(int64_t p, int64_t precision, map< Vec<int64_t>, ZZ_p, vi64less> f, bool verbose) {
    dR = dR_ND = make_shared<de_Rham_non_degenerate_local>(p, precision, f, verbose);
    init_after_dR();
}

hypersurface_non_degenerate::hypersurface_non_degenerate(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<zz_p> fbar_vector, bool verbose) {
    dR = dR_ND = make_shared<de_Rham_non_degenerate_local>(p, precision, n, d, fbar_vector, verbose);
    init_after_dR();
}


hypersurface_non_degenerate::hypersurface_non_degenerate(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose) {
    dR = dR_ND = make_shared<de_Rham_non_degenerate_local>(p, precision, fbar, verbose);
    init_after_dR();
}


hypersurface_non_degenerate::hypersurface_non_degenerate(const char * filename) {
    dR = dR_ND = make_shared<de_Rham_non_degenerate_local>(filename);
    init_after_dR();
}

bool hypersurface_non_degenerate::save(const char * filename)
{
    return (dR_ND->save)(filename);
}

Vec<ZZ_p> hypersurface_non_degenerate::frob_ND(const int64_t coordinate, const int64_t N)
{
    //return frob_ND_ZZ(coordinate, N, 0);
    return frob_ND_flint(coordinate,N);
}

Vec<ZZ_p> hypersurface_non_degenerate::frob_ND_ZZ(const int64_t coordinate, const int64_t N, const int64_t loop)
{
    assert(N > 0);
    int64_t e, i, j, m, sum, end, val, dpowern;
    ZZ_p fact;
    //ZZ factZZ;
    ZZ tmp, remainder;
    ZZ bin, Djm;
    Vec<int64_t> ei, monomial, ones, src, dest, v, final_dest;
    Vec<ZZ> G, Gswap;
    Vec<ZZ_p> F, G0, G_ZZ_p, H_last, *H_temp;
    Mat<ZZ> Mred;
    map< Vec<int64_t>, Vec<ZZ_p>, vi64less > H, Hnew;
    map< Vec<int64_t>, Vec<ZZ_p>, vi64less >::const_iterator Hit;
    map< Vec<int64_t>, Vec<ZZ_p>, vi64less >::iterator Hnewit;
    map< Vec<int64_t>, ZZ_p, vi64less>::const_iterator it;

    set(fact);

    v.SetLength(n+1);

    ones.SetLength(n+1);
    for( i = 0; i <= n; i++)
        ones[i] = 1;

    ei = dR_ND->coKernels_J_basis[coordinate];

    F.SetLength(dR_ND->coKernels_J_basis.length());

    compute_fpow(N-1);

    if( loop == 0)
    {
        Vec<long> fpow_length;
        fpow_length.SetLength(N);
        for( i = 0; i < N ; i++ )
        {
            fpow_length[i] = fpow.get(i).size();
        }
        //cout << fpow_length << endl;
    }

    sum = 0;
    for(i = 0; i <= n; i++)
        sum += ei[i];

    m = (sum + n + 1)/d;

    dpowern = dR_ND->coKernels_ND_basis.length();

    H_last.SetLength(dpowern);

    for( e = m + N - 1; e != 0; e--)
    {
        Hnew.clear();

        if(verbose)
            cout<<"e = "<<e<<endl<<"\tadding new terms\n";

        if( e >= m)
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
                H_temp->SetLength( dpowern );
                (*H_temp)[0] += conv<ZZ_p>(Djm * rep(it->second) * rep(fact) );
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
            {
                cout<<"\treducing "<< Hit->first <<endl;
            }


            G = conv< Vec<ZZ> >(Hit->second);
            // choosing the main reduction direction

            if( e > 1)
            {
                sum = 0;
                for(i = 0; i<= n; i++)
                    v[i] = 0;

                while( sum < d)
                {
                    for( i = 0; i <= n; i++)
                    {
                       if( Hit->first[i] > v[i] + 1 )
                       {
                           v[i]++;
                           sum++;
                           break;
                       }
                    }
                }
            }
            else
                v = Hit->first;

            dest = p*Hit->first;

            end = (e > 1)? p: p - 1;

            switch(loop)
            {
                case 0:
                    dest =  p * Hit->first - end * v;
                    dR_ND->reduce_vector_ND_poly_ZZ(G, dest, v, end, G);
                    break;
                case 1:
                    for( i = 0; i < end; i++)
                    {
                        dest = dest - v;
                        dR_ND->reduce_vector_ND_ZZ(Gswap, dest, v, G);
                        for( j = 0; j < (int64_t) G.length(); j++)
                            rem(G[j], Gswap[j], ZZ_p::modulus() );
                    }
                    break;
                case 2:
                    for( i = 0; i < end; i++)
                    {
                        dest = dest - v;
                        Mred = dR_ND->get_reduction_matrix_ND_ZZ(dest, v);
                        mul(G, Mred , G);
                    }
                    break;
            }
            G_ZZ_p = conv< Vec<ZZ_p> >(G);
            if( ! IsZero(G_ZZ_p) )
            {
                if( e > 1 )
                {
                    final_dest = Hit->first - v;
                    assert( p*final_dest == dest);
                    Hnewit = Hnew.find(final_dest);
                    if( Hnewit == Hnew.end() )
                        Hnew[final_dest] = G_ZZ_p;
                    else
                        Hnewit->second += G_ZZ_p;
                }
                else
                {
                    assert(dest == Hit->first);
                    F += (*(dR_ND->get_coKernels_ND_to_basis(dest))) * G_ZZ_p;
                }
            }
        }
        // reduced all terms from p*e to p*(e-1) and if e = 1 to basis
        H.swap(Hnew);
    }
    //Perform factorial division
    //assert( fact == factorial<ZZ_p>( p*(m + N - 1) - 1) );
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

Vec<ZZ_p> hypersurface_non_degenerate::frob_ND_ZZ_p(const int64_t coordinate, const int64_t N, const int64_t loop)
{
    assert(N > 0);
    int64_t bin, e, i, j, m, sum, end, val, dpowern;
    ZZ_p fact;
    //ZZ factZZ;
    ZZ tmp, remainder;
    int64_t Djm;
    Vec<int64_t> ei, monomial, ones, src, dest, v, final_dest;
    Vec<ZZ_p> F, G, G0, H_last, *H_temp;
    Mat<ZZ_p> Mred;
    map< Vec<int64_t>, Vec<ZZ_p>, vi64less > H, Hnew;
    map< Vec<int64_t>, Vec<ZZ_p>, vi64less >::const_iterator Hit;
    map< Vec<int64_t>, Vec<ZZ_p>, vi64less >::iterator Hnewit;
    map< Vec<int64_t>, ZZ_p, vi64less>::const_iterator it;

    set(fact);

    v.SetLength(n+1);

    ones.SetLength(n+1);
    for( i = 0; i <= n; i++)
        ones[i] = 1;

    ei = dR_ND->coKernels_J_basis[coordinate];

    F.SetLength(dR_ND->coKernels_J_basis.length());

    compute_fpow(N-1);

    sum = 0;
    for(i = 0; i <= n; i++)
        sum += ei[i];

    m = (sum + n + 1)/d;

    dpowern = dR_ND->coKernels_ND_basis.length();

    H_last.SetLength(dpowern);

    for( e = m + N - 1; e != 0; e--)
    {
        Hnew.clear();

        if(verbose)
            cout<<"e = "<<e<<endl<<"\tadding new terms\n";

        if( e >= m)
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
                H_temp->SetLength( dpowern );
                (*H_temp)[0] += Djm * it->second * fact;
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
            {
                cout<<"\treducing "<< Hit->first <<endl;
            }

            G = Hit->second;
            // choosing the main reduction direction

            if( e > 1)
            {
                sum = 0;
                for(i = 0; i<= n; i++)
                    v[i] = 0;

                while( sum < d)
                {
                    for( i = 0; i <= n; i++)
                    {
                       if( Hit->first[i] > v[i] + 1 )
                       {
                           v[i]++;
                           sum++;
                           break;
                       }
                    }
                }
            }
            else
                v = Hit->first;

            dest = p*Hit->first;

            //end = (e > 1)? p: p - 1;

            switch(loop)
            {
                case 0:
                    dest =  p * Hit->first - end * v;
                    dR_ND->reduce_vector_ND_poly(G, dest, v, end, G);
                    break;
                case 1:
                    for( i = 0; i < end; i++)
                    {
                        dest = dest - v;
                        Mred = dR_ND->get_reduction_matrix_ND(dest, v);
                        mul(G, Mred , G);
                    }
                    break;
            }
            if(! IsZero(G) )
            {
                if( e > 1 )
                {
                    final_dest = Hit->first - v;
                    assert( p*final_dest == dest);
                    Hnewit = Hnew.find(final_dest);
                    if( Hnewit == Hnew.end() )
                        Hnew[final_dest] = G;
                    else
                        Hnewit->second += G;
                }
                else
                {
                    assert(dest == Hit->first);
                    F += (*(dR_ND->get_coKernels_ND_to_basis(dest))) * G;
                }
            }
        }
        // reduced all terms from p*e to p*(e-1) and if e = 1 to basis
        H.swap(Hnew);
    }
    //Perform factorial division
    //assert( fact == factorial<ZZ_p>( p*(m + N - 1) - 1) );
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

Vec<ZZ_p> hypersurface_non_degenerate::frob_ND_flint(const int64_t coordinate, const int64_t N)
{
    assert(N > 0);
    int64_t e, i, j, m, sum, end, val, dpowern;
    ZZ_p fact;
    //ZZ factZZ;
    ZZ tmp, remainder;
    ZZ bin, Djm;
    Vec<int64_t> ei, monomial, ones, src, dest, u, v, final_dest, shift;
    Vec<ZZ> G;
    Vec<ZZ_p> F, G0, G_ZZ_p, H_last;
    Mat<ZZ> Mred;
    map< Vec<int64_t>, Vec<ZZ_p>, vi64less > H, Hnew;
    map< Vec<int64_t>, Vec<ZZ_p>, vi64less >::const_iterator Hit;
    map< Vec<int64_t>, Vec<ZZ_p>, vi64less >::iterator Hnewit;
    map< Vec<int64_t>, ZZ_p, vi64less>::const_iterator it;

    set(fact);

    v.SetLength(n+1);
    shift.SetLength(n+1);

    ones.SetLength(n+1);
    for( i = 0; i <= n; i++)
        ones[i] = 1;

    ei = dR_ND->coKernels_J_basis[coordinate];

    F.SetLength(dR_ND->coKernels_J_basis.length());

    compute_fpow(N-1);


    sum = 0;
    for(i = 0; i <= n; i++)
        sum += ei[i];

    m = (sum + n + 1)/d;

    dpowern = dR_ND->coKernels_ND_basis.length();

    H_last.SetLength(dpowern);

    for( e = m + N - 1; e != 0; e--)
    {
        Hnew.clear();

        if(verbose)
            cout<<"e = "<<e<<endl<<"\tadding new terms\n";

        if( e >= m)
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
                ZZ_p tmp;
                tmp = conv<ZZ_p>(Djm * rep(it->second) * rep(fact) );
                if(not IsZero(tmp))
                {
                    if(H.find(monomial) == H.end())
                        H[monomial].SetLength( dpowern, ZZ_p(0) );
                    H[monomial][0] += tmp;
                }
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

        int64_t Hlen;
        Hlen = H.size();
        int64_t * u_list;// Hlen*(n+1) matrix
        int64_t * v_list;// Hlen*(n+1) matrix
        fmpz_mat_struct * poly_list;// Hlen*(n+2) matrix
        fmpz * G_list; // Hlen*dpowern matrix
        fmpz_t modulus;
        fmpz_init(modulus);
        conv(modulus, ZZ_p::modulus());
#ifdef _OPENMP
        if(omp_get_max_threads() > 1)
        {
            u_list = (int64_t *)flint_malloc(sizeof(int64_t)*Hlen*(n+1));
            v_list = (int64_t *)flint_malloc(sizeof(int64_t)*Hlen*(n+1));
            poly_list = (fmpz_mat_struct *)flint_malloc(sizeof(fmpz_mat_struct)*Hlen*(n+2));
            G_list = _fmpz_vec_init(Hlen*dpowern);
            for(i = 0; i < Hlen*(n+2); i++)
                fmpz_mat_init(poly_list + i, dpowern, dpowern);
            j = 0;
            for( Hit = H.begin(); Hit != H.end(); Hit++)
            {
                if( verbose )
                {
                    cout<<"\treducing "<< Hit->first <<endl;
                }

                // choosing the main reduction direction

                if( e > 1)
                {
                    sum = 0;
                    for(i = 0; i<= n; i++)
                        v[i] = 0;
                    if( e > m)
                        shift = ei;
                    else
                        shift = v; //=0;

                    while( sum < d)
                    {
                        for( i = 0; i <= n; i++)
                        {
                           if( Hit->first[i] > v[i] + 1 + shift[i] )
                           {
                               v[i]++;
                               sum++;
                               break;
                           }
                        }
                    }
                }
                else
                    v = Hit->first;


                //end = (e > 1)? p: p - 1;
                dest =  p * Hit->first - end * v;

                conv(u_list + j*(n+1), Hit->first);
                conv(v_list + j*(n+1), v);
                conv(G_list + j*dpowern,Hit->second);
                dR_ND->get_ND_poly_flint(poly_list + j*(n+2), dest, v);
                j++;
            }

            #pragma omp parallel for shared(G_list, poly_list, end, dpowern, modulus) private(i)
            for( i = 0; i < Hlen; i++)
                dR_ND->reduce_vector_ND_poly_flint(G_list + i*dpowern, poly_list + i*(n+2), end, G_list + i*dpowern, modulus);

            for(i = 0; i < Hlen*(n+2); i++)
                fmpz_mat_clear(poly_list + i);

            for(i = 0; i < Hlen; i++)
            {
                conv(G_ZZ_p, G_list + i*dpowern, dpowern);
                conv(u, u_list + i*(n+1), n + 1);
                conv(v, v_list + i*(n+1), n + 1);

                if( ! IsZero(G_ZZ_p) )
                {
                    if( e > 1 )
                    {
                        final_dest = u - v;
                        assert( p*final_dest == p * u - end * v);
                        Hnewit = Hnew.find(final_dest);
                        if( Hnewit == Hnew.end() )
                            Hnew[final_dest] = G_ZZ_p;
                        else
                            Hnewit->second += G_ZZ_p;
                    }
                    else
                    {
                        F += (*(dR_ND->get_coKernels_ND_to_basis(u))) * G_ZZ_p;
                    }
                }
            }
            flint_free(u_list);
            flint_free(v_list);
            flint_free(poly_list);
            _fmpz_vec_clear(G_list,Hlen*dpowern);
            fmpz_clear(modulus);
            // reduced all terms from p*e to p*(e-1) and if e = 1 to basis
            H.swap(Hnew);
        }
        else
#endif
        {
            //Hlen = H.size()
            Hlen = 1;
            u_list = (int64_t *)flint_malloc(sizeof(int64_t)*Hlen*(n+1));
            v_list = (int64_t *)flint_malloc(sizeof(int64_t)*Hlen*(n+1));
            poly_list = (fmpz_mat_struct *)flint_malloc(sizeof(fmpz_mat_struct)*Hlen*(n+2));
            G_list = _fmpz_vec_init(Hlen*dpowern);
            for(i = 0; i < Hlen*(n+2); i++)
                fmpz_mat_init(poly_list + i, dpowern, dpowern);
            j = 0;
            for( Hit = H.begin(); Hit != H.end(); Hit++)
            {
                if( verbose )
                {
                    cout<<"e = "<<e<<endl<<"\tadding new terms\n";
                    cout<<"\treducing "<< Hit->first <<endl;
                }

                // choosing the main reduction direction
                if( e > 1)
                {
                    sum = 0;
                    for(i = 0; i<= n; i++)
                        v[i] = 0;

                    if( e > m)
                        shift = ei;
                    else
                        shift = v; //=0;

                    while( sum < d)
                    {
                        for( i = 0; i <= n; i++)
                        {
                           if( Hit->first[i]  > shift[i] + v[i] + 1 )
                           {
                               v[i]++;
                               sum++;
                               break;
                           }
                        }
                    }
                }
                else
                    v = Hit->first;

                //end = (e > 1)? p: p - 1;
                dest =  p * Hit->first - end * v;

                conv(u_list + j*(n+1), Hit->first);
                conv(v_list + j*(n+1), v);
                conv(G_list + j*dpowern, Hit->second);
                dR_ND->get_ND_poly_flint(poly_list + j*(n+2), dest, v);
                //j++;

                dR_ND->reduce_vector_ND_poly_flint(G_list, poly_list, end, G_list, modulus);

                conv(G_ZZ_p, G_list, dpowern);

                conv(u, u_list, n + 1);
                conv(v, v_list, n + 1);

                if( ! IsZero(G_ZZ_p) )
                {
                    if( e > 1 )
                    {
                        final_dest = u - v;
                        assert( p*final_dest == p * u - end * v);
                        Hnewit = Hnew.find(final_dest);
                        if( Hnewit == Hnew.end() )
                            Hnew[final_dest] = G_ZZ_p;
                        else
                            Hnewit->second += G_ZZ_p;
                    }
                    else
                    {
                        F += (*(dR_ND->get_coKernels_ND_to_basis(u))) * G_ZZ_p;
                    }
                }
            }
            for(i = 0; i < Hlen*(n+2); i++)
                fmpz_mat_clear(poly_list + i);
            flint_free(u_list);
            flint_free(v_list);
            flint_free(poly_list);
            _fmpz_vec_clear(G_list,Hlen*dpowern);
            fmpz_clear(modulus);
            // reduced all terms from p*e to p*(e-1) and if e = 1 to basis
            H.swap(Hnew);

        }
    }
    //Perform factorial division
    /*
    assert( fact == factorial<ZZ_p>( p*(m + N - 1) - 1) );
    val = valuation_of_factorial(p*(m + N - 1) - 1, p);
    DivRem(tmp, remainder,  factorial<ZZ>( p*(m + N - 1) - 1 ) , power_ZZ(p,val));
    assert(IsZero(remainder));
    fact = to_ZZ_p(tmp);
    */
    val = factorial_p_adic(fact,  p*(m + N - 1) - 1, p);
    inv(fact, fact );
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

Mat<ZZ_p> hypersurface_non_degenerate::frob_matrix_ND(Vec<int64_t> N)
{
    assert( n == (int64_t) N.length() );
    Mat<ZZ_p> F;
    F.SetDims( dR->coKernels_J_basis.length(), dR->coKernels_J_basis.length() );
    #if defined _OPENMP && defined NTL_THREADS
    if(omp_get_max_threads() > 1) {
      dR_ND->compute_everything_ND(false, false);
      compute_fpow(max(N) - 1);
    }
    ZZ_pContext context;
    context.save();
    #pragma omp parallel for schedule(dynamic)
    #endif
    for(int64_t i = 0; i < (int64_t) dR->coKernels_J_basis.length(); i++)
    {
        #if defined _OPENMP && defined NTL_THREADS
        context.restore();
        #endif
        int64_t j, m, sum;
        sum = 0;
        for( j = 0; j<= n; j++)
            sum += dR->coKernels_J_basis[i][j];

        m = (sum + n + 1)/d;

        timestamp_type wtime1, wtime2;
        double wall_time, user_time;
        user_time = get_cpu_time();
        get_timestamp(&wtime1);

        if(verbose)
            cout<<"Computing F("<<dR->coKernels_J_basis[i]<<") m = "<<m<<" N = "<<N[m-1]<<endl;
        if(N[m-1]>0)
            F[i] = frob_ND(i, N[m-1]);
        get_timestamp(&wtime2);
        wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
        user_time = get_cpu_time() - user_time;

        if(verbose)
        {
            cout <<"Computing F("<<dR->coKernels_J_basis[i]<<", N = "<<N[m-1]<<" = "<<F[i]<<endl;
            printf("Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );
        }
    }
    return transpose(F);
}

