// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
//
// dr_nd.cpp: routines for the class de_Rham_non_degenerate_local

#include "dr_nd.h"
#include "matrix.h"
#include "timing.h"
#ifdef _OPENMP
#include <omp.h>
#endif
# if __FLINT_RELEASE > 20800
#include <flint/fmpz_vec.h>
#include <flint/nmod.h>
#endif

void finitediff_flint_nmod(fmpz * result, fmpz_mat_struct * M_fmpz, const int64_t Mlength, const int64_t k, const fmpz * G, const fmpz_t &modulus);

//constructors
de_Rham_non_degenerate_local::de_Rham_non_degenerate_local(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<ZZ_p> f_vector, bool verbose, bool save_memory) {
    Vec< Vec<int64_t> > tuple_d;
    map< Vec<int64_t>, ZZ_p, vi64less> f_map;
    int64_t i;
    tuple_list_generator(tuple_d, d, (int64_t) (n + 1) );
    assert(tuple_d.length() == f_vector.length() );
    for( i = 0; i < (int64_t) f_vector.length() ; i++)
        f_map[ tuple_d[i] ] = f_vector[i];

    init_ND(p, precision, f_map, verbose,save_memory);

}
de_Rham_non_degenerate_local::de_Rham_non_degenerate_local(int64_t p, int64_t precision, map< Vec<int64_t>, ZZ_p, vi64less> f, bool verbose, bool save_memory) {
    init_ND(p, precision, f, verbose, save_memory);
};
de_Rham_non_degenerate_local::de_Rham_non_degenerate_local(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose, bool save_memory) {
    // naively lift fbar to f
    map< Vec<int64_t>, zz_p, vi64less>::iterator it;
    for(it = fbar.begin(); it != fbar.end(); ++it) {
        f[it->first] = rep(it->second);
    }
    init_ND(p, precision, f, verbose, save_memory);
}


de_Rham_non_degenerate_local::de_Rham_non_degenerate_local(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<zz_p> fbar_vector, bool verbose, bool save_memory)
{
    Vec< Vec<int64_t> > tuple_d;
    map< Vec<int64_t>, ZZ_p, vi64less> f_map;
    int64_t i;
    tuple_list_generator(tuple_d, d, (int64_t) (n + 1) );
    assert(tuple_d.length() == fbar_vector.length() );
    for( i = 0; i < (int64_t) fbar_vector.length() ; i++)
        f_map[ tuple_d[i] ] = rep(fbar_vector[i]);

    init_ND(p, precision, f_map, verbose,save_memory);

}

de_Rham_non_degenerate_local::de_Rham_non_degenerate_local(int64_t p, int64_t precision, int64_t n, int64_t d, bool verbose, bool save_memory)
{
    int64_t i, j;
    Vec<int64_t> B;
    zz_p x;
    map< Vec<int64_t>, zz_p, vi64less>::iterator it;

    bool boolean = true;
    this->precision = 1;
    this->verbose = verbose;
    this->n = n;
    this->d = d;

    tuple_list.SetLength( (n + 1) * d + 1);
    tuple_dict.SetLength( (n + 1) * d + 1);
    for(i = 0 ; i < (n + 1) * d + 1 ; i++)
    {
        tuple_list_generator( tuple_list[i], i, (int64_t) (n+1));
        for( j = 0; j < (int64_t)tuple_list[i].length(); j++)
        {
            tuple_dict[i][ ((tuple_list[i])[j]) ] = j;
        }
    }

    while(boolean)
    {
        solve_ND.clear();
        for( i = 0; i < (int64_t) tuple_list[d].length(); i++ )
        {
            x = random_zz_p();
            fbar[ tuple_list[d][i] ] = x;
            f[ tuple_list[d][i] ] = rep( x );

        }
        B = (get_solve_ND((d-1)*(n+1)+1))->first;
        boolean = bool( B.length() != 0 );
    }

    solve_ND.clear();
    this->precision = precision;


    init_ND(p, precision, f, verbose, save_memory);
}

de_Rham_non_degenerate_local::de_Rham_non_degenerate_local(const char* filename)
{
    int64_t i, j;
    ifstream file;
    Vec<zz_p> fbar_vector;
    pair< Vec<int64_t>, Mat<ZZ_p> > *solve_pair;
    file.open(filename);

    if(file.is_open())
    {
        file >> p;
        assert( p == (int64_t) zz_p::modulus() );
        file >> precision;
        assert( power_ZZ(p, precision) == ZZ_p::modulus() );
        file >> verbose;
        file >> save_memory;
        file >> n;
        file >> d;
        file >> Hilbert_J;
        file >> fbar_vector;

        tuple_list.SetLength( (n + 1) * d + 1);
        tuple_dict.SetLength( (n + 1) * d + 1);
        for(i = 0 ; i < (n + 1) * d + 1 ; i++)
        {
            tuple_list_generator( tuple_list[i], i, (int64_t) (n+1));
            for( j = 0; j < (int64_t)tuple_list[i].length(); j++)
            {
                tuple_dict[i][ ((tuple_list[i])[j]) ] = j;
            }
        }
        assert( tuple_list[d].length() == fbar_vector.length() );
        for(i = 0 ; i < (int64_t) fbar_vector.length() ; i++)
        {
            if( ! IsZero(fbar_vector[i]) )
            {
                f[ tuple_list[d][i] ] = rep( fbar_vector[i] );
                fbar[ tuple_list[d][i] ] = fbar_vector[i];
            }
        }

        for(i = 0; i < n; i++)
        {
            solve_pair = &( solve_J[(i+1)*d - (n+1)]);
            file >> solve_pair->first;
            file >> solve_pair->second;
            assert( Hilbert_J[(i+1)*d-(n+1)] == (int64_t) (solve_pair->first).length() );
            for( j = 0; j < (int64_t)(solve_pair->first).length(); j++)
            {
                append( coKernels_J_basis, tuple_list[ (i + 1) * d - (n + 1) ][ (solve_pair->first)[j] ] );
            }
        }
        solve_pair = &( solve_J[(d - 2) *  (n + 1) + 1]);
        file >> solve_pair->first;
        file >> solve_pair->second;

        for( i = 0; i < (int64_t)coKernels_J_basis.length(); i++ )
        {
            coKernels_J_basis_dict[ coKernels_J_basis[i] ] = i;
        }

        file >> Hilbert_ND;

        for( i = 0; i <= n; i++)
        {
            solve_pair = &( solve_ND[ i * d ]);
            file >> solve_pair->first;
            file >> solve_pair->second;
            assert( Hilbert_ND[i*d] == (int64_t) (solve_pair->first).length() );
            for( j = 0; j < (int64_t) (solve_pair->first).length(); j++ )
            {
                append( coKernels_ND_basis, tuple_list[ i * d ][ (solve_pair->first)[j] ] );
            }
        }
        for( i = 0; i < (int64_t) coKernels_ND_basis.length(); i++)
            coKernels_ND_basis_dict[ coKernels_ND_basis[i] ] = i;

        solve_pair = &( solve_ND[ (d - 1) * (n + 1) + 1 ] );
        file >> solve_pair->first;
        file >> solve_pair->second;
        assert( coKernels_ND_basis.length() == power_long(d, n) );
        file.close();
    }
    else
    {
        cout << "Could not open \'"<<filename<<"\'"<<endl;
        abort();
    }
}

void de_Rham_non_degenerate_local::init_ND(int64_t p, int64_t precision, map< Vec<int64_t>, ZZ_p, vi64less> f, bool verbose, bool save_memory)
{
    init(p, precision, f, verbose, save_memory);

    int64_t i, j;
    ZZX pol, H;
    pair< Vec<int64_t>, Mat<ZZ_p> >* solve_pair;
    //Hilbert_ND = (1 + x + x^2 + ... x^(d-1) )^(n+1)
    Hilbert_ND.SetLength((d-1)*(n+1)+1);

    for(i = 0; i < d ; i++)
        SetCoeff(pol,i);

    H = 1;

    for(i = 0; i <= n; i++)
        mul(H, H, pol);

    for(i = 0; i <= (d-1)*(n+1); i++)
        Hilbert_ND[i] = (int64_t) to_ulong(coeff(H,i));

    if(verbose)
        cout <<"Hilbert_ND = "<<Hilbert_ND<<endl;

    if(verbose)
        cout << "Asserting that f is non degenerate"<<endl;

    solve_pair = get_solve_ND( (d - 1)*(n+1) + 1 );
    assert( (solve_pair->first).length() == 0);

    for( i = 0; i <= n; i++)
    {
        solve_pair = get_solve_ND( i * d );
        assert( Hilbert_ND[i*d] == (int64_t) (solve_pair->first).length() );
        for( j = 0 ; j < (int64_t) (solve_pair->first).length(); j++)
            append(coKernels_ND_basis, tuple_list[ i*d ][ (solve_pair->first)[j] ] );
    }

    assert( coKernels_ND_basis.length() == power_long(d, n) );

    for( i = 0; i < (int64_t) coKernels_ND_basis.length(); i++)
        coKernels_ND_basis_dict[ coKernels_ND_basis[i] ] = i;

    if(verbose)
    {
        cout << "coKernels_ND_basis = ";
        cout << coKernels_ND_basis;
        cout <<endl;
    }
}

bool de_Rham_non_degenerate_local::save( const char* filename)
{
    if(verbose)
        cout <<"Saving de_Rham_non_degenerate_local to: "<< filename << endl;

    ofstream fileswp;
    char fileswpname[sizeof(filename)+4];
    strcpy(fileswpname, filename);
    strcat(fileswpname,".swp");
    fileswp.open(fileswpname);

    int64_t i;
    ofstream file;
    Vec<zz_p> fbar_vector;
    map< Vec<int64_t>, zz_p, vi64less>::const_iterator it;
    pair< Vec<int64_t> , Mat<ZZ_p> > *solve_pair;

    fbar_vector.SetLength( tuple_list[d].length() );
    for( it = fbar.begin(); it != fbar.end(); it++)
        fbar_vector[ tuple_dict[d][it->first] ] = it->second;

    file.open(filename, ios::out | ios::trunc);



    if(file.is_open())
    {
        file << p <<endl;
        file << precision <<endl;
        file << verbose <<endl;
        file << save_memory << endl;
        file << n <<endl;
        file << d <<endl;
        file << Hilbert_J <<endl;
        file << fbar_vector << endl;

        for( i = 0; i < n ; i++)
        {
            solve_pair = get_solve_J( (i + 1) * d - (n + 1) );
            file << solve_pair->first <<endl;
            file << solve_pair->second <<endl;
        }
        solve_pair = get_solve_J((d - 2)*(n + 1) + 1);
        file << solve_pair->first <<endl;
        file << solve_pair->second <<endl;

        file << Hilbert_ND << endl;

        for( i = 0; i <= n ; i++ )
        {
            solve_pair = get_solve_ND( i * d );
            file << solve_pair->first <<endl;
            file << solve_pair->second <<endl;
        }

        solve_pair = get_solve_ND( (d - 1) * (n + 1) + 1 );
        file << solve_pair->first <<endl;
        file << solve_pair->second <<endl;



        file.close();
        fileswp.close();
        remove(fileswpname);
        return true;
    }
    else
    {
        cout << "Wasn't able to open "<<filename<<" or "<<fileswpname<<endl;
        return false;
    }
}


void de_Rham_non_degenerate_local::compute_everything_ND(bool J, bool ND_ZZ)
{
    if(J)
      compute_everything_J();
    int64_t i,j;
    Vec<int64_t> u;
    for(i = 0; i < (int64_t) tuple_list[d].length(); i++)
    {
        u = tuple_list[d][i];
        compute_reduction_matrix_ND(u);
        // Not used anywhere by default
        if(ND_ZZ)
          compute_reduction_matrix_ND_ZZ(u);
        compute_reduction_matrix_ND_poly( u );
        compute_reduction_matrix_ND_poly_ZZ(u);



        bool positive = true;
        for(j = 0; j <= n; j++)
            positive = positive && (u[j] > 0);
        if( positive  && (inclusion_matrix_ND_dict.find(u) == inclusion_matrix_ND_dict.end()) )
        {
            if(verbose)
                cout << "Computing the inclusion matrix ND for u = "<<u<<endl;
            compute_inclusion_matrix_ND(u);
        }
        if( positive && ( coKernels_ND_to_basis_dict.find(u) == coKernels_ND_to_basis_dict.end() ) )
        {
            if(verbose)
                cout << "Computing coKernels ND to basis for u = "<<u<<endl;
            compute_coKernels_ND_to_basis(u);
        }


    }

}

pair< Vec<int64_t>, Mat<ZZ_p> >* de_Rham_non_degenerate_local::get_solve_ND(const int64_t level)
{
    map< int64_t, pair< Vec<int64_t>, Mat<ZZ_p> > >::iterator it;
    it = solve_ND.find(level);
    if( it != solve_ND.end() )
        return &(it->second);

    timestamp_type wtime1, wtime2;
    double wall_time, user_time;
    Mat<ZZ_p> MND;
    matrix_ND(MND, level);
    if(verbose) {
        user_time = get_cpu_time();
        get_timestamp(&wtime1);
        cout<<"Computing and solving the ND matrix of relations at degree = "<<level<<" ( "<<MND.NumRows()<<"x"<<MND.NumCols()<<" )."<<endl;
    }
    solve_system_padic( solve_ND[level].first, solve_ND[level].second, MND, precision);
    if (verbose) {
        get_timestamp(&wtime2);
        wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
        user_time = get_cpu_time() - user_time;
        printf("Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );
    }
    return &(solve_ND[level]);
}

Mat<ZZ_p> de_Rham_non_degenerate_local::get_reduction_matrix_ND(const Vec<int64_t> u, const Vec<int64_t> v)
{

    int64_t i, sum_v, sum_u;
    ZZ monomial_evaluated;
    Mat<ZZ_p> result;
    sum_u = 0;
    sum_v = 0;
    for(i = 0; i <= n ; i++)
    {
        sum_v += v[i];
        sum_u += u[i];
    }
    assert(sum_v == d);
    assert(sum_u % d == 0);

    map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less>::const_iterator it;

    it = compute_reduction_matrix_ND(v);

    result.SetDims( coKernels_ND_basis.length() , coKernels_ND_basis.length());

    map< Vec<int64_t>, Mat<ZZ_p>, vi64less>::const_iterator itM;

    for(itM = (it->second).begin(); itM != (it->second).end(); itM++)
    {
        NTL::set(monomial_evaluated);
        for( i = 0; i <= n ; i++)
            mul(monomial_evaluated, monomial_evaluated, power_ZZ(u[i],(itM->first)[i]) );

        result += conv<ZZ_p>(monomial_evaluated) * itM->second;
    }

    return result;
}

Mat<ZZ> de_Rham_non_degenerate_local::get_reduction_matrix_ND_ZZ(const Vec<int64_t> u, const Vec<int64_t> v)
{

    int64_t i, sum_v, sum_u;
    ZZ monomial_evaluated;
    Mat<ZZ> result;
    sum_u = 0;
    sum_v = 0;
    for(i = 0; i <= n ; i++)
    {
        sum_v += v[i];
        sum_u += u[i];
    }
    assert(sum_v == d);
    assert(sum_u % d == 0);

    map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ>, vi64less> , vi64less>::const_iterator it;

    it = compute_reduction_matrix_ND_ZZ(v);

    result.SetDims( coKernels_ND_basis.length() , coKernels_ND_basis.length());

    map< Vec<int64_t>, Mat<ZZ>, vi64less>::const_iterator itM;

    for(itM = (it->second).begin(); itM != (it->second).end(); itM++)
    {
        NTL::set(monomial_evaluated);
        for( i = 0; i <= n ; i++)
            mul(monomial_evaluated, monomial_evaluated, power_ZZ(u[i],(itM->first)[i]) );
        rem(monomial_evaluated, monomial_evaluated, ZZ_p::modulus());
        result += monomial_evaluated * itM->second;
    }
    return result;
}

void de_Rham_non_degenerate_local::reduce_vector_ND(Vec<ZZ_p> &result, const Vec<int64_t> u, const Vec<int64_t> v, const Vec<ZZ_p> G)
{
    int64_t i, sum_v, sum_u;
    ZZ monomial_evaluated;
    sum_u = 0;
    sum_v = 0;
    for(i = 0; i <= n ; i++)
    {
        sum_v += v[i];
        sum_u += u[i];
    }
    assert(sum_v == d);
    assert(sum_u % d == 0);


    map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less>::const_iterator it;

    it = compute_reduction_matrix_ND(v);

    result.SetLength( coKernels_ND_basis.length());
    clear(result);

    map< Vec<int64_t>, Mat<ZZ_p>, vi64less>::const_iterator itM;

    for(itM = (it->second).begin(); itM != (it->second).end(); itM++)
    {
        NTL::set(monomial_evaluated);
        for( i = 0; i <= n ; i++)
            mul(monomial_evaluated, monomial_evaluated, power_ZZ(u[i],(itM->first)[i]) );

        result +=  itM->second * (conv<ZZ_p>(monomial_evaluated) * G);
    }

}

void de_Rham_non_degenerate_local::reduce_vector_ND_ZZ(Vec<ZZ> &result, const Vec<int64_t> u, const Vec<int64_t> v, const Vec<ZZ> G)
{
    int64_t i, sum_v, sum_u;
    ZZ monomial_evaluated;
    sum_u = 0;
    sum_v = 0;
    for(i = 0; i <= n ; i++)
    {
        sum_v += v[i];
        sum_u += u[i];
    }
    assert(sum_v == d);
    assert(sum_u % d == 0);

    map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ>, vi64less> , vi64less>::const_iterator it;

    it = compute_reduction_matrix_ND_ZZ(v);

    map< Vec<int64_t>, Mat<ZZ>, vi64less>::const_iterator itM;
    itM = (it->second).begin();
    NTL::set(monomial_evaluated);
    for( i = 0; i <= n ; i++)
        mul(monomial_evaluated, monomial_evaluated, power_ZZ(u[i],(itM->first)[i]) );
    result = itM->second * ( monomial_evaluated * G);
    itM++;
    for(; itM != (it->second).end(); itM++)
    {
        NTL::set(monomial_evaluated);
        for( i = 0; i <= n ; i++)
            mul(monomial_evaluated, monomial_evaluated, power_ZZ(u[i],(itM->first)[i]) );
        result += itM->second * ( monomial_evaluated * G );
    }
}

void de_Rham_non_degenerate_local::reduce_vector_ND_poly(Vec<ZZ_p> &result, const Vec<int64_t> u, const Vec<int64_t> v,const int64_t iterations, const Vec<ZZ_p> G)
{
    int64_t i, sum_v, sum_u, dpowern;
    ZZ_p xpower;
    int64_t x;
    ZZ monomial_evaluated;
    sum_u = 0;
    sum_v = 0;
    for(i = 0; i <= n ; i++)
    {
        sum_v += v[i];
        sum_u += u[i];
    }
    assert(sum_v == d);
    assert(sum_u % d == 0);

    dpowern = coKernels_ND_basis.length();
    map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less>::const_iterator it;
    map< Vec<int64_t>, Mat<ZZ_p>, vi64less>::const_iterator itM;
    Vec<ZZ_p> Gout, Gin, Gtmp;

    Vec< Mat<ZZ_p> > poly;
    poly.SetLength(n+2);
    for( i = 0; i < n+2; i++)
        poly[i].SetDims(dpowern,dpowern);


    it = compute_reduction_matrix_ND_poly(v);


    for( itM = (it->second).begin(); itM != (it->second).end() ; itM++)
    {
        NTL::set(monomial_evaluated);
        sum_u = (itM->first)[n+1];
        for( i = 0; i <= n ; i++)
        {
            sum_u += (itM->first)[i];
            mul(monomial_evaluated, monomial_evaluated, power_ZZ(u[i], (itM->first)[i]) );
        }
        if(sum_u == n+1)
            (poly[ itM->first[n+1] ])[0] += conv<ZZ_p>(monomial_evaluated) * (itM->second)[0];
        else
            poly[ (itM->first)[n+1] ] += itM->second * conv<ZZ_p>(monomial_evaluated);

    }

    //assert( poly[0] ==  get_reduction_matrix_ND(u,v) );
    Gin = G;
    //Mat<ZZ_p> M;
    for(x = iterations -1 ; x != (int64_t) -1 ; x--)
    {
        xpower = 1;
        //M = poly[0];
        Gout = poly[0] * Gin;
        for(i = 1; i <= n; i++)
        {
            xpower *= x; // xpower = x^i
            mul(Gtmp, poly[i], Gin);
            mul(Gtmp, xpower, Gtmp);
            add(Gout, Gout, Gtmp);
            //Gout += xpower * ( poly[i] * Gin );
            //M += xpower * poly[i];
        }
        xpower *= x;
        Gout[0] += xpower * (poly[n+1][0] * Gin);
        swap(Gin,Gout);
        //M += xpower * poly[n+1];
        //assert( M == get_reduction_matrix_ND(u + x*v ,v) );
    }

    result = Gin;
}




void de_Rham_non_degenerate_local::reduce_vector_ND_poly_ZZ(Vec<ZZ> &result, const Vec<int64_t> u, const Vec<int64_t> v, const int64_t iterations, const Vec<ZZ> G)
{
    int64_t i, j, k, sum_v, sum_u, dpowern;
    int64_t x;
    ZZ monomial_evaluated, tmp;
    sum_u = 0;
    sum_v = 0;
    assert( (int64_t) v.length() == n + 1);
    assert( (int64_t) u.length() == n + 1);
    for(i = 0; i <= n ; i++)
    {
        sum_v += v[i];
        sum_u += u[i];
    }
    assert(sum_v == d);
    assert(sum_u % d == 0);

    dpowern = coKernels_ND_basis.length();
    map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ>, vi64less> , vi64less>::const_iterator it;
    map< Vec<int64_t>, Mat<ZZ>, vi64less>::const_iterator itM;
    Vec<ZZ> Gout,Gin,Gtmp;
    Vec< Mat<ZZ> > poly, poly_unreduced;
    poly.SetLength(n+2);
    for( i = 0; i < n+2; i++)
        poly[i].SetDims(dpowern,dpowern);
    poly_unreduced = poly;

    it = compute_reduction_matrix_ND_poly_ZZ(v);


    for( itM = (it->second).begin(); itM != (it->second).end() ; itM++)
    {
        NTL::set(monomial_evaluated);
        sum_u = (itM->first)[n+1];
        for( i = 0; i <= n ; i++)
        {
            sum_u += (itM->first)[i];
            mul(monomial_evaluated, monomial_evaluated, power_ZZ(u[i], (itM->first)[i]) );
        }
        if(sum_u == n+1)
            (poly_unreduced[ itM->first[n+1] ])[0] += monomial_evaluated * (itM->second)[0];
        else
            poly_unreduced[ (itM->first)[n+1] ] += itM->second * monomial_evaluated;
    }

    for(i = 0; i < n+2; i++)
    {
        for(j = 0; j < dpowern; j++)
        {
            for(k = 0; k < dpowern; k++)
            {
                rem(poly[i][j][k], poly_unreduced[i][j][k], ZZ_p::modulus() );
            }
        }
    }

    //for(i = 1; i < dpowern; i++)
    //    assert(IsZero(poly[n+1][i]));
    //assert( poly[0] ==  get_reduction_matrix_ND(u,v) );
    Gin = G;
    ZZ xpower;

    //Mat<ZZ_p> M;
    for(x = iterations -1 ; x != (int64_t) -1 ; x--)
    {
        xpower = 1;
        mul(Gout, poly[0], Gin);
        //M = poly[0];
        for(i = 1; i <= n; i++)
        {
            xpower *= x; // xpower = x^i
            mul(Gtmp, poly[i], Gin);
            Gout += xpower * Gtmp;
            //M += xpower * poly[i];
        }
        xpower *= x;
        InnerProduct(tmp, poly[n+1][0], Gin);
        Gout[0] += xpower * tmp;
        //M += xpower * poly[n+1];
        //assert( M == get_reduction_matrix_ND(u + x*v ,v) );

        for(j = 0; j < dpowern; j++)
            rem(Gin[j], Gout[j], ZZ_p::modulus() );
    }

    result = Gin;
}

void de_Rham_non_degenerate_local::get_ND_poly_flint(fmpz_mat_struct * result, const Vec<int64_t> u, const Vec<int64_t> v)
{
    int64_t i, j, k, sum_v, sum_u, dpowern;
    ZZ monomial_evaluated, tmp;
    sum_u = 0;
    sum_v = 0;
    assert( (int64_t) v.length() == n + 1);
    assert( (int64_t) u.length() == n + 1);

    for(i = 0; i <= n ; i++)
    {
        sum_v += v[i];
        sum_u += u[i];
    }
    assert(sum_v == d);
    assert(sum_u % d == 0);

    dpowern = coKernels_ND_basis.length();
    map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ>, vi64less> , vi64less>::const_iterator it;
    map< Vec<int64_t>, Mat<ZZ>, vi64less>::const_iterator itM;
    Vec< Mat<ZZ> > poly, poly_unreduced;
    poly.SetLength(n+2);
    for( i = 0; i < n+2; i++)
    {
        poly[i].SetDims(dpowern,dpowern);
    }
    poly_unreduced = poly;

    it = compute_reduction_matrix_ND_poly_ZZ(v);


    for( itM = (it->second).begin(); itM != (it->second).end() ; itM++)
    {
        NTL::set(monomial_evaluated);
        sum_u = (itM->first)[n+1];
        for( i = 0; i <= n ; i++)
        {
            sum_u += (itM->first)[i];
            mul(monomial_evaluated, monomial_evaluated, power_ZZ(u[i], (itM->first)[i]) );
        }
        if(sum_u == n+1)
            (poly_unreduced[ itM->first[n+1] ])[0] += monomial_evaluated * (itM->second)[0];
        else
            poly_unreduced[ (itM->first)[n+1] ] += itM->second * monomial_evaluated;
    }
    for(i = 0; i < n+2; i++)
    {
        for(j = 0; j < dpowern; j++)
        {
            for(k = 0; k < dpowern; k++)
            {
                rem(poly[i][j][k], poly_unreduced[i][j][k], ZZ_p::modulus() );
            }
        }
        conv(result + i,poly[i]);
    }
}




void de_Rham_non_degenerate_local::reduce_vector_ND_poly_flint(fmpz * result, fmpz_mat_struct * poly, const int64_t iterations, const fmpz * G, const fmpz_t& modulus)
{
    int64_t dpowern = fmpz_mat_ncols(poly + 0);

    if(fmpz_size(modulus) > 1)
    {
        fmpz_t xpower;
        fmpz_t tmp;
        fmpz_init(xpower);
        fmpz_init(tmp);
        fmpz * Gin;
        fmpz * Gout;
        fmpz * Gtmp;

        Gin = _fmpz_vec_init(dpowern);
        _fmpz_vec_set(Gin, G, dpowern);
        Gout = _fmpz_vec_init(dpowern);
        Gtmp =  _fmpz_vec_init(dpowern);

        for(int64_t x = iterations -1 ; x >= 0 ; x--)
        {
            fmpz_one(xpower);
            mul(Gout, poly + 0, Gin);
            for(int64_t i = 1; i <= n; i++)
            {
                fmpz_mul_ui(xpower, xpower, x); // xpower = x^i
                mul(Gtmp, poly + i, Gin);
                _fmpz_vec_scalar_mul_fmpz(Gtmp, Gtmp, dpowern, xpower);
                _fmpz_vec_add(Gout, Gout, Gtmp, dpowern);
            }
            fmpz_mul_ui(xpower, xpower, x);
            mul(tmp,(poly + n+1)->rows[0], Gin,dpowern);
            fmpz_addmul(Gout + 0, xpower, tmp);
            _fmpz_vec_scalar_mod_fmpz(Gin, Gout, dpowern, modulus);
        }
        fmpz_clear(xpower);
        fmpz_clear(tmp);
        _fmpz_vec_set(result, Gin, dpowern);
        _fmpz_vec_clear(Gin, dpowern);
        _fmpz_vec_clear(Gout, dpowern);
        _fmpz_vec_clear(Gtmp, dpowern);
    }
    else
    {
        if(iterations >  n + 1)
            finitediff_flint_nmod(result, poly, n + 2, iterations, G, modulus);
        else
        {
            nmod_t mod;
            int nlimbs;
            nmod_init(&mod, fmpz_get_ui(modulus));
            int64_t j,k;

            mp_limb_t xpower;
            mp_limb_t tmp;

            mp_ptr Gout;
            mp_ptr Gin;
            mp_ptr Gtmp;


            Gin = _nmod_vec_init(dpowern);
            for(int64_t i = 0; i < dpowern; i++)
            {
                Gin[i] = fmpz_get_ui(G + i);
            }
            Gout = _nmod_vec_init(dpowern);
            Gtmp =  _nmod_vec_init(dpowern);

            nlimbs = _nmod_vec_dot_bound_limbs(dpowern, mod);

            nmod_mat_struct * poly_nmod;
            poly_nmod = (nmod_mat_struct *)flint_malloc(sizeof(nmod_mat_struct) * (n+2));

            for(int64_t i = 0; i < n + 2; i++)
            {
                nmod_mat_init(poly_nmod + i, dpowern, dpowern, mod.n);
                for(j = 0; j < dpowern; j++)
                    for(k = 0; k < dpowern; k++)
                        nmod_mat_entry(poly_nmod + i, j, k) = fmpz_get_ui( fmpz_mat_entry(poly + i, j, k));
            }
            for(int64_t x = iterations -1 ; x >= 0; x--)
            {
                xpower = 1;
                mul(Gout, poly_nmod + 0, Gin, nlimbs);
                for(int64_t i = 1; i <= n; i++)
                {
                    xpower = nmod_mul(xpower, x, mod); // xpower = x^i
                    mul(Gtmp, poly_nmod + i, Gin, nlimbs);
                    _nmod_vec_scalar_addmul_nmod(Gout, Gtmp, dpowern, xpower, mod);
                }
                xpower = nmod_mul(xpower, x, mod);
                tmp =_nmod_vec_dot( (poly_nmod + n+1)->rows[0], Gin, dpowern, mod, nlimbs);
                //NMOD_VEC_DOT(tmp, i, dpowern, (poly_nmod + n+1)->rows[0], Gin, mod, nlimbs);
                tmp = nmod_mul(tmp, xpower, mod);
                Gout[0] = nmod_add(Gout[0], tmp, mod);
                swap(Gin,Gout);
            }

            for(int64_t i = 0; i < n + 2; i++)
                nmod_mat_clear(poly_nmod + i);
            flint_free(poly_nmod);
            for(int64_t i = 0; i < dpowern; i++)
                fmpz_set_ui(result + i, Gin[i]);

            _nmod_vec_clear(Gout);
            _nmod_vec_clear(Gin);
            _nmod_vec_clear(Gtmp);
        }
    }
}



map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less>::const_iterator de_Rham_non_degenerate_local::compute_reduction_matrix_ND(const Vec<int64_t> v)
{
    map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less>::const_iterator it;
    it = reduction_matrix_ND_dict.find(v);
    if( it != reduction_matrix_ND_dict.end() )
        return it;
    else
    {
        if(verbose)
            cout<<"Computing the reduction matrix ND for v = "<<v<<endl;
        map< Vec<int64_t>, Mat<ZZ_p>, vi64less>* result;

        result = &(reduction_matrix_ND_dict[v]);

        int64_t dsum = 0;
        for(int64_t i = 0; i <= n; i++)
            dsum += v[i];

        assert(dsum == d);

        Vec<int64_t> zero;
        zero.SetLength(n+1);
        for(int64_t i = 0 ; i <= n; i++)
            zero[i] = 0;

        Vec<Vec<int64_t> > U;
        U.SetLength(n+1);
        for(int64_t i = 0; i <= n; i++)
        {
            zero[i]++;
            U[i] = zero;
            zero[i]--;
        }

        int64_t dpowern = coKernels_ND_basis.length(); //power_long(d, n);


        Vec< pair< Vec<int64_t>, Mat<ZZ_p> >* > solve;
        solve.SetLength(n + 1);
        for(int64_t i = 0; i < n; i++)
            solve[i] = get_solve_ND( (i+1)*d );
        solve[n] = get_solve_ND( (d-1) * (n + 1) + 1);


        #if defined _OPENMP && defined NTL_THREADS
        ZZ_pContext context;
        context.save();
        #pragma omp parallel for
        #endif
        for(int64_t coordinate_of_monomial = 0; coordinate_of_monomial < dpowern ; ++coordinate_of_monomial) {
            #if defined _OPENMP && defined NTL_THREADS
            context.restore();
            #endif
            map< Vec<int64_t>, Vec<ZZ_p>, vi64less > H, Hnew;
            Vec<ZZ_p>* Hnew_zero;
            Vec<ZZ_p>* Hnew_Ui;
            int64_t reduction_steps;

            int64_t sum = 0;
            for(int64_t i = 0 ; i <=n ; i++)
                sum += coKernels_ND_basis[ coordinate_of_monomial ][i];

            assert( sum % d == 0);
            int64_t l = sum/d;

            if(l == n)
            {
                // d- n=  (d- 1) * (n + 1) + 1  - n*d
                Vec<int64_t> r = tweak(v, d - n);
                Vec<int64_t> s = v - r; // sum(s) = d - n

                // for now  H is just a monomial
                int64_t coordinate_of_H = tuple_dict[ (d-1)*(n+1)+ 1][ coKernels_ND_basis[ coordinate_of_monomial] + s ];

                int64_t dim_Fi = tuple_list[ d * n -n ].length(); // = (d-1)*(n+1) + 1 - d
                int64_t dim_Hnew = tuple_list[n*d].length();

                // F = solve[n] * H
                Vec<ZZ_p> F;
                F.SetLength( (solve[n]->second).NumRows() );
                for(int64_t i = 0 ; i < (int64_t) (solve[n]->second).NumRows(); i++)
                    F[i] = (solve[n]->second)[i][ coordinate_of_H ];


                Hnew_zero =  &(Hnew[zero]);
                Hnew_zero->SetLength(dim_Hnew);

                for(int64_t i = 0; i <= n; i++)
                {
                    Hnew_Ui = &(Hnew[U[i]]);
                    Hnew_Ui->SetLength(dim_Hnew);
                    for(int64_t coordinate_of_monomial_of_Fi = 0; coordinate_of_monomial_of_Fi < dim_Fi; ++coordinate_of_monomial_of_Fi )
                    {
                        // x^(u+r) H / x0 ... xn = x^(u+r) * \sum Fi xi di f / x0 ... xn
                        // ps: there are no coKernels at this level
                        // ~ \sum di( x^(u+r) * xi * Fi / x0 ... xn) = x^u * newH / x0 ... xn
                        // Let x^w be a monomial of Fi
                        // x^(u + r) * x^w * xi * di f / x0 ... xn ~ (u[i] - 1 + r[i] + w[i] + 1) x^u * x^r * x^w
                        int64_t coordinate_of_Hnew = tuple_dict[n * d][ r + tuple_list[d * n - n][coordinate_of_monomial_of_Fi] ];
                        (*Hnew_Ui)[ coordinate_of_Hnew ] += F[ i * dim_Fi + coordinate_of_monomial_of_Fi ];
                        (*Hnew_zero)[ coordinate_of_Hnew ] += (r[i] + tuple_list[d * n - n][coordinate_of_monomial_of_Fi][i]) *  F[ i * dim_Fi + coordinate_of_monomial_of_Fi ];
                    }
                }
                H.swap(Hnew);
                //we already did one reduction
                reduction_steps = l;
            }
            else
            {
                reduction_steps = l + 1;
                H.clear();
                Hnew_zero = &(H[zero]);
                Hnew_zero->SetLength( tuple_list[(l + 1) * d].length() );
                NTL::set((*Hnew_zero)[ tuple_dict[(l + 1) * d][ coKernels_ND_basis[ coordinate_of_monomial] + v ]]);
                // x^v * monomial_of_Gl, deg = l * d + d
            }

            for(int64_t k = reduction_steps-1; k >= 0; k-- )
            {
                /*
                 * H is polynomial with vector coefficients
                 * each vector represents a polynomial of degree (k+1)*d that we need to reduce to degree k*d
                 * Hnew will be the result of this reduction
                 *
                 * we will process each monomial of H individually
                 */
                int64_t dim_Hnew = tuple_list[ k * d ].length();
                int64_t dim_Fi = dim_Hnew;
                Hnew.clear();
                for(map< Vec<int64_t>, Vec<ZZ_p>, vi64less >::const_iterator itH = H.begin();  itH != H.end(); itH++)
                {
                    // How to write monomial = \sum Fi xi di f + coKernls
                    // deg Fi = k * d
                    Vec<ZZ_p> F;
                    F = solve[k]->second * itH->second;

                    Hnew_zero = &(Hnew[itH->first]);
                    Hnew_zero->SetLength(dim_Hnew);
                    for(int64_t  i = 0; i <= n; i++)
                    {
                        Hnew_Ui = &(Hnew[itH->first + U[i] ]);
                        Hnew_Ui->SetLength(dim_Hnew);

                        for(int64_t coordinate_of_monomial_of_Fi = 0; coordinate_of_monomial_of_Fi < dim_Fi; coordinate_of_monomial_of_Fi++)
                        {
                            /*
                             * x^u H / x0...xn = x^u * (\sum Fi xi di f )/ x0... xn
                             * ~ \sum di( x^u * xi * Fi / x0 ... xn) + x^u * coKernels / x0 ... xn
                             * = x^u * (newH + coKernels) / x0 ... xn
                             * let x^w be a monomial of Fi
                             * x^u * x^w *  xi * di f / x0 ... xn ~ (u[i] - 1 + w[i] + 1) x^u * x^w/ x0 ... xn
                             */
                            (*Hnew_Ui)[ coordinate_of_monomial_of_Fi] += F[i * dim_Fi + coordinate_of_monomial_of_Fi];
                            (*Hnew_zero)[ coordinate_of_monomial_of_Fi] += tuple_list[k * d][coordinate_of_monomial_of_Fi][i] * F[i * dim_Fi + coordinate_of_monomial_of_Fi];
                        }
                    }

                    if( (int64_t) (k+1) <= n )
                    {
                        #pragma omp critical
                        for(int64_t coordinate_of_coKernel_element = 0; coordinate_of_coKernel_element < (int64_t) solve[k]->first.length(); coordinate_of_coKernel_element++ )
                        {
                            map< Vec<int64_t>, Mat<ZZ_p>, vi64less>::iterator it_result = result->find(itH->first);
                            if( it_result == result->end() )
                            {
                                (*result)[itH->first].SetDims(dpowern, dpowern);
                                it_result = result->find(itH->first);
                                assert(it_result != result->end());
                            }
                            (it_result->second)[ coKernels_ND_basis_dict[ tuple_list[(k+1)*d][ (solve[k]->first)[coordinate_of_coKernel_element] ] ] ][ coordinate_of_monomial ] += F[(n+1)*dim_Fi + coordinate_of_coKernel_element];

                        }
                    }
                }

                H.swap(Hnew);

                if( k == 0)
                {
                    #pragma omp critical
                    // no more reductions needed, only the constant term is left on H
                    for(map< Vec<int64_t>, Vec<ZZ_p>, vi64less >::const_iterator itH = H.begin(); itH != H.end(); itH++ )
                    {
                        map< Vec<int64_t>, Mat<ZZ_p>, vi64less>::iterator it_result = result->find(itH->first);
                        if( it_result == result->end() )
                        {
                            ((*result)[itH->first]).SetDims(dpowern, dpowern);
                            it_result = result->find(itH->first);
                            assert(it_result != result->end());
                        }
                        (it_result->second)[0][coordinate_of_monomial] += (itH->second)[0];
                    }
                }
            }//loop over k ends
        }
        #pragma omp critical
        for(map< Vec<int64_t>, Mat<ZZ_p>, vi64less>::iterator it_result = result->begin(); it_result != result->end(); )
        {
            if( IsZero(it_result->second) )
            {
                result->erase(it_result++);
            }
            else
            {
                it_result++;
            }
        }
        return reduction_matrix_ND_dict.find(v);
    }
}

map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ>, vi64less> , vi64less>::const_iterator de_Rham_non_degenerate_local::compute_reduction_matrix_ND_ZZ(const Vec<int64_t> v)
{
    map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ>, vi64less> , vi64less>::const_iterator it_ZZ;
    it_ZZ = reduction_matrix_ND_ZZ_dict.find(v);
    if( it_ZZ != reduction_matrix_ND_ZZ_dict.end() )
        return it_ZZ;
    else
    {
        if(verbose)
            cout<<"Computing the reduction matrix ND (ZZ) for v = "<<v<<endl;
        map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less>::const_iterator it;
        map< Vec<int64_t>, Mat<ZZ_p>, vi64less>::const_iterator it_result;
        bool existQ;
        map< Vec<int64_t>, Mat<ZZ>, vi64less> &result_zz = reduction_matrix_ND_ZZ_dict[v];
        it = reduction_matrix_ND_dict.find(v);
        existQ = (it != reduction_matrix_ND_dict.end());
        if( !existQ)
            it = compute_reduction_matrix_ND(v);

        for(it_result = (it->second).begin(); it_result != (it->second).end(); it_result++)
            result_zz[it_result->first] = conv< Mat<ZZ> >(it_result->second);
        if( !existQ && save_memory )
            reduction_matrix_ND_dict.erase(v);

        return reduction_matrix_ND_ZZ_dict.find(v);
    }
}





map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less>::const_iterator de_Rham_non_degenerate_local::compute_reduction_matrix_ND_poly(const Vec<int64_t> v)
{
    map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less>::const_iterator it;
    it = reduction_matrix_ND_poly_dict.find(v);
    if( it != reduction_matrix_ND_poly_dict.end() )
        return it;
    else
    {
        if(verbose)
            cout<<"Computing the reduction matrix poly ND for v = "<<v<<endl;

        bool existQ;
        it = reduction_matrix_ND_dict.find(v);
        existQ = (it != reduction_matrix_ND_dict.end());

        if(!existQ)
            it = compute_reduction_matrix_ND(v);

        map< Vec<int64_t>, Mat<ZZ_p>, vi64less>::const_iterator it_matrix;
        map< Vec<int64_t>, Mat<ZZ_p>, vi64less>::iterator it_result;
        map< Vec<int64_t>, Mat<ZZ_p> , vi64less>* result;
        map< Vec<int64_t>, int64_t, vi64less> monomial_expanded;
        map< Vec<int64_t>, int64_t, vi64less>::const_iterator it_monomial;
        int64_t i, sum, dpowern;
        dpowern = coKernels_ND_basis.length();

        result = &(reduction_matrix_ND_poly_dict[v]);

        for( it_matrix =  (it->second).begin(); it_matrix != (it->second).end(); it_matrix++)
        {
            sum = 0;
            for(i = 0; i <= n; i++)
                sum += it_matrix->first[i];

            monomial_expanded = change_of_variables_monomial( it_matrix->first, v);

            for( it_monomial = monomial_expanded.begin(); it_monomial != monomial_expanded.end(); it_monomial++)
            {
               it_result = result->find(it_monomial->first);
               if( it_result == result->end() )
               {
                    (*result)[it_monomial->first].SetDims( dpowern, dpowern);
                    it_result = result->find(it_monomial->first);
               }

               if( sum == n+1)
               {
                    (it_result->second)[0] += it_monomial->second * it_matrix->second[0];
               }
               else
               {
                   it_result->second += it_monomial->second * it_matrix->second;
               }
            }
        }

        for(it_result = result->begin(); it_result != result->end(); )
        {
            if( IsZero(it_result->second) )
            {
                result->erase(it_result++);
            }
            else
            {
                it_result++;
            }
        }
        if(!existQ && save_memory)
            reduction_matrix_ND_dict.erase(v);

        return reduction_matrix_ND_poly_dict.find(v);
    }
}

map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ>, vi64less> , vi64less>::const_iterator de_Rham_non_degenerate_local::compute_reduction_matrix_ND_poly_ZZ(const Vec<int64_t> v)
{
    map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ>, vi64less> , vi64less>::const_iterator it_ZZ;
    it_ZZ = reduction_matrix_ND_poly_ZZ_dict.find(v);
    if( it_ZZ != reduction_matrix_ND_poly_ZZ_dict.end() )
        return it_ZZ;
    else
    {
         if(verbose)
            cout<<"Computing the reduction matrix poly ND (ZZ) for v = "<<v<<endl;
        bool existQ;
        map< Vec<int64_t>, map< Vec<int64_t>, Mat<ZZ_p>, vi64less> , vi64less>::const_iterator it;
        map< Vec<int64_t>, Mat<ZZ_p>, vi64less>::const_iterator it_result;

        it = reduction_matrix_ND_poly_dict.find(v);
        existQ = (it != reduction_matrix_ND_poly_dict.end());
        if(!existQ)
            it = compute_reduction_matrix_ND_poly(v);

        //map< Vec<int64_t>, Mat<ZZ>, vi64less> result_zz;

        for(it_result = (it->second).begin(); it_result != (it->second).end(); it_result++)
            conv(reduction_matrix_ND_poly_ZZ_dict[v][it_result->first], it_result->second);

        if( !existQ && save_memory )
            reduction_matrix_ND_poly_dict.erase(v);

        //reduction_matrix_ND_poly_ZZ_dict[v] = result_zz;


        return reduction_matrix_ND_poly_ZZ_dict.find(v);
    }
}




Mat<ZZ_p>* de_Rham_non_degenerate_local::get_inclusion_matrix_ND(const Vec<int64_t> u, const int64_t l)
{
    int64_t i, sum;
    map< Vec<int64_t>, Vec< Mat<ZZ_p> >, vi64less >::iterator it;
    sum = 0;
    for(i = 0; i <= n; i++)
    {
        sum += u[i];
        assert(u[i] > 0);
    }
    assert( sum == d);
    assert( l < n );

    it = inclusion_matrix_ND_dict.find(u);

    if( it == inclusion_matrix_ND_dict.end() )
    {
        if( verbose )
            cout << "Computing inclusion matrix ND for u = "<<u<<endl;
        compute_inclusion_matrix_ND(u);
        it = inclusion_matrix_ND_dict.find(u);
    }

    return &( (it->second)[l] );
}



void de_Rham_non_degenerate_local::compute_inclusion_matrix_ND(const Vec<int64_t> u)
{
    assert( inclusion_matrix_ND_dict.find(u) == inclusion_matrix_ND_dict.end() );


    int64_t i, l, coordinate_of_monomial;
    Vec<int64_t> u_minus_ones;
    Vec< Mat<ZZ_p> >* result;
    pair< Vec<int64_t>, Mat<ZZ_p> >* solve;

    result = &inclusion_matrix_ND_dict[u];
    result->SetLength(n);

    u_minus_ones.SetLength(n+1);
    for(i = 0; i <= n ; i++)
        u_minus_ones[i] = u[i] - 1;

    for(l = 0; l < n ; l++)
    {
        /*
         * maps x^u l! Gl / f^(l+1) (x0 ... xn) -> l! H / f^(l+1) (the factorial is implicit)
         * deg Gl = l*d
         * and Gl \in coKernel of J_\empptyset at level l*d
         * deg H = (l+1) * d - (n + 1)
         */
        solve = get_solve_ND(l*d);
        (*result)[l].SetDims( tuple_list[(l + 1) * d - (n + 1)].length() , (solve->first).length() );
        for( coordinate_of_monomial = 0; coordinate_of_monomial < (int64_t) (solve->first).length(); coordinate_of_monomial++ )
            NTL::set( (*result)[l][ tuple_dict[ (l + 1) * d - (n + 1) ][ tuple_list[ l*d ][ (solve->first)[coordinate_of_monomial] ] + u_minus_ones ] ][coordinate_of_monomial] );
    }
}

Mat<ZZ_p>* de_Rham_non_degenerate_local::get_coKernels_ND_to_basis(const Vec<int64_t> u)
{
    int64_t i, sum;
    map< Vec<int64_t>, Mat<ZZ_p> , vi64less >::iterator it;

    sum = 0;
    for(i = 0; i <= n; i++)
    {
        sum += u[i];
        assert(u[i] > 0);
    }
    assert( sum == d);

    it = coKernels_ND_to_basis_dict.find(u);

    if( it  == coKernels_ND_to_basis_dict.end() )
    {
        if( verbose )
            cout << "Computing coKernels ND to basis for u = "<<u<<endl;
        compute_coKernels_ND_to_basis(u);
        it = coKernels_ND_to_basis_dict.find(u);
    }

    return &(it->second);
}


void de_Rham_non_degenerate_local::compute_coKernels_ND_to_basis(const Vec<int64_t> u)
{
    assert( coKernels_ND_to_basis_dict.find(u) == coKernels_ND_to_basis_dict.end() );
    /*
     * Computes the matrix from W_u to the cohomology basis
     * sum(u) = d and ui > 0
     *
     * For x^u * ( monomial of Gn ) Omega / x0 ... xn
     * write it as
     * x^a x^w where sum(w) = (d - 2) * (n + 1) + 1
     * hence x^w = \sum Fi di f
     * reduce x^a \sum Fi di f \Omega / f^(n+1) to H \Omega / f^n and apply final_reduction
     *
     * for all the rest use the inclusion and final_reductions
     */
    int64_t coordinate_of_monomial_w, coordinate_of_monomial_Fi, dim_Fi, i, j, l, pos;
    pair<Vec<int64_t>, Mat<ZZ_p> > *solve_top_J, *solve_n_ND;
    Mat<ZZ_p>* result;
    Mat<ZZ_p> tmp;
    Vec<ZZ_p> tmp_vec;
    Vec<int64_t> monomial, tt, a, w, ones_vector, ei;

    result = &(coKernels_ND_to_basis_dict[u]);
    result->SetDims( coKernels_J_basis.length(), coKernels_ND_basis.length() );

    ones_vector.SetLength(n+1);
    ei.SetLength(n+1);
    for(i = 0; i <= n; i++)
    {
        ones_vector[i] = 1;
        ei[i] = 0;
    }

    pos = 0;

    for(l = 0; l < n; l++)
    {
        /*
         * get_inclusion_matrix_ND(u, l)
         * x^u l! G_l/f^(l+1) (x0...xn) to -> l! H / f^(l+1) ( G_l \in W_u  and the factorial is implicit)
         * deg H = (l + 1) * d - (n + 1)
         *
         * get_final_reduction_matrix_J(l)
         * (l -1)! G \Omega / f^l to the cohomology basis elements
         */

        mul(tmp, *( get_final_reduction_matrix_J(l+1) ), *( get_inclusion_matrix_ND(u, l) ) );
        for( i = 0 ; i < (int64_t) tmp.NumRows(); i++)
        {
            for( j = 0; j < (int64_t) tmp.NumCols(); j++)
            {
                result->put(i,pos+j,tmp[i][j]);
            }
        }
        pos += tmp.NumCols();
    }

    //we take care of Gn
    solve_top_J = get_solve_J( (d - 2)*(n + 1) + 1 );
    solve_n_ND =  get_solve_ND( d * n );
    dim_Fi = tuple_list[ d * n - 2 * n ].length();

    for(l = 0 ; l < (int64_t) (solve_n_ND->first).length() ; l++)
    {
        Vec<ZZ_p> H, F;
        //first reduce the monomial to degree n
        monomial = tuple_list[d*n][ (solve_n_ND->first)[l] ];
        tt = tweak(monomial, n);
        a = monomial - tt;
        w = tt + u - ones_vector;
        //sum(w) = d*n - n + d - n - 1 = (d -2) * (n + 1) + 1

        H.SetLength( tuple_list[ d * n - (n + 1)].length() );
        F.SetLength( (solve_top_J->second).NumRows() );
        //F = U.column( dict_w[ immutable(w) ] )
        coordinate_of_monomial_w = tuple_dict[ (d - 2)*(n + 1) + 1 ][w];
        for(i = 0; i < (int64_t) F.length() ; i++)
            F[i] = (solve_top_J->second)[i][ coordinate_of_monomial_w ];

        for(i = 0; i <= n; i++)
        {
            ei[i]++;
            for( coordinate_of_monomial_Fi = 0; coordinate_of_monomial_Fi < dim_Fi; coordinate_of_monomial_Fi++)
            {
                if(a[i] + tuple_list[d * n - 2 * n][coordinate_of_monomial_Fi][i] >0)
                    H[ tuple_dict[ d * n - (n + 1) ][ a + tuple_list[d * n - 2 * n][coordinate_of_monomial_Fi] - ei] ] += (a[i] + tuple_list[d * n - 2 * n][coordinate_of_monomial_Fi][i]) * F[i*dim_Fi + coordinate_of_monomial_Fi];
            }
            ei[i]--;
        }
        mul(tmp_vec, *( get_final_reduction_matrix_J(n) ), H);
        for(i = 0; i < (int64_t) tmp_vec.length(); i++)
            result->put(i, pos, tmp_vec[i]);
        pos++;
    }
}

void de_Rham_non_degenerate_local::matrix_ND(Mat<ZZ_p>& result, int64_t l)
{
    /*
     * computes the matrix of the map
     *  (H0, ..., Hn) -> H0 * x0 * df/dx0 + ... + Hn * xn * df/dxn
     *  where deg Hi = l - d.
     */

    result.kill();
    if( l < d)
    {
        result.SetDims( tuple_list[l].length(), 0);
        return;
    }
    Vec< Vec<int64_t> >  *list_l, *list_lminusd;
    map< Vec<int64_t>, int64_t, vi64less > *dict_l;
    map< Vec<int64_t>, ZZ_p, vi64less >::const_iterator itf;
    int64_t i, j, len_lminusd;

    list_l = &tuple_list[l];
    list_lminusd = &tuple_list[l - d];
    dict_l = &tuple_dict[ l ];

    len_lminusd = list_lminusd->length();

    result.SetDims( list_l->length(), (n+1) * len_lminusd  );

    for(i = 0; i <= n ; i++)
    {
        for( j = 0 ; j < len_lminusd ; j++ )
        {
            for( itf = f.begin(); itf != f.end() ; itf++ )
            {
                result.put( (*dict_l)[ itf->first + (*list_lminusd)[j] ], i * len_lminusd + j, itf->first[i] * itf->second);
            }
        }
    }
}





Vec<ZZ_p> de_Rham_non_degenerate_local::monomial_to_basis_ND(Vec<int64_t> u)
{
    // Computes the coordinates of (m-1)! x^u \Omega / x0 .. xn f^m in the cohomology basis

    int64_t i, m, e, sum;
    Vec<ZZ_p> G, Gtmp;
    Vec<int64_t> src, dest, v;
    Mat<ZZ_p> M;

    sum = 0;
    for(i = 0; i <= n ; i++)
    {
        assert( u[i] > 0 );
        sum += u[i];
    }

    assert( sum % d == 0 );

    m = sum/d;

    G.SetLength( coKernels_ND_basis.length() );
    NTL::set(G[0]);

    v.SetLength(n+1);
    for( e = m; e > 1; e--)
    {
        for(i = 0; i <= n ; i++)
            v[i] = 0;

        sum = 0;

        while( sum < d)
        {
            for( i = 0; i <= n; i++ )
            {
                if( u[i] > v[i] + 1 )
                {
                    v[i]++;
                    sum++;
                    break;
                }
            }
        }

        u = u - v;
        reduce_vector_ND(Gtmp,u,v,G);
        swap(G,Gtmp);
        //M = get_reduction_matrix_ND(u, v);
        //mul(G, M, G);
    }
    sum = 0;
    for(i = 0; i <= n; i++)
        sum += u[i];

    assert( sum == d);

    return *(get_coKernels_ND_to_basis(u)) * G;
}


bool de_Rham_non_degenerate_local::test_monomial_to_basis_ND(int64_t N)
{
    if(verbose)
        cout<<"de_Rham_non_degenerate_local::test_monomial_to_basis_ND("<<N<<")"<<endl;
    Vec< map< Vec<int64_t>, ZZ_p, vi64less > > fpow;
    map< Vec<int64_t>, ZZ_p, vi64less>::iterator it1, it2;
    Vec<int64_t> u;
    Vec<ZZ_p> v;
    int64_t i, j, k, m, fact, sum;

    v.SetLength( coKernels_J_basis.length() );

    fpow.SetLength(N+1);

    u.SetLength(n+1);
    for( i = 0; i <= n; i++)
        u[i] = 0;

    NTL::set( fpow[0][u] );

    for( i = 0; i < N; i++)
    {
        for( it1 = fpow[i].begin() ; it1 != fpow[i].end() ; it1++ )
        {
            for( it2 = f.begin() ; it2 != f.end() ; it2++ )
            {
                u = it1->first + it2->first;
                fpow[i+1][u] += it1->second * it2->second;
            }
        }
    }
    for( i = 0; i <= N ; i++)
    {
        if(verbose)
            cout<<"i = " <<i<<endl;
        for( j = 0; j < (int64_t) coKernels_J_basis.length() ; j++ )
        {
            sum = 0;
            for(k = 0 ; k <= n; k++)
                sum += coKernels_J_basis[j][k];
            sum += (n + 1);
            m = sum/d + i;
            fact = 1;
            for( k = 1; k < m; k++)
            {
                fact *= k;
            }

            for( k = 0; k < (int64_t) coKernels_J_basis.length() ; k++ )
                v[k] = 0;

            for( it1 = fpow[i].begin(); it1 != fpow[i].end(); it1++ )
            {
                u = it1->first + coKernels_J_basis[j];
                for( k = 0 ; k <=n; k++ )
                    u[k]++;
                v += it1->second * monomial_to_basis_ND(u);
            }
            for( k = 0 ; k < (int64_t) coKernels_J_basis.length() ; k++ )
            {
                if ( k == j and v[k] != to_ZZ_p( fact ) )
                    return false;
                if ( k !=j and not IsZero(v[k]) )
                    return false;
            }
        }
    }
    return true;
}


bool de_Rham_non_degenerate_local::test_paths_ND(int64_t trials, int64_t paths)
{
    if(verbose)
        cout << "de_Rham_local::test_paths_J("<<trials<<", "<<paths<<")\n";
    assert(trials != 0);
    assert(paths != 0);

    int64_t attempt, pathi, i, k, sum, sum_u;
    Vec<int64_t> u_orig, u, v;
    Vec< Mat<ZZ_p> > M_saved;
    Mat<ZZ_p> M;
    Mat<ZZ_p> Mred;
    u_orig.SetLength(n+1);
    v.SetLength(n+1);

    for(attempt = 0; attempt < trials; attempt++)
    {
        if(verbose)
            cout << "attempt = "<<attempt<<endl;
        M_saved.kill();
        M_saved.SetLength(paths);
        sum = 0;
        for( i = 0 ; i <=n ; i++)
        {
            u_orig[i] = 1+RandomBnd(20);
            sum += u_orig[i];
        }

        while( sum % d != 0 )
        {
            k = RandomBnd(n+1);
            u_orig[k]++;
            sum++;
        }


        for( pathi = 0 ; pathi < paths ; pathi++)
        {
            if(verbose)
                cout << "path = "<<pathi <<endl;
            M.kill();
            M.SetDims( coKernels_ND_basis.length(), coKernels_ND_basis.length());
            for( i = 0; i < (int64_t) coKernels_ND_basis.length(); i++)
                NTL::set( M[i][i] );

            sum_u = 0;

            u = u_orig;
            for( i = 0; i < n+1; i++)
                sum_u += u[i];

            while( sum_u > d)
            {
                sum = 0;
                for( i = 0; i <= n; i++)
                    v[i] = 0;

                while( sum < d)
                {
                    k = RandomBnd(n + 1);
                    if( u[k] > v[k] + 1 )
                    {
                        v[k]++;
                        sum++;
                    }
                }
                u = u - v;
                sum_u -= d;

                Mred = get_reduction_matrix_ND(u, v);
                mul(M, Mred, M);
            }

            mul(M, *(get_coKernels_ND_to_basis(u)), M);

            M_saved[pathi] = M;
            if(pathi > 0)
                if(M_saved[pathi-1] != M_saved[pathi])
                    return false;
        }
    }
    return true;
}


bool isND(const map< Vec<int64_t>, zz_p, vi64less> &f)
{
    int64_t i, j;
    de_Rham_non_degenerate_local D;
    D = de_Rham_non_degenerate_local();
    D.precision = 1;
    D.verbose = false;
    map< Vec<int64_t>, zz_p, vi64less>::const_iterator it;
    it = f.begin();
    D.n = (it->first).length() - 1;
    D.d = 0;
    for(i = 0; i <=D.n ; i++)
        D.d += (it->first)[i];

    for(it = f.begin(); it != f.end() ; it++ )
        D.f[it->first] = rep( it->second);


    D.tuple_list.SetLength( (D.n + 1) * D.d + 1);
    D.tuple_dict.SetLength( (D.n + 1) * D.d + 1);
    for(i = 0 ; i < (D.n + 1) * D.d + 1 ; i++)
    {
        tuple_list_generator( D.tuple_list[i], i, (int64_t) (D.n+1));
        for( j = 0; j < (int64_t)D.tuple_list[i].length(); j++)
        {
            D.tuple_dict[i][ ((D.tuple_list[i])[j]) ] = j;
        }
    }

    Mat<ZZ_p> M;
    D.matrix_ND(M, (D.d-1)*(D.n+1) + 1 );
    Mat<zz_p> M_zz_p;
    M_zz_p.SetDims(M.NumRows(), M.NumCols());
    for(i = 0; i < (int64_t) M.NumRows(); i++)
    {
        for(j = 0; j < (int64_t) M.NumCols(); j++)
        {
            M_zz_p[i][j] = conv<zz_p>(rep(M[i][j]));
        }
    }
    nmod_mat_t M_flint;
    nmod_mat_init(M_flint,M.NumRows(), M.NumCols(), zz_p::modulus());
    long * P;
    bool result;
    P = (long *)flint_malloc(sizeof(long) * M.NumRows());
    conv(M_flint, M_zz_p);
    result = bool (M.NumRows() == nmod_mat_lu(P, M_flint, false));
    nmod_mat_clear(M_flint);
    flint_free(P);

    return result;
}


Mat<zz_p> find_change_of_variables( map< Vec<int64_t>, zz_p, vi64less> f, int64_t N)
{
    bool not_nd;
    int64_t i, n;
    Mat<zz_p> M;
    map< Vec<int64_t>, zz_p, vi64less> new_f;
    n = ((f.begin())->first).length();
    not_nd = !( isND(f) );
    ident(M, n);
    i = 0;
    while( i < N && not_nd)
    {
        M = random_SL_matrix(n);
        new_f = change_of_variables<zz_p>(f, M);
        not_nd = !( isND(new_f) );
        i++;
    }

    if(not_nd)
    {
        clear(M);
    }
    return M;
}



static void  nmod_mat_sub_submatrix(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, const int64_t max_row, const int64_t min_col)
{
    slong i;

    if (C->c == 0)
       return;

    for (i = 0; i < max_row; i++)
       _nmod_vec_sub(C->rows[i] + min_col, A->rows[i] + min_col, B->rows[i] + min_col, B->c - min_col, C->mod);
/*
    slong d = C->c - min_col;
    const slong left = max_row % 8;
    for(i = 0; i + 7 < max_row ; i+=8)
    {
        _nmod_vec_sub(C->rows[i] + min_col, A->rows[i] + min_col, B->rows[i] + min_col, d, C->mod);
        _nmod_vec_sub(C->rows[i + 1] + min_col, A->rows[i + 1] + min_col, B->rows[i + 1] + min_col, d, C->mod);
        _nmod_vec_sub(C->rows[i + 2] + min_col, A->rows[i + 2] + min_col, B->rows[i + 2] + min_col, d, C->mod);
        _nmod_vec_sub(C->rows[i + 3] + min_col, A->rows[i + 3] + min_col, B->rows[i + 3] + min_col, d, C->mod);
        _nmod_vec_sub(C->rows[i + 4] + min_col, A->rows[i + 4] + min_col, B->rows[i + 4] + min_col, d, C->mod);
        _nmod_vec_sub(C->rows[i + 5] + min_col, A->rows[i + 5] + min_col, B->rows[i + 5] + min_col, d, C->mod);
        _nmod_vec_sub(C->rows[i + 6] + min_col, A->rows[i + 6] + min_col, B->rows[i + 6] + min_col, d, C->mod);
        _nmod_vec_sub(C->rows[i + 7] + min_col, A->rows[i + 7] + min_col, B->rows[i + 7] + min_col, d, C->mod);
    }
    switch (left)
    {
        case 7 : _nmod_vec_sub(C->rows[i + 6] + min_col, A->rows[i + 6] + min_col, B->rows[i + 6] + min_col, d, C->mod);
        case 6 : _nmod_vec_sub(C->rows[i + 5] + min_col, A->rows[i + 5] + min_col, B->rows[i + 5] + min_col, d, C->mod);
        case 5  : _nmod_vec_sub(C->rows[i + 4] + min_col, A->rows[i + 4] + min_col, B->rows[i + 4] + min_col, d, C->mod);
        case 4  : _nmod_vec_sub(C->rows[i + 3] + min_col, A->rows[i + 3] + min_col, B->rows[i + 3] + min_col, d, C->mod);
        case 3  : _nmod_vec_sub(C->rows[i + 2] + min_col, A->rows[i + 2] + min_col, B->rows[i + 2] + min_col, d, C->mod);
        case 2  : _nmod_vec_sub(C->rows[i + 1] + min_col, A->rows[i + 1] + min_col, B->rows[i + 1] + min_col, d, C->mod);
        case 1  : _nmod_vec_sub(C->rows[i] + min_col, A->rows[i] + min_col, B->rows[i] + min_col, d, C->mod);
        case 0 :  ;
    }

    if(C->mod.norm)
    {
        switch(d)
        {
            case 64:
                for (i = 0; i < max_row; i++)
                    _NMOD_VEC_SUB64(C->rows[i] + min_col,  A->rows[i] + min_col, B->rows[i] + min_col, C->mod.n);
                break;

            case 63:
                for (i = 0; i < max_row; i++)
                    _NMOD_VEC_SUB63(C->rows[i] + min_col,  A->rows[i] + min_col, B->rows[i] + min_col, C->mod.n);
                break;

            case 32:
                for (i = 0; i < max_row; i++)
                    _NMOD_VEC_SUB32(C->rows[i] + min_col,  A->rows[i] + min_col, B->rows[i] + min_col, C->mod.n);
                break;

            case 1:
                for (i = 0; i < max_row; i++)
                    _NMOD_VEC_SUB1(*(C->rows[i] + min_col),  *(A->rows[i] + min_col), *(B->rows[i] + min_col), C->mod.n);
                break;

            default:
                for (i = 0; i < max_row; i++)
                    _nmod_vec_sub(C->rows[i] + min_col, A->rows[i] + min_col, B->rows[i] + min_col, C->c - min_col, C->mod);
                break;
        }
    }
    else
    {
        switch(d)
        {
            case 64:
                for (i = 0; i < max_row; i++)
                    NMOD_VEC_SUB64(C->rows[i] + min_col,  A->rows[i] + min_col, B->rows[i] + min_col, C->mod.n);
                break;

            case 63:
                for (i = 0; i < max_row; i++)
                    NMOD_VEC_SUB63(C->rows[i] + min_col,  A->rows[i] + min_col, B->rows[i] + min_col, C->mod.n);
                break;

            case 32:
                for (i = 0; i < max_row; i++)
                    NMOD_VEC_SUB32(C->rows[i] + min_col,  A->rows[i] + min_col, B->rows[i] + min_col, C->mod.n);
                break;

            case 1:
                for (i = 0; i < max_row; i++)
                    NMOD_VEC_SUB1(*(C->rows[i] + min_col),  *(A->rows[i] + min_col), *(B->rows[i] + min_col), C->mod.n);
                break;

            default:
                for (i = 0; i < max_row; i++)
                    _nmod_vec_sub(C->rows[i] + min_col, A->rows[i] + min_col, B->rows[i] + min_col, C->c - min_col, C->mod);
                break;
        }
    }
    */
}









/*
 * Input:
 * - k, a positive integer > n + 1
 * - G, a vector
 * - poly = M, a vector with (n + 1) matrices, representing  M(Y) = M_0 + M_1 * Y + M_2 Y^2 + ... + M_ * Y^n + M_{n+1} *Y^(n+1)
 *
 * Output:
 * - H = M(0) M(1) ... M(k-1) G
*/
void finitediff_flint_nmod(fmpz * result, fmpz_mat_struct * M_fmpz, const int64_t Mlength, const int64_t k, const fmpz * G, const fmpz_t &modulus)
{
    int64_t l, i, j;
    int64_t n = Mlength - 2;
    int64_t d = fmpz_mat_ncols(M_fmpz + 0);
    assert(k > n + 1);


    nmod_t mod;
    int nlimbs;
    nmod_init(&mod, fmpz_get_ui(modulus));
    mp_ptr H, Hout;
    nlimbs = _nmod_vec_dot_bound_limbs(d, mod);


    H =  _nmod_vec_init(d);
    for(i = 0; i < d; i++)
    {
        H[i] = fmpz_get_ui(G + i);
    }
    Hout = _nmod_vec_init(d);


    nmod_mat_struct * M = new nmod_mat_struct[Mlength];

    // we usually deal with matrices that the only nonzero entries are on the top right corner
    // we should take advantage of that
    int64_t * max_row = new int64_t[Mlength];
    int64_t * min_col = new int64_t[Mlength];
    // recall Mlength = n + 2
    for(i = 0; i < n + 2; i++)
    {
        max_row[i] = 0;
        min_col[i] = d;

        nmod_mat_init(M + i, d, d, mod.n);
        for(j = 0; j < d; j++)
        {
            for(l = 0; l < d; l++)
            {
                nmod_mat_entry(M + i, j, l) = fmpz_get_ui( fmpz_mat_entry(M_fmpz + i, j, l));
                if(nmod_mat_entry(M + i, j, l) != 0)
                {
                    if(max_row[i] <= j)
                        max_row[i] = j + 1;
                    if(min_col[i] > l)
                        min_col[i] = l;
                }
            }
        }
    }
    i = n + 1;
    while( max_row[i] == 0 and min_col[i] == d and i >= 0)
    {
        //M[i] = 0;
        n--;
        i--;
    }




    // Mfd[l] = M(k - 1 - (self.n + 1) + l)
    // equivalently
    // Mfd[n + 1 -l] = M(k - 1 - l)
    // Mfd = [ sum( ( self.R(k - 1 - (self.n + 1) + l) ** i) * Mi for i, Mi in  enumerate(M) ) for l in range(self.n + 2) ];
    nmod_mat_struct * Mfd = new nmod_mat_struct[n + 2];
    for(l = 0; l < n + 2; l++)
    {
        //compute initial table of differences
        nmod_mat_init(Mfd + l, d, d, mod.n);
        for(i = 0; i < n + 2; i++)
        {
            mp_limb_t tmp = nmod_pow_ui(k - 1 - (n + 1) + l, i, mod);
# if __FLINT_RELEASE > 20700
            nmod_mat_scalar_addmul_ui(Mfd + l, Mfd + l, M + i, tmp);
# else
            nmod_mat_scalar_mul_add(Mfd + l, Mfd + l, tmp, M + i);
# endif
        }
    }
    for(l = 0; l < n + 2; l++)
    {
        //u v^{(k - 1 - l)  + 1} --> u v^{(k - 1 - l)}
        // Mfd[n + 1 -l] = M(k - 1 - l)
        mul(Hout, Mfd + n + 1 - l, H, nlimbs);
        swap(H, Hout);
    }
    //
    // make Mfd[l] = M[a, a - 1, ..., a - l]
    // where a = k - 1 - (n + 1);
    for(l = 1; l < n + 2; l++)
        for(j = n + 1; j >= l; j--)
            nmod_mat_sub(Mfd + j, Mfd + j, Mfd + j - 1);

    for(l = 0; l < (k - 1 - (n + 1)); l++)
    {
        // Mfd[0] =  M(k - 1 - (n + 1) - l)
        // update Mfd vector
        //
        // deg(M) = n + 1 ==> Mfd[n+1] is constant
        //
        // in the setting of controlled reduction for projective hypersurfaces
        // only the first row of M can have degree n + 1 (indeed we only really need to subtract the last entries)
        for(j = n; j >= 0; j--)         {
            //nmod_mat_sub(Mfd + j, Mfd + j, Mfd + j + 1);
            nmod_mat_sub_submatrix(Mfd + j, Mfd + j,  Mfd + j + 1, max_row[j+1], min_col[j + 1]);
            //cout << "j = "<< j << endl;
            //nmod_mat_print_pretty(Mfd + j);
            //cout << endl;
        }
        // after
        // Mfd[0] =  M(k - 1 - (n + 1) - l - 1)
        mul(Hout, Mfd + 0, H, nlimbs);
        swap(H, Hout);
    }
    //check that we computed M[0] in S
    assert( nmod_mat_equal(Mfd + 0, M + 0));

    for(i = 0; i < Mlength; ++i)
        nmod_mat_clear(M + i);
    for(i = 0; i < n + 2; ++i)
        nmod_mat_clear(Mfd + i);
    delete[] Mfd;
    delete[] M;

    for(i = 0; i < d; i++)
        fmpz_set_ui(result + i, H[i]);
    _nmod_vec_clear(Hout);
    _nmod_vec_clear(H);

    delete[] max_row;
    delete[] min_col;
}





