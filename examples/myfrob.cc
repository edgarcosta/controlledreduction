// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// frob.cpp: computes frob and it's charpoly, and possibly the picard rank

/*
 * Input file:
 * p
 * precision
 * n
 * d
 * f_vector
 * N
 * absolute precision
 *
 * Output file:
 * p
 * precision
 * n
 * d
 * f_original
 * N
 * absolute precision
 * change of variables
 * f
 * Frob
 * charpoly_corrected
 * wall time
 * user time
 */
#include "hypersurface_nd.h"
#include "matrix.h"
#include "tools.h"
#include "timing.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        cout<<"Needs 3 arguments: "<<argv[0]<<" #threads input_file output_file "<<endl;
        cout<<"Input file format (for a generic hypersurface):"<<endl;
        cout<<"\tp\n\tprecision\n\tn\n\td\n\tf = [a_0 ... a_k] with f = a0 * x_n ^ d + a1 * x_{n-1}*x_n^(d-1) + a2* x_{n-1} ^2 * x_n^(d-2) + ... + ad x_{n-1} ^d + ...\n";
        cout <<"\tN = [N_1 ... N_n]\n\tabsolute precision of the characteristic polynomial"<<endl<<endl;
        cout<<"OR for a K3 surface, ie, n = 3 and d =4:\n";
        cout<<"\tp\n\tf = [a_0 ... a_34] with f = a0 * w^4 + a1 * z*w^3 ...\n";
        cout<<endl;
        cout<<"Output file format:"<<endl;
        cout<<"\tp\n\tprecision\n\tn\n\td\n\tf_original\n\tN\n\tabsolute precision of the characteristic polynomial\n";
        cout<<"\tchange of variables (st f is non degenerate)\n\tf\n";
        cout<<"\tFrob matrix\n\tcharacteristic polynomial of Frob\n\tgeometric picard rank (for  n = 3 and d = 4) otherwise 0\n";
        cout<<"\twall time\n\tuser time\n";
        return 0;
    }
    SetSeed(to_ZZ(42));
#ifdef _OPENMP
    omp_set_num_threads(atoi(argv[1]));
#endif
    int64_t p, precision, n , d;
    int64_t i, numberoflines;
    timestamp_type wtime1, wtime2;
    get_timestamp(&wtime1);
    double wall_time, user_time, user_time1, user_time2;
    user_time1 = get_cpu_time();
    Vec<int64_t> N;
    Vec<int64_t> absolute_precision;
    Mat<ZZ> Frob_ZZ;
    Vec<ZZ> cp;
    bool is_ND;
    




    ifstream file;
    file.open(argv[2]);
    if(! file.is_open())
    {
        cout << "Could not open \'"<<argv[2]<<"\'"<<endl;
        abort();
    }

    ifstream file2;
    file2.open(argv[2]);
    numberoflines = count(istreambuf_iterator<char>(file2), istreambuf_iterator<char>(), '\n');
    file2.close();

    file >> p;

    int64_t ND_tries = p*100;
    if(argc > 4)
      ND_tries = atoi(argv[4]);

    stringstream buffer;
    if(numberoflines  < 7)
    {
        n = 2;
        d = 4;
        if( p == 3)
        {
            precision = 7;
            buffer << "[4 4]\n[5 4 3 3 3 3 2]";
        }
        else if( p == 5)
        {
            precision = 5;
            buffer << "[3 3]\n[5 4 3 3 3 3 2]";

        }
        else if( 5 < p && p < 16)
        {
            precision = 5;
            buffer << "[3 3]\n[5 4 3 3 3 3 3]";
        }
        else if( 16 < p && p < 300)
        {
            precision = 3;
            buffer << "[2 2]\n[4 3 3 2 2 2 2]";
        }
        else
        {
            assert(p >= 300);
            precision = 3;
            buffer << "[0 1]\n[4 3 2 2 2 2 2]";

        }
    }
    else
    {
        file >> precision;
        file >> n;
        file >> d;
    }

    zz_p::init(p);
    ZZ_p::init( power_ZZ(p,precision) );

    hypersurface_non_degenerate hs_ND;
    hypersurface hs;

    Vec<zz_p> f_original_vector;
    map< Vec<int64_t>, zz_p, vi64less> f_original_map;
    Vec<zz_p> f_vector;
    map< Vec<int64_t>, zz_p, vi64less> f_map;
    Mat<zz_p> M;
    Mat<ZZ_p> Frob;

    file >> f_original_vector;
    if(numberoflines < 7)
    {
        buffer >> N;
        buffer >> absolute_precision;
    }
    else
    {
        file >> N;
        file >> absolute_precision;
    }

    file.close();

    Vec< Vec<int64_t> > tuple_d;
    map< Vec<int64_t>, int64_t, vi64less>  tuple_d_dict;
    tuple_list_generator(tuple_d, d, (int64_t) (n +1) );
    assert(tuple_d.length() == f_original_vector.length() );

    for( i = 0; i < (int64_t) f_original_vector.length() ; i++)
    {
        if( !IsZero(f_original_vector[i]) )
            f_original_map[ tuple_d[i] ] = f_original_vector[i];
        tuple_d_dict[ tuple_d[i] ] = i;
    }

    if(!isSmooth(f_original_map))
    {
        cout << "f is not smooth!" <<endl;
        abort();
    }
    M = find_change_of_variables(f_original_map, ND_tries);
    is_ND = !IsZero(M);
    if( is_ND )
    {
        cout<<"Found a change of variables!"<<endl;
        f_map = change_of_variables<zz_p>(f_original_map, M);
        hs_ND = hypersurface_non_degenerate(p, precision, f_map, true);
        assert(hs_ND.dR->coKernels_J_basis.length() + 1 == absolute_precision.length() );
    }
    else
    {
        cout<<"Wasn't able to find a suitable change of variables to make it non degenerate!"<<endl;
        f_map = f_original_map;
        hs = hypersurface(p, precision, f_map, true);
        assert(hs.dR->coKernels_J_basis.length() + 1 == absolute_precision.length() );
    }


    f_vector.SetLength(f_original_vector.length());
    map< Vec<int64_t>, zz_p, vi64less>::const_iterator fit;
    for( fit = f_map.begin(); fit != f_map.end(); fit++)
        f_vector[ tuple_d_dict[ fit->first ] ] = fit->second;

    if(is_ND)
        Frob = hs_ND.frob_matrix_ND(N);
    else
        Frob = hs.frob_matrix_J(N);

    cout << "Frob = "<<Frob<<endl;

    Frob_ZZ = conv<Mat<ZZ> >(Frob);

    cp = charpoly_frob(Frob_ZZ, absolute_precision, p, n - 1);

    cout <<"Characteristic polynomial = "<< cp <<endl;

    size_t rank;
    if(n == 3 && d == 4)
    {
        rank = geometric_picard(cp, p, n - 1);
        cout <<"Geometric rank = "<<rank<<endl;
    }
    else
    {
        rank = 0;
    }

    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time2 = get_cpu_time();
    user_time = user_time2 - user_time1;
    /*
     * Output file:
     * p
     * precision
     * n
     * d
     * f_original
     * N
     * absolute precision
     * change of variables
     * f
     * Frob
     * charpoly_corrected
     * rank
     * wall time
     * user time
     */
    ofstream output;
    stringstream ss;

    //ss << p <<endl;
    //ss << precision <<endl;
    //ss << n <<endl;
    //ss << d <<endl;
    //ss << f_original_vector <<endl;
    //ss << N <<endl;
    //ss << absolute_precision <<endl;
    //ss << M <<endl;
    //ss << f_vector <<endl;
    //ss << Frob <<endl;
    ss << cp <<endl;
    //ss << rank <<endl;
    //ss << wall_time <<endl;
    //ss << user_time <<endl;


    output.open(argv[3]);
    if(output.is_open())
    {
        output << ss.str() ;
        output.close();
        cout << "Output file saved on \'"<<argv[3]<<"\'"<<endl;
    }
    else
    {
        cout << "Could not open \'"<<argv[3]<<"\'"<<endl; 
    }
    cout << ss.str() ;
    cout << "EOF"<<endl;

    /*
    //cleaning memory
    if(is_ND)
    {
        delete hs_ND.dR_ND;
    }
    else
    {
        delete hs.dR;
    }
    */

    return 0;
}
