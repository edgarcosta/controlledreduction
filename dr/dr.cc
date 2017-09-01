// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// dr.cpp: routines for the class de_Rham_local
   
#include "dr.h"
//constructors

de_Rham_local::de_Rham_local(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose, bool save_memory)
{
    init(p, precision, fbar, verbose, save_memory);
}

de_Rham_local::de_Rham_local(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<zz_p> fbar_vector, bool verbose, bool save_memory)
{
    Vec< Vec<int64_t> > tuple_d;
    map< Vec<int64_t>, zz_p, vi64less> fbar_map;
    int64_t i;
    tuple_list_generator(tuple_d, d, (int64_t) (n + 1) );
    assert(tuple_d.length() == fbar_vector.length() );
    for( i = 0; i < (int64_t) fbar_vector.length() ; i++)
        fbar_map[ tuple_d[i] ] = fbar_vector[i];

    init(p, precision, fbar_map, verbose, save_memory);

}

de_Rham_local::de_Rham_local(int64_t p, int64_t precision, int64_t n, int64_t d, bool verbose, bool save_memory)
{
    int64_t i ,j;
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
        solve_J.clear();
        for( i = 0; i < (int64_t) tuple_list[d].length(); i++ )
        {
            x = random_zz_p();
            fbar[ tuple_list[d][i] ] = x;
            f[ tuple_list[d][i] ] = rep( x );
        
        }
        B = (get_solve_J((d-2)*(n+1)+1))->first;
        boolean = bool( B.length() != 0 );
    }

    solve_J.clear();
    this->precision = precision;
    
    init(p, precision, fbar, verbose, save_memory);

    
}

de_Rham_local::de_Rham_local(const char* filename)
{
    int64_t i, j;
    ifstream file;
    Vec<zz_p> fbar_vector;
    pair< Vec<int64_t>, Mat<ZZ_p> > *solve_pair;
    file.open(filename);

    if(file.is_open())
    {
        file >> p;
        assert( p == (int64_t)  zz_p::modulus() );
        file >> precision;
        assert(  power_ZZ(p, precision) == ZZ_p::modulus() );
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
                append(coKernels_J_basis, tuple_list[ (i + 1) * d - (n + 1) ][ (solve_pair->first)[j] ]);
            }
        }
        solve_pair = &( solve_J[(d - 2) * (n + 1) + 1]);
        file >> solve_pair->first;
        file >> solve_pair->second;

        for( i = 0; i < (int64_t)coKernels_J_basis.length(); i++ )
        {
            coKernels_J_basis_dict[ coKernels_J_basis[i] ] = i;
        }
        
    }
    else
    {
        cout << "Could not open \'"<<filename<<"\'"<<endl;
        abort();
    }
}



void de_Rham_local::init(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose, bool save_memory)
{
    int64_t i, j;
    ZZX pol, H;
    pair< Vec<int64_t>, Mat<ZZ_p> >* solve_pair;
    
    map< Vec<int64_t>, zz_p, vi64less>::iterator it;
    

    this->p = p;
    this->precision = precision;
    this->verbose = verbose;
    this->save_memory = save_memory;
    this->fbar = fbar;
    for(it = fbar.begin(); it != fbar.end(); it++)
    {
        f[it->first] = rep(it->second);
    }

    it = fbar.begin();
    n = (int64_t) (it->first).length() - 1;
    d = 0;
    for(i = 0; i <= n; i++)
        d += (it->first)[i];

    if(verbose)
    {
        cout<<"n = "<<n;
        cout<<" d = "<<d;
        cout<<" p = "<<p;
        cout<<" precision = "<<precision;
        cout<<endl;
        cout <<"fbar = \n";
        cout <<= fbar;
        cout <<endl;
    }

    
    assert(d>n);

    assert(d%p != 0);
    
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

    // Hilbert_J = (1 + x + x^2 + ... + x^(d-2) ) ^(n+1)
    Hilbert_J.SetLength((d - 2) * (n + 1) + 1);

    for(i = 0; i < d -1; i++)
        SetCoeff(pol, i);

    H = 1;

    for(i = 0; i< n + 1; i++)
        mul(H, H, pol);

    for(i = 0; i <= (d - 2)*(n + 1) ; i++)
        Hilbert_J[i] = (int64_t) to_ulong(coeff(H, i));

    if(verbose)
    {
        cout<<"Hilbert_J = "<<Hilbert_J<<endl;
    }

    if(verbose)
        cout << "Asserting that f is smooth.\n";

    solve_pair = get_solve_J((d - 2)*(n + 1) + 1);
    assert((solve_pair->first).length() == 0);
 
    for( i = 0 ; i < n; i++)
    {
        solve_pair = get_solve_J((i+1)*d-(n+1));
        assert( Hilbert_J[(i+1)*d-(n+1)] == (int64_t) (solve_pair->first).length() );
        for( j = 0; j < (int64_t)(solve_pair->first).length(); j++)
        {
            append(coKernels_J_basis, tuple_list[ (i + 1) * d - (n + 1) ][ (solve_pair->first)[j] ]);
        }

    }
    for( i = 0; i < (int64_t)coKernels_J_basis.length(); i++ )
    {
        coKernels_J_basis_dict[ coKernels_J_basis[i] ] = i;
    }    

}

bool de_Rham_local::save(const char * filename)
{
    if(verbose)
        cout <<"Saving de_Rham_local to: "<< filename << endl;

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

    

    if(file.is_open() && fileswp.is_open())
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
            solve_pair = get_solve_J( (i+1)*d - (n+1) );
            file << solve_pair->first <<endl;
            file << solve_pair->second <<endl;
        }
        solve_pair = get_solve_J((d - 2)*(n + 1) + 1);
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

void de_Rham_local::compute_everything_J()
{
    int64_t i;
    Vec<int64_t> u;
    for(i = 0; i < (int64_t) tuple_list[d].length(); i++)
    {
        u = tuple_list[d][i];
        compute_reduction_matrix_J( u );
        compute_reduction_matrix_J_ZZ( u );
    }

    for(i = 0; i  < (int64_t) tuple_list[n].length(); i++)
    {
        u = tuple_list[n][i];
        get_inclusion_matrix_J(u);
    }
    get_final_reduction_matrix_J(n); 
}

pair< Vec<int64_t>, Mat<ZZ_p> >* de_Rham_local::get_solve_J(int64_t level)
{
    map< int64_t, pair< Vec<int64_t>, Mat<ZZ_p> > >::iterator it;
    it = solve_J.find(level);
    if( it != solve_J.end() )
        return &(it->second);
    
    
    Mat<ZZ_p> MJ;
    matrix_J(MJ, level);

    if(verbose)
        cout << "Computing and solving matrix of relations at degree = "<<level<<" ( "<<MJ.NumRows()<<"x"<<MJ.NumCols()<<" )."<<endl;
    solve_system_padic(solve_J[level].first, solve_J[level].second, MJ, precision);
    return &(solve_J[ level ]);
}

Mat<ZZ_p> de_Rham_local::get_reduction_matrix_J(const Vec<int64_t> u, const Vec<int64_t> v)
{
    int64_t i , sum;
    Mat<ZZ_p> result;
    sum = 0;
    for(i = 0; i <= n; i++)
    {
        sum += v[i];
        assert(u[i]==0 || v[i] > 0 );
    }
    assert( sum == d);
    
    Mat<ZZ_p> temp;
    map< Vec<int64_t> ,  Vec<Mat<ZZ_p> >, vi64less>::const_iterator it;
    it = compute_reduction_matrix_J(v);
    result = it->second[0];
    for( i = 0; i <= n ; i++)
    {
        mul(temp, it->second[i+1], u[i]);
        add(result,result,temp);
    }
    return result;
}

Mat<ZZ> de_Rham_local::get_reduction_matrix_J_ZZ(const Vec<int64_t> u, const Vec<int64_t> v)
{
    int64_t i , sum;
    Mat<ZZ> result;
    sum = 0;
    for(i = 0; i <= n; i++)
    {
        sum += v[i];
        assert(u[i]==0 || v[i] > 0 );
    }
    assert( sum == d);
    
    map< Vec<int64_t>, Vec<Mat<ZZ> >, vi64less>::const_iterator it;
    it =  compute_reduction_matrix_J_ZZ(v);
    result = it->second[0];
    for( i = 0; i <= n ; i++)
    {
        result += it->second[i+1]*u[i];
    }
    return result;
}



void de_Rham_local::reduce_vector_J_ZZ(Vec<ZZ> &result, const Vec<int64_t> u, const Vec<int64_t> v, const Vec<ZZ> G)
{
    int64_t i, sum;
    sum = 0;
    for(i = 0; i <= n; i++)
    {
        sum += v[i];
        assert(u[i]==0 || v[i] > 0 );
        //u_zz_p[i] = conv<ZZ_p>(u[i]);
    }
    assert( sum == d);

        
    map< Vec<int64_t> ,  Vec<Mat<ZZ> >, vi64less>::const_iterator it;
    it = compute_reduction_matrix_J_ZZ(v);
    
    result = it->second[0] * G;
    for( i = 0; i <= n ; i++)
    {
        result += u[i] * (it->second[i+1] * G);
    }
}






Vec<ZZ_p> de_Rham_local::reduce_vector_J(const Vec<int64_t> u, const Vec<int64_t> v, const Vec<ZZ_p> G)
{
    int64_t i, sum;
    Vec<ZZ_p> result;
    sum = 0;
    for(i = 0; i <= n; i++)
    {
        sum += v[i];
        assert(u[i]==0 || v[i] > 0 );
        //u_zz_p[i] = conv<ZZ_p>(u[i]);
    }
    assert( sum == d);

        
    map< Vec<int64_t> ,  Vec<Mat<ZZ_p> >, vi64less>::const_iterator it;
    it = compute_reduction_matrix_J(v);
    result = it->second[0] * G;
    for( i = 0; i <= n ; i++)
    {
        result += u[i] * (it->second[i+1] * G);
    }
    return result;
}


void de_Rham_local::reduce_vector_J_poly(Vec<ZZ_p> &result, const Vec<int64_t> u, const Vec<int64_t> v, const int64_t iterations, const Vec<ZZ_p> G)
{
    int64_t i, sum;
    sum = 0;
    for(i = 0; i <= n; i++)
    {
        sum += v[i];
        assert(u[i]==0 || v[i] > 0 );
        //u_zz_p[i] = conv<ZZ_p>(u[i]);
    }
    assert( sum == d);
    map< Vec<int64_t> ,  Vec<Mat<ZZ_p> >, vi64less>::const_iterator it;
    it = compute_reduction_matrix_J(v);

    Mat<ZZ_p>  M0, M1;
    Vec<ZZ_p> Gin,  Gout0, Gout1;
    //Vec<ZZ_p> Gout;
    int64_t x;

    M0 = it->second[0];
    M1.SetDims( G.length(), G.length() );
    for( i = 0; i <= n ; i++)
    {
        M0 += u[i] * it->second[i+1];
        M1 += v[i] * it->second[i+1];
    }

    Gin = G;
    for(x = iterations - 1; x != (int64_t)-1 ; x--) //couting on overflow
    {
        //Gout = M0 * Gin;
        //Gout += x * (M1 * Gin);
        mul(Gout0, M0, Gin);
        mul(Gout1, M1, Gin);
        mul(Gout1,x,Gout1);
        add(Gin,Gout0,Gout1);
        //swap(Gin,Gout);
    }
    result = Gin;
}



void de_Rham_local::reduce_vector_J_poly_ZZ(Vec<ZZ> &result, const Vec<int64_t> u, const Vec<int64_t> v, const int64_t iterations, const Vec<ZZ> G)
{
    int64_t i, sum;
    sum = 0;
    for(i = 0; i <= n; i++)
    {
        sum += v[i];
        assert(u[i]==0 || v[i] > 0 );
        //u_zz_p[i] = conv<ZZ_p>(u[i]);
    }
    assert( sum == d);
    map< Vec<int64_t> ,  Vec<Mat<ZZ> >, vi64less>::const_iterator it;
    it = compute_reduction_matrix_J_ZZ(v);

    Mat<ZZ> M0, M1;
    Vec<ZZ> Gin, Gout, Gout0, Gout1;
    int64_t x;

    M0 = it->second[0];
    M1.SetDims( G.length(), G.length() );
    for( i = 0; i <= n ; i++)
    {
        M0 += u[i] * it->second[i+1];
        M1 += v[i] * it->second[i+1];
    }

    Gin = G;
    for(x = iterations - 1; x != (int64_t)-1; x--)//counting on overflow
    {
        mul(Gout0, M0, Gin);
        mul(Gout1, M1, Gin);
        mul(Gout1,x,Gout1);
        add(Gout,Gout0,Gout1);
        for(i = 0; i < (int64_t) G.length(); i++)
            rem(Gin[i], Gout[i], ZZ_p::modulus() );
    }
    result = Gin;
}


/*
 * computes the reduction matrix from V_{u+v} to V_u as a polynomial in u0,...,un
 * and automatically adds it to the dictionary
 * Recall V_u = {  (m-1)! x^u G Omega / x0 ... xn  f^m}
 * G dense polynomial of degree d*n-n and x0...xn | x^u G
 * and we assume that u_i == 0 for every i such that v_i == 0
*/
map< Vec<int64_t>, Vec<Mat<ZZ_p> >, vi64less>::const_iterator de_Rham_local::compute_reduction_matrix_J(const Vec<int64_t> v)
{
    map< Vec<int64_t> ,  Vec<Mat<ZZ_p> >, vi64less>::const_iterator it;
    it = reduction_matrix_J_dict.find(v);
    if(it != reduction_matrix_J_dict.end())
    {
        return it;
    }
    else
    {
        if(verbose)
            cout<<"Computing the reduction matrix J for v = "<<v<<endl;

        Mat<ZZ_p>* solve_top;
        Vec< Vec<int64_t> >* list_G;
        Vec< Vec<int64_t> >* list_F;
        map< Vec<int64_t>, int64_t, vi64less>* dict_G;
        map< Vec<int64_t>, int64_t, vi64less>* dict_w;
        Vec<int64_t> w, one_minus_ei;
        Vec< Mat<ZZ_p> >* M;
        Vec<ZZ_p> F;
        int64_t i, dim_Fi, coordinate_of_monomial, coordinate_of_w, coordinate_of_monomial_of_Fi, row;
        bool boolean;
        

        solve_top = &( ( get_solve_J( (d - 2)*(n + 1) + 1) )->second );

        F.SetLength(solve_top->NumRows());

        list_G = &tuple_list[ d * n - n ];
        dict_G = &tuple_dict[ d * n - n ];
        dict_w = &tuple_dict[ (d * n - n) + d - (n+1) ]; //deg G + sum(v) - (n+1) = (d-2)*(n+1) + 1
        list_F = &tuple_list[ d * n - 2 * n ]; //deg G + d -(n+1) - (d-1) = (d-2)*(n+1) + 1 - ( d - 1)

        
        dim_Fi = list_F->length();

        M = &(reduction_matrix_J_dict[v]);
        M->SetLength(n+2);

        for(i = 0; i < n + 2; i++)
        {
            (*M)[i].SetDims(list_G->length(),list_G->length());
        }

        one_minus_ei.SetLength(n+1);
        for(i = 0; i<= n ; i++)
            one_minus_ei[i] = 1;

        for( coordinate_of_monomial = 0; coordinate_of_monomial < (int64_t) list_G->length() ; coordinate_of_monomial++ )
        {
            w = (*list_G)[coordinate_of_monomial] + v; // sum(w) = d * n - n + d
            boolean = true;
            for(i = 0; i <= n; i++)
            {
                if( w[i] > 0)
                    w[i]--;
                else
                {
                    boolean = false;
                    break;
                }
            }
            if(boolean) // sum(w) = d * n - 2n + d - 1
            {
                int64_t sum=0;
                for( i =0 ; i <= n ; i++)
                    sum += w[i];
                assert(sum == d * n - 2*n + d - 1);
                coordinate_of_w = (*dict_w)[w];
                // F = solve_top.column(coordinate_of_w)
                for( i = 0; i < (int64_t) solve_top->NumRows(); i++)
                {
                    F[i] = solve_top->get(i, coordinate_of_w);
                }
                for( i = 0; i <= n; i++)
                {
                    one_minus_ei[i]--;
                    for(coordinate_of_monomial_of_Fi = 0; coordinate_of_monomial_of_Fi < dim_Fi; coordinate_of_monomial_of_Fi++)
                    {
                        // x^(u+v) * G / x0 ... xn = \sum c_w x^u * x^w
                        // x^u * x^w = x^u \sum F_i di f ~ \sum di( x^u F_i )
                        // x^u * x^s * di f ~ (u[i]+s[i]) x^u * x^s / xi
                        // (u[i]+s[i]) * x^u * x^s x0 * ... * xi-1 * xi+1 * ... xn / x0 ... xn
                        row = (*dict_G)[ (*list_F)[coordinate_of_monomial_of_Fi] + one_minus_ei ];
                        (*M)[0][row][coordinate_of_monomial] += (*list_F)[coordinate_of_monomial_of_Fi][i] * F[i*dim_Fi + coordinate_of_monomial_of_Fi];
                        (*M)[i+1][row][coordinate_of_monomial] += F[i*dim_Fi + coordinate_of_monomial_of_Fi];
                    }
                    one_minus_ei[i]++;
                }
            }
        }
        return reduction_matrix_J_dict.find(v);
    }
}
map< Vec<int64_t>, Vec<Mat<ZZ> >, vi64less>::const_iterator de_Rham_local::compute_reduction_matrix_J_ZZ(const Vec<int64_t> v)
{
    map< Vec<int64_t> ,  Vec<Mat<ZZ> >, vi64less>::const_iterator it_ZZ;
    it_ZZ = reduction_matrix_J_ZZ_dict.find(v);
    if(it_ZZ != reduction_matrix_J_ZZ_dict.end())
    {
        return it_ZZ;
    }
    else
    {
        map< Vec<int64_t> ,  Vec<Mat<ZZ_p> >, vi64less>::const_iterator it;
        it = reduction_matrix_J_dict.find(v);
        if( it == reduction_matrix_J_dict.end() )
        {
            it = compute_reduction_matrix_J(v);
            reduction_matrix_J_ZZ_dict[v] = conv< Vec< Mat<ZZ> > >(it->second);
            if( save_memory )
                //intel compiler dislikes
                //reduction_matrix_J_dict.erase( it );
                reduction_matrix_J_dict.erase( v );
        }
        else
            reduction_matrix_J_ZZ_dict[v] = conv< Vec< Mat<ZZ> > >(it->second);
        
        return reduction_matrix_J_ZZ_dict.find(v);
    }
}

Mat<ZZ_p>* de_Rham_local::get_inclusion_matrix_J(Vec<int64_t> u)
{
    int64_t i, sum;
    map< Vec<int64_t>, Mat<ZZ_p>, vi64less>::iterator it;

    sum = 0;
    for( i = 0; i <=n ; i++)
        sum += u[i];
    assert(sum == n);

    it = inclusion_matrix_J_dict.find(u);
    if( it == inclusion_matrix_J_dict.end() )
    {
        if( verbose )
            cout << "Computing the inclusion matrix J for u = "<<u<<endl;
        compute_inclusion_matrix_J(u);
        it = inclusion_matrix_J_dict.find(u);
        assert( it !=  inclusion_matrix_J_dict.end() ); 
    }
    return &(it->second);
}

/*
 * computes inclusion_matrix_J(u) is the matrix that maps
* V_u - > {(n-1)! H Omega/f^n : deg H = d*n - n - 1}
* and automatically adds it to the dict
*/
void de_Rham_local::compute_inclusion_matrix_J(Vec<int64_t> u)
{
    assert( inclusion_matrix_J_dict.find(u) == inclusion_matrix_J_dict.end() );
    Vec< Vec<int64_t> >* list_G;
    map< Vec<int64_t>, int64_t, vi64less >* dict_H;
    Mat<ZZ_p>* M;
    Vec<int64_t> w;
    int64_t i, coordinate_of_monomial;
    bool boolean;

    list_G = &tuple_list[ d * n - n ];
    dict_H = &tuple_dict[ d * n - n - 1 ];

    M = &(inclusion_matrix_J_dict[u]);
    M->SetDims( dict_H->size(), list_G->length());
    for(coordinate_of_monomial = 0; coordinate_of_monomial < (int64_t) list_G->length(); coordinate_of_monomial++)
    {
        w = (*list_G)[coordinate_of_monomial] + u;
        boolean = true;
        for( i = 0; i <= n ; i++)
        {
            if(w[i] > 0)
                w[i]--;
            else
            {
                boolean = false;
                break;
            }
        } 
        if(boolean)
            set( (*M)[  (*dict_H)[w] ][ coordinate_of_monomial ] );
    }
}

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
                    set((*M)[i][j]);
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



Mat<ZZ_p>* de_Rham_local::get_final_reduction_matrix_J(int64_t k)
{
    assert(k <= n);
    assert(k > 0);
    int64_t i;
    if( k > (int64_t)final_reduction_matrix_J_dict.length() )
    {
        i = final_reduction_matrix_J_dict.length() + 1;
        final_reduction_matrix_J_dict.SetLength(k);
        for(  ; i <= k ; i++)
        {
            if( verbose )
                cout <<"Computing the final reduction matrix J for pole = "<<i<<endl;
           compute_final_reduction_matrix_J(i);
        }
    }
    return &(final_reduction_matrix_J_dict[k-1]);
}

void de_Rham_local::matrix_J(Mat<ZZ_p> &result, int64_t l)
{
    //free storage and make 0 x 0
    result.kill();

    if( d - 1 > l)
    {
        result.SetDims( tuple_list[l].length(),0);
        return;
    }

    int64_t i, j;
    int64_t  len_list_lminusd1;
    

    Vec< Vec<int64_t> > *tuple_list_l;
    Vec< Vec<int64_t> > *tuple_list_lminusd1;

    map< Vec<int64_t>, int64_t, vi64less> *tuple_dict_l;
    map< Vec<int64_t>, ZZ_p, vi64less>::const_iterator itf;
   
    tuple_list_l = &tuple_list[l];
    tuple_list_lminusd1 = &tuple_list[l - (d - 1)];
    
    tuple_dict_l = &tuple_dict[l];
    
    len_list_lminusd1 = tuple_list_lminusd1->length();

    
    result.SetDims( tuple_list_l->length() , (n+1) * len_list_lminusd1);
    
    for(i = 0 ; i <= n ; i++)
    {
        for(j = 0; j < len_list_lminusd1; j++)
        {
            for(itf = f.begin(); itf != f.end(); itf++)
            {
                if( (itf->first)[i] > 0 )
                {
                    result[ (*tuple_dict_l)[ diff(itf->first, i) + (*tuple_list_lminusd1)[j] ]][ i*len_list_lminusd1 + j ] = (itf->first)[i] * itf->second;
                }
            }
         }
    }
}

Vec<ZZ_p> de_Rham_local::monomial_to_basis_J(const Vec<int64_t> u)
{
    assert( (int64_t) u.length() == n+1 );
    int64_t i, m, sum;

    sum = 0;
    for( i = 0; i <= n ; i++)
    {
        sum += u[i];
        assert( u[i] > 0 );
    }
    assert( sum%d == 0);

    m = sum / d;
    
    if( m > n)
    {
        int64_t e;
        Vec<int64_t> s, v, w;
        Mat<ZZ_p> M_red;
        Mat<ZZ_p>* M;
        Vec<ZZ_p> G, G_new;
        G.SetLength( tuple_list[d*n-n].length() );
        s.SetLength(n+1);
        v.SetLength(n+1);
        w = tweak( u, d * n - n );
        sub(s, u, w);
        // u = w + s
        set( G[ tuple_dict[d*n - n][s] ] );

        for( e = m; e > n; e--)
        {
            sum = 0;
            for( i = 0; i <= n; i++)
            {
                if(w[i] > 0)
                {
                    v[i] = 1;
                    sum++;
                }
                else
                {
                    v[i] = 0;
                }
            }
            while( sum < d )
            {
                for( i = 0; i <= n; i++)
                {
                    if( w[i] > 0 && v[i] < w[i] )
                    {
                        v[i]++;
                        sum++;
                        break;
                    }
                }
            }
            //w = w - v;
            sub(w,w,v);
            M_red = get_reduction_matrix_J(w, v);
            mul(G, M_red, G);
            //mul(G_new, M_red, G);
            //swap(G_new, G);
            
        }
        sum = 0;
        for(i = 0; i <= n; i++)
            sum += w[i];

        assert(sum == n);
        // Now we have:
        // (m-1)! x^u \Omega / x0...xn ~ (n-1)! x^w \Omega / x0 ... xn f^n
        M = get_inclusion_matrix_J(w);

        //mul(G_new, *M, G);
        mul(G, *M, G);


        M = get_final_reduction_matrix_J(n);
        //mul(G, *M, G_new);
        mul(G, *M, G);

        return G;
    }
    else
    {
        Vec<ZZ_p> G, G_new;
        Vec<int64_t> w;
        Mat<ZZ_p>* M;

        G.SetLength( tuple_list[ sum - (n+1)].length() );
        w = u;
        for( i = 0; i <= n ; i++)
            w[i]--;

        set( G[ tuple_dict[ sum - (n+1)][w] ] );
        M = get_final_reduction_matrix_J(m);
        mul(G_new, *M, G);

        return G_new;
    }
}

bool de_Rham_local::test_monomial_to_basis_J(int64_t N)
{
    if(verbose)
        cout <<"de_Rham_local::test_monomial_to_basis_J("<<N<<")"<<endl;

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

    set( fpow[0][u] );

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
            cout<<"i = "<<i<<endl;
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
                v += it1->second * monomial_to_basis_J(u);
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

bool de_Rham_local::test_paths_J(int64_t trials, int64_t paths)
{
    if(verbose)
        cout << "de_Rham_local::test_paths_J("<<trials<<", "<<paths<<")\n";
    assert(trials > 0);
    assert(paths > 0);

    int64_t attempt, pathi, i, k, sum, sum_u;
    Vec<int64_t> u_orig, u, v;
    Vec< Mat<ZZ_p> > M_saved;
    Mat<ZZ_p> M;
    Mat<ZZ_p>* Mfinalred;
    Mat<ZZ_p>* Minclusion;
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
            u_orig[i] = RandomBnd(20);
            sum += u_orig[i];
        }
        
        while( sum % d != n )
        {
            k = RandomBnd(n+1);
            u_orig[k]++;
            sum++;
        }


        for( pathi = 0 ; pathi < paths ; pathi++)
        {
            if(verbose)
                cout << "path = "<<pathi<<endl;
            M.kill();
            M.SetDims( tuple_list[d*n-n].length(), tuple_list[d*n-n].length());
            for( i = 0; i < (int64_t) tuple_list[d*n-n].length(); i++)
                set( M[i][i] );
            sum_u = 0;
            u = u_orig;
            for( i = 0; i < n+1; i++)
                sum_u += u[i];
            while( sum_u > n)
            {
                sum = 0;
                for( i = 0; i <= n; i++)
                {
                    if( u[i] > 0 )
                    {
                        v[i] = 1;
                        sum++;
                    }
                    else
                        v[i] = 0;
                }
                while( sum < d)
                {
                    k = RandomBnd(n + 1);
                    if( u[k] > 0 and v[k]<u[k])
                    {
                        v[k]++;
                        sum++;
                    }
                }
                u = u - v;
                sum_u -= d;

                Mred = get_reduction_matrix_J(u,v);
                mul(M, Mred, M);
            }
            Minclusion = get_inclusion_matrix_J(u);
            mul(M, *Minclusion, M);

            Mfinalred = get_final_reduction_matrix_J(n);
            mul(M, *Mfinalred, M);
            M_saved[pathi] = M;
            if(pathi > 0)
                if (M_saved[pathi-1] != M_saved[pathi])
                    return false;
        }
    }
    return true;
}




bool isSmooth(const map< Vec<int64_t>, zz_p, vi64less> &f )
{
    int64_t i, j;
    de_Rham_local D;
    D = de_Rham_local();
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
    D.matrix_J(M, (D.d-2)*(D.n+1) + 1 );
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
    flint_free(P);

    return result;
}
