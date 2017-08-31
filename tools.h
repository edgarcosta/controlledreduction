// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// tools.h: header file for miscellaneous routines

#ifndef TOOLS_H_
#define TOOLS_H_

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#pragma once
#ifdef _OPENMP
# include <omp.h>
#endif

#include <flint/arith.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>



#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vector.h>
#include <NTL/ZZ.h>

#include <cstddef>
#include <map>

#include "vec_int64.h"


#ifndef NDEBUG
#define assert_print(left, operator, right) \
{ \
    if ( !( (left) operator (right) ) ){ \
        std::cerr << "ASSERT FAILED: "; \
        std::cerr << #left << " " << #operator << " " << #right; \
        std::cerr << " @ " << __FILE__ << ":" << __LINE__  << endl; \
        std::cerr << #left << " = " << (left) << "; "; \
        std::cerr << #right << " = " << (right) << endl; \
        abort(); \
    } \
}
#else
#define assert_print(condition, statement) ((void)0)
#endif


#define print(var) { cout << #var << " = " << (var) << endl;}




using namespace NTL;



// returns binomial(n,k)
int64_t binomial(const int64_t n, const int64_t k);
ZZ binomial_ZZ(const int64_t n, const int64_t k);


// change of variables
template<class T>
map< Vec<int64_t>, T, vi64less>
change_of_variables( map< Vec<int64_t>, T, vi64less> f, Mat<T> M) {
    int64_t n = M.NumRows();
    assert(n == M.NumCols() );

    map< Vec<int64_t>, T, vi64less> result;
    Vec< map< Vec<int64_t>, T, vi64less> > Mvec;
    typename map< Vec<int64_t>, T, vi64less >::iterator itf, itH, itHnew, itresult, itX;

    Vec<int64_t> zero, monomial;
    Vec< Vec<int64_t> > X;

    zero.SetLength(n);
    X.SetLength(n);
    for (int64_t i = 0; i < n; ++i)
        zero[i] = 0;

    monomial = zero;
    // X[i] = ei
    for (int64_t i = 0; i < n; ++i) {
        ++monomial[i];
        X[i] = monomial;
        --monomial[i];
    }

    Mvec.SetLength(n);

    // X[i] mapped to Mvec[i]
    for (int64_t i = 0; i < n; ++i)
        for (int64_t j = 0; j < n; ++j)
            Mvec[i][X[j]] = M[i][j];

    for (itf = f.begin(); itf != f.end(); ++itf) {
        assert((int64_t) itf->first.length() == n);
        map< Vec<int64_t>, T, vi64less> H, Hnew;

        H[zero] = itf->second;

        for (int64_t i = 0; i < n ; ++i) {
            for (int64_t j = 0; j < itf->first[i]; ++j) {
                Hnew.clear();
                for (itH = H.begin(); itH != H.end(); ++itH) {
                    for (itX = Mvec[i].begin(); itX != Mvec[i].end(); ++itX) {
                        monomial = itX->first + itH->first;
                        itHnew = Hnew.find(monomial);
                        if ( itHnew == Hnew.end() )
                            Hnew[ monomial ] = itX->second * itH->second;
                        else
                            itHnew->second += itX->second * itH->second;
                    }
                }
                H.swap(Hnew);
            }
        }

        for (itH = H.begin(); itH != H.end(); ++itH) {
            itresult = result.find(itH->first);
            if (itresult == result.end() )
                result[itH->first] = itH->second;
            else
                itresult->second += itH->second;
        }
    }
    return result;
}

// how to write (x0 + x v0)^u0 ... (xn + x vn)^un in monomials in x0...xn and x
map< Vec<int64_t>, int64_t, vi64less>
change_of_variables_monomial(const Vec<int64_t> u, const Vec<int64_t> v);


//conversion between ZZ, mpz and fmpz
//b -> a

void conv(mpz_t a, const ZZ &b);
void conv(fmpz_t a, const ZZ &b);
void conv(ZZ &a, const mpz_t b);
void conv(ZZ &a, const fmpz_t b);

//conversion between Vec/Mat<T> and it's equivalents over fmpz_t, mpz_t and "nmod"
void conv(mpz_t* A, const Mat<ZZ> &B);

template<class T>
void conv(fmpz_mat_t A, const Mat<T> &B)
{
    int64_t rowsB = B.NumRows();
    int64_t colsB = B.NumCols();
    assert( rowsB = (int64_t) fmpz_mat_nrows(A));
    assert( colsB = (int64_t) fmpz_mat_ncols(A));
    int64_t i, j;
    for(i = 0; i < rowsB; ++i)
    {
        for(j = 0; j < colsB; j++)
        {
            conv(fmpz_mat_entry(A, i, j), conv<ZZ>(B[i][j]));
        }
    }
}

template<class T>
void conv(fmpz * a, const Vec<T> &b)
{
    int64_t len = b.length();
    int64_t i;
    ZZ tmp;
    for(i = 0; i < len; ++i)
        conv(a + i,conv<ZZ>(b[i]));
}

template<class T>
void conv(nmod_mat_t A, const Mat<T> &B)
{
    int64_t rowsB = B.NumRows();
    int64_t colsB = B.NumCols();
    assert( rowsB = (int64_t) nmod_mat_nrows(A));
    assert( colsB = (int64_t) nmod_mat_ncols(A));
    int64_t i,j;
    for(i = 0; i < rowsB; ++i)
    {
        for(j = 0; j < colsB; j++)
        {
            nmod_mat_entry(A, i, j) = conv<ulong>(rep(B[i][j])% A->mod.n );
        }
    }
}

void conv(Mat<ZZ> &A, mpz_t* B, int64_t rowsB, int64_t colsB);

template<class T>
void conv(Mat<T> &A, const fmpz_mat_t B)
{
    int64_t rowsA =  fmpz_mat_nrows(B);
    int64_t colsA =  fmpz_mat_ncols(B);
    A.SetDims(rowsA, colsA);
    int64_t i,j;
    ZZ tmp;
    for(i = 0; i < rowsA; ++i)
    {
        for(j = 0; j < colsA; j++)
        {
            clear(tmp);
            conv(tmp, fmpz_mat_entry(B, i, j));
            conv(A[i][j], tmp);
        }
    }
}

template<class T>
void conv(Vec<T> &a, const fmpz * b, int64_t len)
{
    int64_t i;
    a.SetLength(len);
    ZZ tmp;
    for(i = 0; i < len; ++i)
    {
        conv(tmp, b + i);
        conv(a[i],tmp);
    }
}


template<class T>
void conv(Mat<T> &A, const nmod_mat_t B)
{
    int64_t rowsA = nmod_mat_nrows(B);
    int64_t colsA = nmod_mat_ncols(B);
    A.SetDims(rowsA, colsA);
    int64_t i,j;
    for(i = 0; i < rowsA; ++i)
    {
        for(j = 0; j < colsA; j++)
        {
            conv(A[i][j], nmod_mat_entry(B, i, j));
        }
    }
}


// returns factorial(n)/factorial(start-1)
template<class T>
T  factorial (int64_t n, int64_t start = 1) {
    T result;
    result = 1;
    for (int64_t i = start; i <=n ; ++i)
        result *= i;
    return result;
}

// returns v_p(n!) and result = (n! / p^{v_p(n!)})
int64_t factorial_p_adic(ZZ_p &result, const int64_t n, const int64_t p, const int64_t start = 1);

// given the characteristic polynomial of Frob acting on the middle cohomology, the field characteristic and the dimension computes the geometric_picard
int64_t geometric_picard(const Vec<ZZ> &charpoly_frob, const int64_t prime, const int64_t dimension);

// generates all the tuples of sum = d and with n variables
void tuple_list_generator(Vec< Vec<int64_t> > &result, const int64_t d, const int64_t n);



// returns v_p(n!)
int64_t valuation_of_factorial(const int64_t n, const int64_t p);


/*
 * extending << and >>
 */
// istream for a map< Vec<T>, R, Compare>
template<typename T, typename R, typename Compare>
NTL_SNS istream & operator>> (NTL_SNS istream& s, map< Vec<T>, R, Compare>& a) {
    Vec< Vec<T> > monomials;
    Vec<R> coefficients;

    s >> monomials;
    s >> coefficients;
    assert_print(monomials.length(), ==, coefficients.length());
    map< Vec<T>, R, Compare> ibuf;
    for (int64_t i = 0; i < coefficients.length(); ++i)
        ibuf[monomials[i]] = coefficients[i];
    a = ibuf;
    return s;
}

// ostream for a map< Vec<T>, R, Compare>
template<typename T, typename R, typename Compare>
NTL_SNS ostream & operator<< (NTL_SNS ostream& s, const  map< Vec<T>, R, Compare>& a) {
    typename map< Vec<T>, R, Compare>::const_iterator it;
    int64_t i;

    Mat<T>  monomials;
    Vec<R> coefficients;
        
    monomials.SetDims(a.size(), a.cbegin()->first.length());
    coefficients.SetLength(a.size());
    for (i = 0, it = a.cbegin(); it != a.cend(); ++it, ++i){
        monomials[i] = it->first;
        coefficients[i] = it->second;
    }
    s << monomials<<endl;
    s << coefficients << endl;
    return s;
}


template<typename S>
NTL_SNS ostream & operator<<= (NTL_SNS ostream& s, const S& a){
    s << a;
    return s;
}

// similar to << but more human readable and compatible with SAGE 
template<typename T>
NTL_SNS ostream & operator<<= (NTL_SNS ostream& s, const Mat<T>& a) {
    s <<"[";
    for (int64_t i = 0; i <  a.NumRows(); ++i) {
        s <<= a[i];
        if (i <  a.NumRows() - 1)
            s << ",\n";
    }
    s << "]";
    return s;
}

// similar to << but more human readable and compatible with SAGE 
template<typename T>
NTL_SNS ostream & operator<<=(NTL_SNS ostream& s, const Vec<T>& a) {
    s << "(";
    for (int64_t i = 0; i < a.length(); ++i) {
        s <<= a[i];
        if ( i < a.length() - 1 )
            s << ", ";
    }
    s << ")";
    return s;
}

// similar to << but more human readable and compatible with SAGE
template<typename T, typename R, typename Compare>
NTL_SNS ostream & operator<<= (NTL_SNS ostream& s, const map< Vec<T>, R, Compare>& a) {
    s << "{";
    int64_t i, n;
    i = 0;
    n = a.size();
    typename  map< Vec<T>, R, Compare>::const_iterator it;
    for (it = a.cbegin() ; a.cend() != it ; ++it) {
        s <<= it->first;
        s << ": ";
        s << it->second;
        if (i < n - 1)
            s <<",\n ";
        /*else
            break;*/
        ++i;
    }
    s << "}";
    return s;
}




#endif  // TOOLS_H_
