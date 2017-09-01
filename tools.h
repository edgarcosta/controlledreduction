// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// tools.h: header file for miscellaneous routines

#ifndef TOOLS_H_
#define TOOLS_H_

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

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
        std::cerr << " @ " << __FILE__ << ":" << __LINE__  << std::endl; \
        std::cerr << #left << " = " << (left) << "; "; \
        std::cerr << #right << " = " << (right) << std::endl; \
        abort(); \
    } \
}
#else
#define assert_print(condition, statement) ((void)0)
#endif


#define print(var) { std::cout << #var << " = " << (var) << std::endl;}

/*
template <class Key, class T, class Compare = std::less<Key>,
              class Allocator = std::allocator<std::pair<const Key, T> > >
using std::map = typename std::std::map<Key, T, Compare, Allocator>;
*/

// returns binomial(n,k)
int64_t binomial(const int64_t n, const int64_t k);
NTL::ZZ binomial_ZZ(const int64_t n, const int64_t k);


// change of variables
template<class T>
std::map< NTL::Vec<int64_t>, T, vi64less> 
change_of_variables( std::map< NTL::Vec<int64_t>, T, vi64less> f, NTL::Mat<T> M) {
    int64_t n = M.NumRows();
    assert(n == M.NumCols() );

    std::map< NTL::Vec<int64_t>, T, vi64less> result;
    NTL::Vec< std::map< NTL::Vec<int64_t>, T, vi64less> > Mvec;
    typename std::map< NTL::Vec<int64_t>, T, vi64less >::iterator itf, itH, itHnew, itresult, itX;

    NTL::Vec<int64_t> zero, monomial;
    NTL::Vec< NTL::Vec<int64_t> > X;

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

    // X[i] std::mapped to Mvec[i]
    for (int64_t i = 0; i < n; ++i)
        for (int64_t j = 0; j < n; ++j)
            Mvec[i][X[j]] = M[i][j];

    for (itf = f.begin(); itf != f.end(); ++itf) {
        assert((int64_t) itf->first.length() == n);
        std::map< NTL::Vec<int64_t>, T, vi64less> H, Hnew;

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
std::map< NTL::Vec<int64_t>, int64_t, vi64less>
change_of_variables_monomial(const NTL::Vec<int64_t> u, const NTL::Vec<int64_t> v);


// a look up table for the triple (precision, N, charpoly_prec)
// for n =  2 or 3 and n < d <= 5
void default_args(int64_t &precision, NTL::Vec<int64_t> &N, NTL::Vec<int64_t> &charpoly_prec, const int64_t &p, const int64_t &n, const int64_t &d);


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
int64_t factorial_p_adic(NTL::ZZ_p &result, const int64_t n, const int64_t p, const int64_t start = 1);

// given the characteristic polynomial of Frob acting on the middle cohomology, the field characteristic and the dimension computes the geometric_picard
int64_t geometric_picard(const NTL::Vec<NTL::ZZ> &charpoly_frob, const int64_t prime, const int64_t dimension);

// generates all the tuples of sum = d and with n variables
void tuple_list_generator(NTL::Vec< NTL::Vec<int64_t> > &result, const int64_t d, const int64_t n);



// returns v_p(n!)
int64_t valuation_of_factorial(const int64_t n, const int64_t p);


/*
 * extending << and >>
 */
// istream for a std::map< NTL::Vec<T>, R, Compare>
template<typename T, typename R, typename Compare>
NTL_SNS istream & operator>> (NTL_SNS istream& s, std::map< NTL::Vec<T>, R, Compare>& a) {
    NTL::Vec< NTL::Vec<T> > monomials;
    NTL::Vec<R> coefficients;

    s >> monomials;
    s >> coefficients;
    assert_print(monomials.length(), ==, coefficients.length());
    std::map< NTL::Vec<T>, R, Compare> ibuf;
    for (int64_t i = 0; i < coefficients.length(); ++i)
        ibuf[monomials[i]] = coefficients[i];
    a = ibuf;
    return s;
}

// ostream for a std::map< NTL::Vec<T>, R, Compare>
template<typename T, typename R, typename Compare>
NTL_SNS ostream & operator<< (NTL_SNS ostream& s, const  std::map< NTL::Vec<T>, R, Compare>& a) {
    typename std::map< NTL::Vec<T>, R, Compare>::const_iterator it;
    int64_t i;

    NTL::Mat<T>  monomials;
    NTL::Vec<R> coefficients;
        
    monomials.SetDims(a.size(), a.cbegin()->first.length());
    coefficients.SetLength(a.size());
    for (i = 0, it = a.cbegin(); it != a.cend(); ++it, ++i){
        monomials[i] = it->first;
        coefficients[i] = it->second;
    }
    s << monomials<<std::endl;
    s << coefficients << std::endl;
    return s;
}


template<typename S>
NTL_SNS ostream & operator<<= (NTL_SNS ostream& s, const S& a){
    s << a;
    return s;
}

// similar to << but more human readable and compatible with SAGE 
template<typename T>
NTL_SNS ostream & operator<<= (NTL_SNS ostream& s, const NTL::Mat<T>& a) {
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
NTL_SNS ostream & operator<<=(NTL_SNS ostream& s, const NTL::Vec<T>& a) {
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
NTL_SNS ostream & operator<<= (NTL_SNS ostream& s, const std::map< NTL::Vec<T>, R, Compare>& a) {
    s << "{";
    int64_t i, n;
    i = 0;
    n = a.size();
    typename  std::map< NTL::Vec<T>, R, Compare>::const_iterator it;
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
