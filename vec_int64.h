// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#ifndef CONTROLLED_REDUCTION_VEC_INT64_H
#define CONTROLLED_REDUCTION_VEC_INT64_H

#include <cstdint>
#include <assert.h>


#include <NTL/vector.h>
#include <NTL/matrix.h>

using namespace std;
using namespace NTL;



//extending Vec<int64_t>

bool operator<(const Vec<int64_t>& a, const Vec<int64_t>& b);

typedef struct vec_int64_t_less
{
    bool operator()(const Vec<int64_t>& a, const Vec<int64_t>& b) const
    {
        return a < b;
    }
} vi64less;

inline bool operator<(const Vec<int64_t>& a, const Vec<int64_t>& b)
{
    if( a.length() != b.length())
        Error("operator<: dimension mismatch");

    for(int64_t i = 0; i < a.length() - 1; ++i)
    {
        if( a[i] > b[i] )
        {
            return false;
        }
        else if( a[i] < b[i] )
        {
            return true;
        }
    }
    return a[a.length() - 1] < b[b.length() - 1];
}





// res = c * v;
void inline mul(Vec<int64_t> &res, const int64_t c, const Vec<int64_t> &v)
{
    res.SetLength(v.length());
    for(int64_t i = 0; i < v.length(); ++i)
        res[i] = c * v[i];
}

// res = v + w
void inline add(Vec<int64_t> &res, const Vec<int64_t> &v, const Vec<int64_t> &w)
{
    assert(v.length() == w.length());
    res.SetLength(v.length());
    for(int64_t i = 0; i < v.length(); ++i)
        res[i] = v[i] + w[i];
}
// res = v - w
void inline sub(Vec<int64_t> &res, const Vec<int64_t> &v, const Vec<int64_t> &w)
{
    assert(v.length() == w.length());
    res.SetLength(v.length());
    for(int64_t i = 0; i < v.length(); ++i)
        res[i] = v[i] - w[i];
}

// res = v*w
void inline mul(int64_t &res, const Vec<int64_t> &v, const Vec<int64_t> &w)
{
    assert(v.length() == w.length());
    res = 0;
    for(int64_t i = 0; i < v.length(); ++i)
        res += v[i]*w[i];
}


// res = A*v
void inline mul(Vec<int64_t> &res, const Mat<int64_t> &A, const Vec<int64_t> &v)
{
    assert(A.NumCols() == v.length());
    res.SetLength(A.NumRows());
    for(int64_t i =0; i < A.NumRows(); ++i)
        mul(res[i], A[i], v);
}

/*
 * Operators
 */

inline Vec<int64_t> operator+(const Vec<int64_t>& a, const Vec<int64_t>& b)
{
    Vec<int64_t> x;
    add(x,a,b);
    return x;
}

inline Vec<int64_t> operator-(const Vec<int64_t>& a, const Vec<int64_t>& b)
{
    Vec<int64_t> x;
    sub(x,a,b);
    return x;
}

inline Vec<int64_t> operator*(const int64_t c, const Vec<int64_t>& v)
{
    Vec<int64_t> x;
    mul(x, c , v);
    return x;
}

inline int64_t operator*(const Vec<int64_t> &v, const Vec<int64_t> &w)
{
    int64_t x;
    mul(x, v , w);
    return x;
}

inline Vec<int64_t> operator*(const Mat<int64_t> &A, const Vec<int64_t> &w)
{
    Vec<int64_t> x;
    mul(x, A , w);
    return x;
}


//converts Vec to array
inline void conv(int64_t* result, const Vec<int64_t> &v)
{
    for(int64_t i = 0; i < v.length(); ++i)
        result[i] = v[i];
}

//converts an array to Vec<int64_t>

inline void conv(Vec<int64_t> &result, const int64_t array[], const int64_t size)
{
    result.SetLength(size);
    for(int64_t i = 0; i < size ; ++i)
    {
        result[i] = array[i];
    }
}

/*
 * returns v-ei
 */
Vec<int64_t> diff(const Vec<int64_t> v, const int64_t i);

// reduces an entry by 1
Vec<int64_t> tweak_step( const Vec<int64_t> v);

// applies tweak r times
Vec<int64_t> tweak( const Vec<int64_t> v, const int64_t r);


#endif


