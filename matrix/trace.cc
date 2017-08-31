// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "matrix.h"

#include <assert.h>
#include <stdint.h>

#include <NTL/matrix.h>
#include <NTL/ZZ.h>

using namespace NTL;

ZZ trace(const Mat<ZZ> M) {
    ZZ result(0);
    assert( M.NumRows() ==  M.NumCols() );

    for(int64_t i = 0; i < M.NumRows(); ++i)
        result += M[i][i];

    return result;
}

