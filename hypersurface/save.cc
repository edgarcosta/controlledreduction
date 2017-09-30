// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// 
// hypersurface.cpp: routines for the classes hypersurface
  
#include "hypersurface.h"

bool hypersurface::save(const char * filename)
{
    return (dR->save)(filename);
}

