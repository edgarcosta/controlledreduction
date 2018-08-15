// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
//
// test_dr.cpp: tests if frob_J == frob_ND, and also checks that the charpoly is the same if we change variables

#include "hypersurface_nd.h"
#include "time.h"
#include "timing.h"
#include "tools.h"

#define NTL_RANGE_CHECK
void test_JND_and_invariance(int64_t p, int64_t n)
{
  int64_t precision, d;

  if(n != 2 && n!=3)
  {
    cout << "n most be 2 or 3."<<endl;
    abort();
  }
  d = 4;
  Vec<int64_t> N;
  N.SetLength(n);
  if(n == 2)
  {
    N[0] = 2;
    N[1] = 2;
  }
  if(n == 3)
  {
    N[0] = 3;
    N[1] = 3;
    N[2] = 2;
  }
  if( n == 2)
  {
    precision = 3;
  }
  if( n == 3)
  {
    precision = 4;
  }
  {
    zz_pPush push(p);
    ZZ_pPush push2(power_ZZ(p,precision));

    hypersurface_non_degenerate hs;
    hs = hypersurface_non_degenerate(p, precision, n, d, false);
    hs.verbose = false;
    Mat<ZZ_p> Fp_J = hs.frob_matrix_J(N);
    Mat<ZZ_p> Fp_ND =  hs.frob_matrix_ND(N);

    if(Fp_J != Fp_ND)
    {
      cout <<"\n\tFAIL\n";
      cout <<"Fp_J != Fp_ND"<<endl;
      cout <<= hs.dR->fbar ;
      cout << endl;
      cout <<"p = " << p << " n = "<<n<<" d = "<<d<<" N = "<< N<<endl;
      cout << Fp_J <<endl<<Fp_J<<endl;

    }

    Mat<ZZ> Fzz;
    conv(Fzz, Fp_J);
    Vec<ZZ> cp = charpoly(Fzz);

    map< Vec<int64_t>, zz_p , vi64less> new_fbar, random_fbar;

    random_fbar =  change_of_variables((hs.dR)->fbar, random_SL_matrix(n+1));

    Mat<zz_p> M = find_change_of_variables( random_fbar, 1000);
    new_fbar = change_of_variables( random_fbar, M);

    hypersurface_non_degenerate hs2;

    hs2 = hypersurface_non_degenerate(p, precision, new_fbar, false);
    hs2.verbose = false;
    //cout << "frob_matrix_ND(N)" <<endl;
    Mat<ZZ_p> Fp_new = hs2.frob_matrix_ND(N);
    Mat<ZZ> M1, M2;
    M1.SetDims(Fp_ND.NumRows(),  Fp_ND.NumCols());
    M2.SetDims(Fp_ND.NumRows(),  Fp_ND.NumCols());
    ZZ px, p2, p3;
    p2 = power_ZZ(p,2);
    p3 = power_ZZ(p,3);
    int64_t i, j;
    for( i = 0; i < (int64_t) Fp_ND.NumRows(); i++)
    {
      for(j = 0 ; j < (int64_t) Fp_ND.NumCols(); j++)
      {
        if(n == 2)
        {
          rem(M1[i][j], rep(Fp_ND[i][j]), p2);
          rem(M2[i][j], rep(Fp_new[i][j]), p2);
        }
        if(n == 3)
        {
          if(j != 20)
          {
            rem(M1[i][j], rep(Fp_ND[i][j]), p3);
            rem(M2[i][j], rep(Fp_new[i][j]), p3);
          }
          else
          {
            rem(M1[i][j], rep(Fp_ND[i][j]), p2);
            rem(M2[i][j], rep(Fp_new[i][j]), p2);
          }
        }
      }
    }
    Vec<ZZ> c1, c2;
    c1 = charpoly(M1);
    c2 = charpoly(M2);
    /*quartic curve*/
    if( n == 2)
    {
      for(i = 0; i < 3; i++)
      {
        px = power_ZZ(p,4-i);
        c1[i] = c1[i] % px;
        c2[i] = c2[i] % px;
      }
      for(;i < 7; i++)
      {
        rem(c1[i], c1[i], p2);
        rem(c2[i], c2[i], p2);
      }
    }
    if( n == 3)
    {
      for(i = 0; i < 21; i++)
      {
        px = power_ZZ(p,22-i);
        c1[i] = c1[i] % px;
        c2[i] = c2[i] % px;
      }
    }

    if (c1 != c2){
      cout <<"\n\tFAIL\n";
      cout <<"cp(Fp_ND(change of variables) != cp(Fp_ND)"<<endl;
      cout <<= hs.dR->fbar ;
      cout << endl;
      cout <<= new_fbar ;
      cout << endl;
      cout <<"change of variables = "<<endl;
      cout <<M<<endl;
      cout <<"p = " << p << " n = "<<n<<" d = "<<d<<" N = "<< N<<endl;
      cout << M1 <<endl<<M2<<endl;
      cout << c1 << endl << c2 <<endl;
      abort();
    }

  }

}


int main()
{
  timestamp_pair pair;
  timestamp_mark(pair);
  SetSeed(to_ZZ(time(NULL)));
  cout<<"test_JND_and_invariance...";
  fflush(stdout);
  int64_t n = 2;
  int64_t runs = 20;
  int64_t plen = 5;
  for (int64_t i = 0; i < runs; ++i) {
    int64_t p = GenPrime_long(plen);
    while ( p <= 4 )
      p = GenPrime_long(plen);
    // cout << p << " " << n << " " << d << endl;
    test_JND_and_invariance(p, n);
  }
  /*
  // it takes too long to be a regular test
  n = 3;
  int64_t d = 4;
  int64_t p = GenPrime_long(5);
  while ( p <= 23 )
  p = GenPrime_long(5);
  cout << p << " " << n << " " << d << endl;
  test_invariance_ND(p, n, d);
  */

  cout << "PASS ";
  timestamp_report(pair);
  cout << endl;
  return 0;
}

