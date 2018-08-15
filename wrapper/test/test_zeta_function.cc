// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "wrapper.h"
#include <vector>
#include <NTL/ZZX.h>
#include <sstream>
#define NTL_RANGE_CHECK

using namespace std;
using namespace NTL;


int test_dwork() {
  vector<ZZ> coef;
  vector< vector<int64_t> > mon;
  int64_t d = 4;
  for(int64_t i = 0; i < 4; ++i) {
    vector<int64_t> v = {0, 0, 0, 0};
    v[i] = d;
    mon.push_back(v);
    coef.push_back(ZZ(1));
  }
  vector<int64_t> v = {1, 1, 1, 1};
  mon.push_back(v);
  coef.push_back(ZZ(1));
  ZZX cp;

  cout<<"zeta_function(ZZX, vector< vector<int64_t> >, vector<ZZ>, 23, false, 1) K3 dwork...";
  fflush(stdout);
  zeta_function(cp, mon, coef, 23, false, 1);
  ZZX expected;

  stringstream buffer;
  buffer << "[-7400249944258160101211 3608386336456458231169 -277995865674611574050 -136470697694809318170 25042602775646827745 1533050078175543917 -641017719088510808 9666359400300680 8818983738586010 -482176723332590 -74181034358860 6743730396260 362266508890 -54758950510 -496037080 271854088 -5373247 -725395 32670 550 -59 1]";
  buffer << expected;
  if(cp == expected) {
    cout <<"FAIL"<<endl;
    cout <<"got = "<<cp<<endl;
    cout <<"expected = "<<expected<<endl;
    return 1;
  }
  cout << "PASS "<<endl;
  return 0;
}

int test_genus3() {
  char input[] = "17\n[[1 2 1] [2 0 2] [3 1 0] [1 3 0] [2 2 0] [0 2 2] [1 1 2] [3 0 1] [0 0 4] [2 1 1] [1 0 3] [0 4 0]]\n[-3 1 1 -4 1 3 -3 1 2 3 -4 2]";
  ZZX out;
  cout<<"test zeta_function(ZZX, string, false, 1) g=3...";
  fflush(stdout);
  zeta_function(out, input, false, 1);
  stringstream buffer;
  buffer << "[4913 -4046 1887 -548 111 -14 1]";
  ZZX expected;
  buffer << expected;
  if(out == expected) {
    cout <<"FAIL"<<endl;
    cout <<"got = "<<out<<endl;
    cout <<"expected = "<<expected<<endl;
    return 1;
  }
  cout << "PASS "<<endl;
  return 0;
}


int test_K3() {
  char input[] = "23\n[[4 0 0 0] [1  3  0  0]  [0  0  1  3]  [3  1  0  0]  [0  0  3  1]]\n[1 1  1  1  1]";
  ZZX out;
  cout<<"test zeta_function(ZZX, string, false, 1) K3...";
  fflush(stdout);
  zeta_function(out, input, false, 1);
  stringstream buffer;
  buffer << "[-39471584120695485887249589623 1716155831334586342923895201 694248294717583133054278966 -30184708465981875350186042 -5562277647022665774607627 241838158566202859765549 26802637582333589558408 -1165332068797112589496 -86212132515567942814 3748353587633388818 193710366435022724 -8422189845000988 -308075416095454 13394583308498 342258975368 -14880825016 -253815787 11035469 113206 -4922 -23 1]";
  ZZX expected;
  buffer << expected;
  if(out == expected) {
    cout <<"FAIL"<<endl;
    cout <<"got = "<<out<<endl;
    cout <<"expected = "<<expected<<endl;
    return 1;
  }
  cout << "PASS"<<endl;
  return 0;
}

int main() {
  int total = 0;
  total += test_genus3();
  total += test_dwork();
  total += test_K3();
  if( total > 0 )
    return 1;
  else
    return 0;
}





