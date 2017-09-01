// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
//
// checks that the charpoly is the same if we change variables

#include <sstream>

#include<NTL/vector.h>

//we can do this much better..... see toric code
void default_args(int64_t &precision, NTL::Vec<int64_t> &N, NTL::Vec<int64_t> &charpoly_prec, const int64_t &p, const int64_t &n, const int64_t &d)
{
    std::stringstream buffer;
    if ( n == 2 && d == 3 ) {
		if ( p < 5) {
			precision = 5;
			buffer << "[2 3]\n[2 2 2]";
		} else if ( 5 <= p && p < 17 ) {
			precision = 3;
			buffer << "[2 2]\n[2 2 2]";
		} else if ( 17 <= p  ) {
			precision = 1;
			buffer << "[0 1]\n[3 1 1]";
        }
	} else if ( n == 2 && d == 4 ) {
		if ( p < 5) {
			precision = 7;
			buffer << "[4 4]\n[5 4 3 3 3 3 3]";
		} else if ( 5 <= p && p < 17 ) {
			precision = 5;
			buffer << "[3 3]\n[5 4 3 3 3 3 3]";
		} else if ( 17 <= p  ) {
			precision = 3;
			buffer << "[2 2]\n[4 3 2 2 2 2 2]";
        }
	} else if ( n == 2 && d == 5 ) {
		if ( p < 5) {
			precision = 12;
			buffer << "[6 6]\n[10 9 8 7 6 5 5 5 5 5 5 5 5]";
		} else if ( 5 <= p && p < 7 ) {
			precision = 9;
			buffer << "[4 5]\n[9 8 7 6 5 4 4 4 4 4 4 4 4]";
		} else if ( 7 <= p  ) {
			precision = 7;
			buffer << "[4 4]\n[9 8 7 6 5 4 4 4 4 4 4 4 4]";
        }
	} else if ( n == 3 && d == 4 ) {
		if ( p < 5) {
			precision = 16;
			buffer << "[7 7 8]\n[24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 5 5]";
		} else if ( 5 <= p && p < 7 ) {
			precision = 9;
			buffer << "[4 5 5]\n[23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 4 4]";
		} else if ( 7 <= p && p < 23 ) {
			precision = 6;
			buffer << "[4 4 3]\n[23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 3]";
		} else if ( 23 <= p && p < 43 ) {
			precision = 5;
			buffer << "[3 3 3]\n[22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 3 3]";
		} else if ( 43 <= p  ) {
			precision = 4;
			buffer << "[3 3 2]\n[22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 2]";
        }
	} else if ( n == 3 && d == 5 ) {
		if ( p < 5) {
			precision = 21;
			buffer << "[11 11 10]\n[58 56 54 52 51 50 49 48 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 7 7 7 7]";
		} else if ( 5 <= p && p < 7 ) {
			precision = 17;
			buffer << "[8 8 9]\n[57 55 53 51 50 49 48 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 7 7 7 7 7]";
		} else if ( 7 <= p && p < 11 ) {
			precision = 14;
			buffer << "[8 8 7]\n[57 55 53 51 50 49 48 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 6 6 6 6]";
		} else if ( 11 <= p && p < 23 ) {
			precision = 12;
			buffer << "[7 7 6]\n[57 55 53 51 50 49 48 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 6 6 6 6]";
		} else if ( 23 <= p && p < 29 ) {
			precision = 11;
			buffer << "[6 6 6]\n[56 54 52 50 49 48 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 6 6 6 6 6]";
		} else if ( 29 <= p  ) {
			precision = 10;
			buffer << "[6 6 5]\n[56 54 52 50 49 48 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 5 5 5 5]";
        }
    } else {
        abort();
    }
    buffer >> N;
    buffer >> charpoly_prec;
}


/*
 * I have manually fixed n = 2 and d == 3
def default_args(n, d):
    p = 2;
    old_pfz = None;
    old_p = None;
    output = "";
    output +="\tif ( n == %d && d == %d ) {\n" % (n, d)
    while p <= 43:
        p = next_prime(p);
        if p == 3:
            pfz = precision_for_zeta(n, d, p, bound = 3)[-1]
        else:
            pfz = precision_for_zeta(n, d, p)[-1]
        #print p, pfz
        if old_pfz is None:
            old_pfz = pfz
        if pfz != old_pfz:
            if old_p is None:
                output += "\t\tif ( p < %d) {\n" % p
            else:
                output += "\t\t} else if ( %d <= p && p < %d ) {\n" % (old_p, p,)
            
            output += "\t\t\tprecision = %d;\n" % old_pfz[1]
            old_pfz[-1][-1] = old_pfz[-1][-2];
            output += '\t\t\tbuffer << \"%s\\n%s\";\n' % (str(old_pfz[3]).replace(",", ""), str(old_pfz[-1]).replace(",", ""))
            old_pfz = pfz
            old_p = p
    output += "\t\t} else if ( %d <= p  ) {\n" % (old_p,)
    output += "\t\t\tprecision = %d;\n" % old_pfz[1]
    old_pfz[-1][-1] = old_pfz[-1][-2];
    output += '\t\t\tbuffer << \"%s\\n%s\";\n' % (str(old_pfz[3]).replace(",", ""), str(old_pfz[-1]).replace(",", ""))
    output += "\t}"
    print output
for n in range(2,4):
    for d in range(n + 1,6):
        default_args(n, d)
*/
