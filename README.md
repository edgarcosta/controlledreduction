# controlled-reduction

[![Build Status](https://travis-ci.org/edgarcosta/controlled-reduction.svg?branch=master)](https://travis-ci.org/edgarcosta/controlled-reduction)

An implementation of the controlled reduction method for computing the
Hasse-Weil zeta functions of smooth projective hypersurfaces over finite
fields. Explicitly, by computing a p-adic approximation of Frobenius
action on p-adic cohomology (Monsky-Washnizter) with sufficient precision
and then lifting its characteristic polynomial to the integers.
An overview of the method can be found in: 
 - "Computing zeta functions of nondegenerate projective hypersurfaces over 
finite fields", (under preparation) by Edgar Costa and David Harvey.
 - "Effective computations of Hasse-Weil zeta functions", by Edgar Costa


It majorly depends on:
 - [NTL: A Library for doing Number Theory](http://www.shoup.net/ntl/)
 - [FLINT: Fast Library for Number Theory](http://flintlib.org/)
 
Which depend on:

 - [GMP: GNU Multiple Precision Arithmetic Library](https://gmplib.org/) (for NTL and FLINT)
 - [MPIR: Multiple Precision Integers and Rationals](mpir.org) (for FLINT)
 - [MPFR: GNU Multiple Precision Floating-Point Reliably](http://www.mpfr.org/) (for FLINT)

[SageMath](http://www.sagemath.org/) comes with all this libraries. 

Alternatively, you might use the script `build_dependencies` to build all these dependencies in a UNIX system. 
You can change the place where these dependencies will be installed by setting the environment variable `$PREFIX`. See the source for more details.

To install:

1. make sure you have the dependencies, if they are installed in a non-standard path, be sure to set  `$LD_LIBRARY_PATH` accordingly.

2. `./configure` to generate the makefile.

   To link against the libaries by SageMath it should be sufficient to run `./configure --with-ntl=<SAGE_DIR>/local/`.
 
   Run `./configure --help` for more options.

3. `make` to build everything

4. `make check` to run some tests. 


