# controlled-reduction

[![Build Status](https://travis-ci.org/edgarcosta/controlledreduction.svg?branch=master)](https://travis-ci.org/edgarcosta/controlledreduction)

## Abstract
An implementation of the controlled reduction method for computing the
Hasse-Weil zeta functions of smooth projective hypersurfaces over finite
fields. Explicitly, by computing a p-adic approximation of Frobenius
action on p-adic cohomology (Monsky-Washnizter) with sufficient precision
and then lifting its characteristic polynomial to the integers.
An overview of the method can be found in: 
 - "Computing zeta functions of nondegenerate projective hypersurfaces over 
finite fields", (under preparation) by Edgar Costa and David Harvey.
 - "Effective computations of Hasse-Weil zeta functions", by Edgar Costa

## Dependencies
It majorly depends on:
 - [NTL: A Library for doing Number Theory](http://www.shoup.net/ntl/)
 - [FLINT: Fast Library for Number Theory](http://flintlib.org/)
 
Which depend on:

 - [GMP: GNU Multiple Precision Arithmetic Library](https://gmplib.org/) (for NTL and FLINT)
 - [MPIR: Multiple Precision Integers and Rationals](mpir.org) (for FLINT)
 - [MPFR: GNU Multiple Precision Floating-Point Reliably](http://www.mpfr.org/) (for FLINT)

However, [SageMath](http://www.sagemath.org/) comes with all this libraries. 


## Installation

There are 2 options:

### Using SageMath to provide the dependencies

1. Figuring out where `SageMath` is installed. 
We recommend doing this and storing in an environmental variable by doing:
```
SAGE_ROOT=$(sage -c "print SAGE_ROOT")
```
Alternatively, in Sage do `print SAGE_ROOT` and on the unix terminal:
`
SAGE_ROOT=<the line printed in Sage>
`

2. Download controlled reduction
```
git clone https://github.com/edgarcosta/controlledreduction.git
```

3. Change your working directory and run the configure file
```
cd controlledreduction && ./configure --with-ntl=$SAGE_ROOT/local
```

4. Compile everything by doing
```
make
```

5. Set up the variable `$LD_LIBRARY_PATH`, so the executables can find the libraries at run time.
One can do this for the current terminal by doing:
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SAGE_ROOT/local/lib
```
and for this line to be ran at the start of every session one can do
```
echo export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SAGE_ROOT/local/lib >> ~/.bashrc
```


6. You can optionaly run some tests by doing
```
make check
```
or check out some of the examples
```
build/examples/K3_dwork
```

### Manually specifying the path for the dependencies 


1. make sure you have the dependencies (see the source of the script `build_dependencies` if you would prefer to build them manually), if they are installed in a non-standard path, be sure to set  `$LD_LIBRARY_PATH` accordingly.

2. `./configure` to generate the makefile.

   To link against the libaries by SageMath it should be sufficient to run `./configure --with-ntl=<SAGE_DIR>/local/`.
 
   Run `./configure --help` for more options.

3. `make` to build everything

4. `make check` to run some tests. 


