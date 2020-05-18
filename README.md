# Healpix

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ziotom78.github.io/Healpix.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ziotom78.github.io/Healpix.jl/dev)
[![Build Status](https://travis-ci.com/ziotom78/Healpix.jl.svg?branch=master)](https://travis-ci.com/ziotom78/Healpix.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/ziotom78/Healpix.jl?svg=true)](https://ci.appveyor.com/project/ziotom78/Healpix-jl)
[![Codecov](https://codecov.io/gh/ziotom78/Healpix.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ziotom78/Healpix.jl)

A set of Julia functions that implement the Healpix spherical
projection.

The purpose of this package is to implement a Julia-only library,
instead of providing bindings to the original implementation of the
[C/C++/Fortran Healpix library](http://healpix.jpl.nasa.gov/). This
should make the package easier to port to those architecture which are
not supported by the original implementation (e.g., Windows).

Of course, this mean that implementing new features is generally not
trivial, as one has to re-implement the algorithm instead of just
figuring out how to bind some C/C++/Fortran function.

The
[mweastwood/LibHealpix.jl](https://github.com/mweastwood/LibHealpix.jl)
library provides straight bindings to the original C++ library. For
the reasons stated above, it is able to provide a wider set of
functions.

## Installation

From the Julia REPL, run

````julia
import Pkg
Pkg.add("Healpix.jl")
````

## Usage examples

Here are some code snippets that show how to use `Healpix.jl`. It is
interesting to have a look at
[test/runtests.jl](https://github.com/ziotom78/Healpix.jl/blob/master/test/runtests.jl)
as well.

Refer to the
[documentation](https://ziotom78.github.io/Healpix.jl/stable) for more
examples.

### Dealing with resolutions

The resolution of a Healpix map is uniquely determined by the `NSIDE`
parameter. Healpix.jl precalculates a number of values from `NSIDE` to
save time during computations; such values are saved in a
`Healpix.Resolution` object:

`````julia
import Healpix
res = Healpix.Resolution(256)
print("The pixel index is $(Healpix.ang2pixNest(res, 0.1, 0.2))\n")
`````

### Reading a map from a FITS file

This snippet loads a map named `planck_70GHz.fits` into an array of
64-bit floating-point numbers:

`````julia
import Healpix

m = Healpix.readMapFromFITS("planck_70GHz.fits", 1, Float64)
print("average: $(mean(m.pixels))\n")
`````

## License

Healpix.jl is released under the GPL license. Versions before 2.3.0
were released under a MIT license, but this was considered
incompatible with the way the code has been written
([#15](https://github.com/ziotom78/Healpix.jl/issues/15)).
