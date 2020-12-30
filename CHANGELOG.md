# HEAD

-   **Breaking change**: `Alm.alm` is now enforced to be a one-dimensional
    array. This fixes type instability ([PR#25](https://github.com/ziotom78/Healpix.jl/pull/25))
    
-   **Breaking change**: Generalize the definition of `Map` and `Alm`
    ([PR#26](https://github.com/ziotom78/Healpix.jl/pull/26)) so that 
    other array types than plain `Array` can be used for these objects.

-   Documentation for new spherical harmonics functions has been added
    ([PR#35](https://github.com/ziotom78/Healpix.jl/pull/35)).
    
-   Add dependency on [Libsharp.jl](https://github.com/ziotom78/libsharp.jl)
    and implement `map2alm`, `alm2map`, and `alm2cl`
    ([#21](https://github.com/ziotom78/Healpix.jl/pull/21), 
    [#23](https://github.com/ziotom78/Healpix.jl/pull/23)). **Caution**: this
    change drops support for Julia 1.0, 1.1, and 1.2, as Libsharp.jl requires
    [Julia's artifacts](https://julialang.org/blog/2019/11/artifacts/), which
    were implemented in Julia 1.3.

- Fix incompatibility with FITSIO 1.0
  ([#32](https://github.com/ziotom78/Healpix.jl/pull/32))

- Add function `interpolate`
  ([#19](https://github.com/ziotom78/Healpix.jl/pull/19))

- Add functions `pix2zphiRing`, `pix2zphiNest`, `ringAbove`
  ([#18](https://github.com/ziotom78/Healpix.jl/pull/18))

# Version 2.3.0

- Relicense the library under GPL2
  ([#16](https://github.com/ziotom78/Healpix.jl/pull/16))

# Version 2.2.0

- Add Appveyor support to test Healpix under Windows
  ([#14](https://github.com/ziotom78/Healpix.jl/pull/14))
- Better structure of the code: now the source code is split into
  smaller files, and tests have been grouped in sets
  ([#13](https://github.com/ziotom78/Healpix.jl/pull/13))
- New type `PolarizedMap` added ([#13](https://github.com/ziotom78/Healpix.jl/pull/13))
- New keyword `write_keywords` added to `saveToFITS`
  ([#13](https://github.com/ziotom78/Healpix.jl/pull/13))
- Make `Map` descend from the abstract type `GenericMap`
  ([#12](https://github.com/ziotom78/Healpix.jl/pull/12))

# Version 2.1.0

- Two new map projections: `orthographic2` and `gnomonic` ([43e90915](https://github.com/ziotom78/Healpix.jl/commit/43e90915dba47577de322970bbc14d58b9830ab5))

# Version 2.0.0

- Use [RecipesBase](https://github.com/JuliaPlots/RecipesBase.jl) and
  [Plots](https://github.com/JuliaPlots/Plots.jl) to display maps. The use of
  Plots is not mandatory and can be avoided for headless terminals.
- Upgrade the documentation to Documenter 0.21

# Version 1.1.1

- Fix a bug in the implementation of `setindex!`

# Version 1.1.0

- Maps implement the iterator interface, so it is possible to treat a
  map like an array. This should make Healpix.jl more similar to healpy.

# Version 1.0.1

- Missing values are treated properly by `project` (using `ismissing`)
- `tod2map` properly skips `missing` values
- `tod2map` should be slightly more performant, as it uses `@inbounds`
- Keyword `numfmt` has been added to `project`, `equirectangular`,
  `mollweide`, `orthographic`. The value of the keyword can be any
  function taking a number and returning a string; the default is `x
  -> @sprintf("%g", x)`.
