# HEAD

-   Add dependency on LazyArtifacts [18fde42](https://github.com/ziotom78/Healpix.jl/commit/18fde420dba069277dfb328076c99a995bc5c275)

# Version 4.2.3

-   Improve the documentation [#124](https://github.com/ziotom78/Healpix.jl/pull/124) [thanks to @abhro]

-   Fix a type stability bug [#125](https://github.com/ziotom78/Healpix.jl/pull/125) [thanks to @hsgg]

# Version 4.2.2

-   Fix bug [#112](https://github.com/ziotom78/Healpix.jl/issues/112)

# Version 4.2.1

-   Fix bug [#108](https://github.com/ziotom78/Healpix.jl/issues/108) in `pix2ringpos` [#109](https://github.com/ziotom78/Healpix.jl/pull/109)

# Version 4.2.0

-   Add overloads of `Base`'s algebraic operators `+`, `-`, `*`, `/` and `LinearAlgebra`'s `dot` product in Alm space; add Alm indexing functions: `each_ell`, `each_ell_idx`, `each_m`, `each_m_idx`, `each_ell_m`; improve `almxfl!`'s performance; add functions `getEquatorIdx`, `ring2theta`, `getRingPixels` [#98](https://github.com/ziotom78/Healpix.jl/pull/98), thanks to [LeeoBianchi](https://github.com/LeeoBianchi)

-   Add `synalm`, `synfast`, `anafast`, `adjoint_alm2map`, `adjoint_map2alm` [#96](https://github.com/ziotom78/Healpix.jl/pull/96), thanks to [LeeoBianchi](https://github.com/LeeoBianchi)

-   Add functions `readClFromFITS`, `writeClToFITS`, `cl2dl`, `dl2cl`, `almxfl`, `almxfl!`, `almExplicitIndex` [#91](https://github.com/ziotom78/Healpix.jl/pull/91), thanks to [LeeoBianchi](https://github.com/LeeoBianchi)

-   Fix a bug in `mollweideproj` and `equiproj` [#97](https://github.com/ziotom78/Healpix.jl/issues/97)

-   Fix a bug in `queryDiscRing` for directions close to the poles [#105](https://github.com/ziotom78/Healpix.jl/issues/105)

-   Re-enable support for Julia 1.9 [#104](https://github.com/ziotom78/Healpix.jl/pull/104)

# Version 4.1.2

-   Use double precision in `ang2pixRing` and `zphi2pixRing` [#94](https://github.com/ziotom78/Healpix.jl/pull/94)

-   Improve CI builds [#89](https://github.com/ziotom78/Healpix.jl/pull/89)

# Version 4.1.1

-   Bugfixes: [#87](https://github.com/ziotom78/Healpix.jl/pull/87)

# Version 4.1.0

- New function `gaussbeam` [#85](https://github.com/ziotom78/Healpix.jl/pull/85)

# Version 4.0.1

-   Add checks for colatitude θ in functions that take angles as input [#84](https://github.com/ziotom78/Healpix.jl/pull/84)

# Version 4.0.0

-   **Breaking change**: `udgrade` now always assumes `pess = false` [#72]()

-   Implement `Base.parent` for `HealpixMap` which returns the underlying array

-   New functions `boundariesRing!` and `boundariesRing` [#81](https://github.com/ziotom78/Healpix.jl/pull/81)

-   New constants/functions: `ORDER_MAX`, `nside2order`, `order2nside`, `numOfRings`, `ring2z`, `max_pixrad`, `queryDiscRing`, `queryStripRing` [#80](https://github.com/ziotom78/Healpix.jl/pull/80)

-   Use a more accurate algorithm for `vec2ang` [#76][https://github.com/ziotom78/Healpix.jl/pull/76]

-   Make `NSIDE_MAX` dependent on the architecture (32/64 bit) of the system [#75](https://github.com/ziotom78/Healpix.jl/pull/75)

-   Documentation improvements [#82](https://github.com/ziotom78/Healpix.jl/pull/82)


# Version 3.0.1

-   Fix the documentation and export `readPolarizedMapFromFITS` [#69](https://github.com/ziotom78/Healpix.jl/issues/67)

-   Bump dependency on Libsharp.jl to fix issue [#67](https://github.com/ziotom78/Healpix.jl/issues/67)

# Version 3.0.0

-   **Breaking change**: rename `Map` → `HealpixMap`, `PolarizedMap` → `HealpixPolarizedMap`, `GenericMap` → `AbstractHealpixMap` ([PR#53](https://github.com/ziotom78/Healpix.jl/pull/53))

-   **Breaking change**: `Alm.alm` is now enforced to be a one-dimensional array. This fixes type instability ([PR#25](https://github.com/ziotom78/Healpix.jl/pull/25))

-   **Breaking change**: Generalize the definition of `HealpixMap` and `Alm` ([PR#26](https://github.com/ziotom78/Healpix.jl/pull/26)) so that other array types than plain `Array` can be used for these objects.

-   Add support for maps whose base type is `Union{Nothing, T}` [PR#63](https://github.com/ziotom78/Healpix.jl/pull/63)

-   Various documentation improvements [#55](https://github.com/ziotom78/Healpix.jl/pull/55)

-   Fix issue [#61](https://github.com/ziotom78/Healpix.jl/pull/61)

-   Fix issue [#59](https://github.com/ziotom78/Healpix.jl/pull/59)

-   Fix issue [#57](https://github.com/ziotom78/Healpix.jl/issues/57)

-   Use [CFITSIO.jl](https://github.com/JuliaAstro/CFITSIO.jl) instead of [FITSIO](https://juliaastro.github.io/FITSIO.jl/stable/) ([PR#50](https://github.com/ziotom78/Healpix.jl/pull/50), it fixes [#47](https://github.com/ziotom78/Healpix.jl/issues/47))

-   Add untyped constructor for `PolarizedMap` ([PR#49](https://github.com/ziotom78/Healpix.jl/pull/49))

-   Implement `udgrade` ([PR#46](https://github.com/ziotom78/Healpix.jl/pull/46), [#51](https://github.com/ziotom78/Healpix.jl/pull/51))

-   Implement `pixwin` ([PR#45](https://github.com/ziotom78/Healpix.jl/pull/45))

-   Implement `applyFullWeights!` ([PR#41](https://github.com/ziotom78/Healpix.jl/pull/41), [#44](https://github.com/ziotom78/Healpix.jl/pull/44)), and [#52](https://github.com/ziotom78/Healpix.jl/pull/52)

-   Implement `ring2nest!` and `nest2ring!` ([PR#40](https://github.com/ziotom78/Healpix.jl/pull/40))

-   Make `ang2vec` return a tuple instead of a list ([PR#38](https://github.com/ziotom78/Healpix.jl/pull/38))

-   Remove `Manifest.toml` ([PR#42](https://github.com/ziotom78/Healpix.jl/pull/42))

-   README has been updated ([PR#39](https://github.com/ziotom78/Healpix.jl/pull/39))

-   Documentation for new spherical harmonics functions has been added ([PR#35](https://github.com/ziotom78/Healpix.jl/pull/35)).

-   Add dependency on [Libsharp.jl](https://github.com/ziotom78/libsharp.jl) and implement `map2alm`, `alm2map`, and `alm2cl` ([#21](https://github.com/ziotom78/Healpix.jl/pull/21), [#23](https://github.com/ziotom78/Healpix.jl/pull/23)). **Caution**: this change drops support for Julia 1.0, 1.1, and 1.2, as Libsharp.jl requires [Julia's artifacts](https://julialang.org/blog/2019/11/artifacts/), which were implemented in Julia 1.3.

-   Fix incompatibility with FITSIO 1.0 ([#32](https://github.com/ziotom78/Healpix.jl/pull/32))

-   Add function `interpolate` ([#19](https://github.com/ziotom78/Healpix.jl/pull/19))

-   Add functions `pix2zphiRing`, `pix2zphiNest`, `ringAbove` ([#18](https://github.com/ziotom78/Healpix.jl/pull/18))

# Version 2.3.0

-   Relicense the library under GPL2 ([#16](https://github.com/ziotom78/Healpix.jl/pull/16))

# Version 2.2.0

-   Add Appveyor support to test Healpix under Windows ([#14](https://github.com/ziotom78/Healpix.jl/pull/14))
-   Better structure of the code: now the source code is split into smaller files, and tests have been grouped in sets ([#13](https://github.com/ziotom78/Healpix.jl/pull/13))
-   New type `PolarizedMap` added ([#13](https://github.com/ziotom78/Healpix.jl/pull/13))
-   New keyword `write_keywords` added to `saveToFITS` ([#13](https://github.com/ziotom78/Healpix.jl/pull/13))
-   Make `Map` descend from the abstract type `GenericMap` ([#12](https://github.com/ziotom78/Healpix.jl/pull/12))

# Version 2.1.0

-   Two new map projections: `orthographic2` and `gnomonic` ([43e90915](https://github.com/ziotom78/Healpix.jl/commit/43e90915dba47577de322970bbc14d58b9830ab5))

# Version 2.0.0

-   Use [RecipesBase](https://github.com/JuliaPlots/RecipesBase.jl) and [Plots](https://github.com/JuliaPlots/Plots.jl) to display maps. The use of Plots is not mandatory and can be avoided for headless terminals.
-   Upgrade the documentation to Documenter 0.21

# Version 1.1.1

- Fix a bug in the implementation of `setindex!`

# Version 1.1.0

- Maps implement the iterator interface, so it is possible to treat a map like an array. This should make Healpix.jl more similar to healpy.

# Version 1.0.1

-   Missing values are treated properly by `project` (using `ismissing`)
-   `tod2map` properly skips `missing` values
-   `tod2map` should be slightly more performant, as it uses `@inbounds`
-   Keyword `numfmt` has been added to `project`, `equirectangular`, `mollweide`, `orthographic`. The value of the keyword can be any function taking a number and returning a string; the default is `x -> @sprintf("%g", x)`.
