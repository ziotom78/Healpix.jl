# HEAD

- Two new map projections: `orthographic2` and `gnomonic`
- Make `Map` descend from the abstract type `GenericMap`
  ([#12](https://github.com/ziotom78/Healpix.jl/pull/12))

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
