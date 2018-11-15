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
