```@meta
DocTestSetup = quote
    using Healpix
end
```

# Healpix.jl: an implementation of the Healpix tessellation scheme in Julia

This is the documentation of the [Healpix.jl](https://github.com/ziotom78/Healpix.jl) package, an implementation of the Healpix spherical tessellation scheme written entirely in Julia.

This package has the main purpose of providing a Julia-only solution, so that it can easily be used on platforms not supported by the Healpix C++ library (e.g., Windows).
This library implements algorithms for converting directions into pixel indices and vice versa. It supports both `RING` and `NESTED` schemes, and it employs Julia's powerful type system to avoid mistaking one scheme in place of the other.

## Index

```@index
```
