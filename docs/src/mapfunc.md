```@meta
DocTestSetup = quote
    using Healpix
end
```

# Map functions

Functions like [`pix2angNest`](@ref) and [`ang2pixNest`](@ref) fully define the Healpix tessellation scheme. They are however extremely impractical in a number of situations. It happens often that a large fraction of pixels in a map need to be processed together. Healpix.jl introduces the [`Map{T, O <: Order}`](@ref) type, which acts as a collection of all the pixels on the sphere. A `Map` type holds the value of all the pixels in its `pixels` field, and it keeps track of the ordering (either `RING` or `NESTED`). Here is an example that shows how to create a map and initialize it:

```julia
nside = 32
m = Map{Float64, RingOrder}(nside)
m.pixels[:] = 1.0  # Set all pixels to 1
```

```@docs
Map
```

## Encoding the order

Healpix.jl distinguishes between `RING` and `NEST` orderings using Julia's typesystem. The abstract type `Order` has two descendeants, `RingOrder` and `NestedOrder`, which are used to instantiate objects of type `Map`.

```@docs
Order
RingOrder
NestedOrder
```

## Pixel functions

When working with maps, it is not needed to pick between [`ang2pixNest`](@ref) and [`ang2pixRing`](@ref) because a `Map` type already encodes the ordering. Functions `pix2ang` and `ang2pix` always choose the correct ordering, but they require a `Map` instead of a [`Resolution`](@ref) as their first argument.

```@docs
pix2ang
ang2pix
```

## Loading and saving maps

Healpix.jl implements a number of functions to save maps in FITS files.

```@docs
saveToFITS
```

Function `savePixelsToFITS` is a low-level function. It knows nothing about the ordering schema used for the pixels, so the caller should manually write the `ORDERING` keyword in the HDU header by itself.

```@docs
savePixelsToFITS
```

To load a map from a FITS file, you can use `readMapFromFITS`.

```@docs
readMapFromFITS
```

## Testing for conformability

It often happens that two Healpix maps need to be combined together: for instance, pixels on a sky map might need to be masked using a sky mask, or one map might need to be subtracted from another one. «Conformability» means that the operation between the two maps can be done directly on the pixels, without oordering or resolution conversions. The function `conformables` checks this.

```@repl
using Healpix # hide
m1 = Map{Float64, RingOrder}(1)
m2 = Map{Float64, RingOrder}(1)
m3 = Map{Float64, NestedOrder}(1)
m4 = Map{Float64, NestedOrder}(2)
conformables(m1, m2)
conformables(m1, m3)
conformables(m1, m4)
```

```@docs
conformables
```

## Map-making

Map-making is the process of converting a time series of measurements into a sky map. The most basic form of map-making is the so-called "binning", where samples in the time stream falling within the same sky pixel are averaged. This map-making algorithm is strictly accurate only if the noise in the time stream is white.

Healpix.jl implements two functions to perform binning, [`tod2map`](@ref) and [`combinemaps`](@ref).

```@docs
tod2map
combinemaps
```