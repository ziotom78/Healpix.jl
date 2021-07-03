```@meta
DocTestSetup = quote
    using Healpix
end
```

# Map functions

Functions like [`pix2angNest`](@ref) and [`ang2pixNest`](@ref) fully
define the Healpix tessellation scheme. They are however extremely
impractical in a number of situations. It happens often that a large
fraction of pixels in a map need to be processed together. Healpix.jl
introduces the [`HealpixMap{T, O <: Order}`](@ref) type, which acts as a
collection of all the pixels on the sphere. A `HealpixMap` type holds the
value of all the pixels in its `pixels` field, and it keeps track of
the ordering (either `RING` or `NESTED`). Here is an example that
shows how to create a map and initialize it:

```julia
nside = 32
m = HealpixMap{Float64, RingOrder}(nside)
m.pixels[:] = 1.0  # Set all pixels to 1
```

Healpix.jl defines the basic operations on maps (sum, subtraction,
multiplication, division). These operations can either combine two
maps or a map and a scalar value:

```julia
mollweide(m * 2.0)
mollweide(m * m)
```

The [`HealpixMap{T, O <: Order}`](@ref) is derived from the abstract
type [`AbstractHealpixMap{T}`](@ref), which does not encode the
ordering. It is useful for functions that can either work on
ring/nested-ordered maps but cannot be executed on plain generic
arrays:

```julia
# Return the number of pixels in the map, regardless of its ordering
maplength(m::AbstractHealpixMap{T}) where T = length(m)

# This returns 12
maplength(HealpixMap{Float64, RingOrder}(1))

# This too returns 12
maplength(HealpixMap{Float64, NestedOrder}(1))

# This fails
maplength(zeros(Float64, 12))
```

Healpix.jl implements the [`PolarizedHealpixMap{T, O <: Order}`](@ref)
type as well, which derives from
[`AbstractPolarizedHealpixMap{T}`](@ref). This encodes three maps
containing the I/Q/U signal: the intensity (I), and the Q and U Stokes
parameters. The three maps must have the same resolution.

```@docs
AbstractHealpixMap
HealpixMap
AbstractPolarizedHealpixMap
PolarizedHealpixMap
```

## Unseen pixels and nothingness

You can use the constant [`UNSEEN`](@ref) to mark unseen pixels, i.e.,
pixels that lack data associated with them, in a way that is
compatible with other Healpix libraries.

```@example
m = Healpix.HealpixMap{Float32, RingOrder}(32)

# Mark all the pixels in the map as «unseen» (missing)
m[:] .= UNSEEN
```

However, Julia provides a sounder way to denote missing pixels through
the use of `Nothing` (type) and `nothing` (value). Whenever you pass a
`Union{Nothing, T}` type to a Healpix map, the map will be initialized
to `nothing`, and you can test if a pixel has been observer or not
using Julia's function `isnothing`:

```@example
m = Healpix.HealpixMap{Union{Int32, Nothing}, Healpix.RingOrder}(1)
m[:] = 1:12
m[5] = nothing
@assert isnothing(m[5])
```

Note that, unlike [`UNSEEN`](@ref), this mechanism permits to signal
«missing» pixels even for maps that do not use floating-point numbers.

!!! warning "Reading/saving maps with `nothing` values"

    You can use [`saveToFITS`](@ref) or [`readMapFromFITS`](@ref) on
    maps whose base type is `Union{Nothing, T}` only if `T` is a
    floating-point number, because for the sake of compatibility with
    other Healpix libraries the FITS file will use [`UNSEEN`](@ref) to
    mark missing values.

## Encoding the order

Healpix.jl distinguishes between `RING` and `NEST` orderings using
Julia's typesystem. The abstract type `Order` has two descendeants,
`RingOrder` and `NestedOrder`, which are used to instantiate objects
of type `HealpixMap`. Applying the functions [`nest2ring`](@ref) and 
[`ring2nest`](@ref) to maps converts those maps to the appropriate orders.
In-place [`nest2ring!`](@ref) and [`ring2nest!`](@ref) versions are also 
available.

```@docs
Order
RingOrder
NestedOrder
nest2ring(m_nest::HealpixMap{T, NestedOrder, AA}) where {T, AA}
ring2nest(m_ring::HealpixMap{T, RingOrder, AA}) where {T, AA}
nest2ring!(m_ring_dst::HealpixMap{T, RingOrder, AAR}, 
  m_nest_src::HealpixMap{T, NestedOrder, AAN}) where {T, AAN, AAR}
ring2nest!(m_nest_dst::HealpixMap{T, NestedOrder, AAN}, 
  m_ring_src::HealpixMap{T, RingOrder, AAR}) where {T, AAR, AAN}
```



## Pixel functions

When working with maps, it is not needed to pick between
[`ang2pixNest`](@ref) and [`ang2pixRing`](@ref) because a `HealpixMap` type
already encodes the ordering. Functions `pix2ang` and `ang2pix` always
choose the correct ordering, but they require a `HealpixMap` instead of a
[`Resolution`](@ref) as their first argument.

```@docs
pix2ang
ang2pix
```

## Loading and saving maps

Healpix.jl implements a number of functions to save maps in FITS files.

```@docs
saveToFITS
```

Function `savePixelsToFITS` is a low-level function. It knows nothing
about the ordering schema used for the pixels, so the caller should
manually write the `ORDERING` keyword in the HDU header by itself.

```@docs
savePixelsToFITS
```

To load a map from a FITS file, you can either use
[`readMapFromFITS`](@ref) or [`readPolarizedMapFromFITS`](@ref).

```@docs
readMapFromFITS
readPolarizedMapFromFITS
```

## Testing for conformability

It often happens that two Healpix maps need to be combined together:
for instance, pixels on a sky map might need to be masked using a sky
mask, or one map might need to be subtracted from another
one. «Conformability» means that the operation between the two maps
can be done directly on the pixels, without oordering or resolution
conversions. The function `conformables` checks this.

```@repl
using Healpix # hide
m1 = HealpixMap{Float64, RingOrder}(1)
m2 = HealpixMap{Float64, RingOrder}(1)
m3 = HealpixMap{Float64, NestedOrder}(1)
m4 = HealpixMap{Float64, NestedOrder}(2)
conformables(m1, m2)
conformables(m1, m3)
conformables(m1, m4)
```

```@docs
conformables
```

## Interpolation

The fact that a Healpix map is, well, pixelized means that there is a
sharp boundary between adjacent pixels. This can lead to undesidable
effects, and therefore Healpix.jl provides a function,
[`interpolate`](@ref), that returns the interpolated value of the map
along some direction in the sky:

- If the direction (θ, ɸ) passes through the center of a pixel,
  `interpolate` returns the value of the pixel itself;
- Otherwise, the value of the pixel and its neighbours will be
  interpolated using a linear function to return a weighted value.

```@docs
interpolate
```

## Upgrading and Downgrading
Changing resolution is done with [`udgrade`](@ref). This is very fast for nested orderings,
but slow for ring ordering.

```@docs
udgrade
```

## Map-making

Map-making is the process of converting a time series of measurements
into a sky map. The most basic form of map-making is the so-called
"binning", where samples in the time stream falling within the same
sky pixel are averaged. This map-making algorithm is strictly accurate
only if the noise in the time stream is white.

Healpix.jl implements two functions to perform binning,
[`tod2map`](@ref) and [`combinemaps!`](@ref).

```@docs
tod2map
combinemaps!
```
