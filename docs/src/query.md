```@meta
DocTestSetup = quote
    using Healpix
end
```

# Query functions

It is often useful to perform calculations on a set of adjacent pixels
in a map. Asking for the indices of pixels within a given region of
the sphere is called a *query*, and this kind of function can be
implemented efficiently using the Healpix scheme.

The functions that query a set of pixels satisfying some proximity criteria are the following:

- [`queryDiscRing`](@ref)
- [`queryStripRing`](@ref)


Let's show an example about how to use `queryDiscRing`, which returns
a list of the indexes of the pixels that fall within some angle from a
direction on the sky sphere:

```@example querydiscexample
using Healpix # hide
resol = Resolution(32)
(theta, phi) = (1.3, 0.7)
radius = deg2rad(10.0)
pixidx = queryDiscRing(resol, theta, phi, radius)
```

We can visualize where these pixels are using a Mollweide projection:

```@example querydiscexample
m = HealpixMap{Float32, RingOrder}(resol.nside);
m[pixidx] .= 1;

m[ang2pix(m, theta, phi)] = 2  # Highlight the pixel at the center

using Plots
gr()  # hide
plot(m)
savefig(joinpath("images", "querydisc.png")) # hide
nothing # hide
```
![](images/querydisc.png)

The function `queryDiscRing` accepts an optional keyword `fact` that
specifies the resolution to be used in computing the result; it can be
any positive integer, and the actual resolution is `fact * NSIDE`.

# Reference

```@docs
queryDiscRing
queryStripRing
```
