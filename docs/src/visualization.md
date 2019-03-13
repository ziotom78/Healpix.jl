# Visualization functions

Healpix.jl uses RecipesBase to display maps. You need to import
`Plots` in order to display maps, using the `plot` functions.  Maps
are internally treated as heatmaps, so your backend should support
this kind of visualization: at the moment, this is true for GR, PlotLy
and PyPlot.


```@example
using Healpix
using Plots
gr()  # Use the GR backend

nside = 8
m = Map{Float64, RingOrder}(nside)
m.pixels[:] = 1:length(m.pixels)
plot(m)
```

A call to `plot` can provide two additional arguments:

1. A carthographic projection (see below).
2. A dictionary containing parameters to be used by the carthographic
   projection.
   
The following example shows the map in orthographic coordinates:

```@example
plot(m, orthographic)
```

## Cartographic projections

```@docs
mollweide
equirectangular
orthographic
```
