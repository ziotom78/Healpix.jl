# Visualization functions

Healpix.jl implements the [`project`](@ref) function, which creates a 2D matrix containing the cartographic projection of a map. A few standard cartographic projections are implemented, but users can provide their own projections. The function can optionally use `heatmap` (from the Plots.jl package) to display the 2D matrix. Two useful wrappers to `project` are [`equirectangular`](@ref) and [`mollweide`](@ref), which employ the equirectangular and Mollweide projections respectively.

```@example
using Healpix
nside = 8
m = Map{Float64, RingOrder}(nside)
m.pixels[:] = 1:length(m.pixels)
mollweide(m)
```

```@docs
project
```

## Cartographic projections

```@docs
equiprojinv
mollweideprojinv
```

## High-level wrappers to `project`

```@docs
equirectangular
mollweide
```
