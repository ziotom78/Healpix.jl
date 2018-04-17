# Map projections

export lat2colat, project, equiprojinv, mollweideprojinv, equirectangular, mollweide

import Plots: heatmap

"""
   lat2colat(x)

Convert latitude into colatitude. Both `x` and the result are expressed in radians.
"""
lat2colat(x) = π / 2 - x

"""
    project(m::Map{T, O}; kwargs...) where {T, O <: Order}

Return a 2D bitmap (array) containing a cartographic projection of the map and a
2D bitmap containing a boolean mask. The size of the bitmap is specified by
`figsize`, which must be a 2-tuple. The function `projfn` must be a function
which accepts as input two parameters `x` and `y` (numbers between -1 and 1).

The following keywords can be used in the call:

- `figsize`: 2-tuple specifying the (height, width) of the bitmap in pixels
- `center`: 2-tuple specifying the location (colatitude, longitude) of the sky
  point that is to be placed in the middle of the image (in radians)
- `show`: Boolean; if true (the default), the bitmap will be shown using
  functions from the "Plots" package.
- `returnmask`: Boolean; if true, the function returns a 2-tuple containing
  the image bitmap and a mask bitmap, which is set to true if the pixel falls
  within the carthographic projection, false otherwise. If `returnmask` is
  false, only the image bitmap is returned.
- `show`: Boolean. If true (the default), the map will be displayed. It has no
  effect if `returnmask` is true.  
"""
function project(invprojfn, m::Map{T,O}; kwargs...) where {T <: AbstractFloat, O <: Order}

    args = Dict(kwargs)
    figsize = get(args, :figsize, (400, 800))
    center = get(args, :center, (0, 0))
    returnmask = get(args, :returnmask, false)
    show = get(args, :show, true)

    image = zeros(eltype(m.pixels), figsize)
    mask = Array{Bool}(figsize)
    
    for j = 1:size(image)[1]
        y = 2(j - 1) / (size(image)[1] - 1) - 1
        for i = 1:size(image)[2]
            x = 2(i - 1) / (size(image)[2] - 1) - 1
            (flag, lat, long) = invprojfn(x, y; kwargs...)
            
            if flag
                image[j, i] = m.pixels[Healpix.ang2pix(m, lat2colat(lat), long)]
            else
                image[j, i] = convert(T, NaN)
            end
    
            mask[j, i] = flag
        end
    end

    if returnmask
        (image, mask)
    else
        if show
            heatmap(image, aspectratio = 1, xaxis = false, yaxis = false)
        else
            image
        end
    end
end

"""
    function equiprojinv(x, y)

Inverse equirectangular projection. Given a point (x, y)
on the plane [-1, 1] × [-1, 1], return a tuple (Bool, Number, Number)
where the first Boolean is a flag telling if the point falls
within the projection (true) or not (false), and the two numbers
are the latitude and colatitude in radians.
"""
function equiprojinv(x, y; kwargs...)
    ((-1 ≤ x ≤ 1) && (-1 ≤ y ≤ 1)) || return (false, 0, 0)
     
    (true, π / 2 * y, π * x)
end

"""
    function mollweideprojinv(x, y)

Inverse Mollweide projection. Given a point (x, y) on the plane,
with x ∈ [-1, 1], y ∈ [-1, 1], return a 3-tuple of type
(Bool, Number, Number). The boolean specifies if (x, y) falls within
the map (true) or not (false), the second and third arguments are
the latitude and longitude in radians.
"""
function mollweideprojinv(x, y; kwargs...)
    # See https://en.wikipedia.org/wiki/Mollweide_projection, we set
    #
    #     R = 1/√2
    #
    # x ∈ [-1, 1], y ∈ [-1, 1]

    if x^2 + y^2 > 1
        return (false, 0, 0)
    end

    sinθ = y
    cosθ = sqrt(1 - sinθ^2)
    θ = asin(sinθ)

    lat = asin((2θ + 2 * sinθ * cosθ) / π)
    long = -2 * π * x / (2cosθ)
    (true, lat, long)

end

"""
    equirectangular(m::Map{T,O}; kwargs...) where {T <: AbstractFloat, O <: Order}

High-level wrapper around `project` for equirectangular projections.
"""
equirectangular(m::Map{T,O}; kwargs...) where {T <: AbstractFloat, O <: Order} = project(equiprojinv, m; kwargs...)

"""
    mollweide(m::Map{T,O}; kwargs...) where {T <: AbstractFloat, O <: Order}

High-level wrapper around `project` for Mollweide projections.
"""
mollweide(m::Map{T,O}; kwargs...) where {T <: AbstractFloat, O <: Order} = project(mollweideprojinv, m; kwargs...)
