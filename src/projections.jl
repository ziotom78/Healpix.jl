# Map projections

export UNSEEN,
    lat2colat,
    colat2lat,
    project,
    equiproj,
    equiprojinv,
    mollweideproj,
    mollweideprojinv,
    orthoinv,
    ortho2inv,
    equirectangular,
    mollweide,
    orthographic,
    orthographic2,
    gnomonic,
    gnominv

import RecipesBase


"""
    project(invprojfn, m::HealpixMap{T, O, AA}, bmpwidth, bmpheight; kwargs...)

Return a 2D bitmap (array) containing a cartographic projection of the
map and a 2D bitmap containing a boolean mask. The size of the bitmap
is `bmpwidth`×`bmpheight` pixels. The function `projfn` must be a
function which accepts as input two parameters `x` and `y` (numbers
between -1 and 1).

The following keywords can be used in the call:

- `center`: 2-tuple specifying the location (colatitude, longitude) of the sky
  point that is to be placed in the middle of the image (in radians)
- `unseen`: by default, Healpix maps use the value -1.6375e+30 to mark
  unseen pixels. You can specify a different value using this
  keyword. This should not be used in common applications.

Return a `Array{Union{Missing, Float32}}` containing the intensity of
each pixel. Pixels falling outside the projection are marked as NaN,
and unseen pixels are marked as `missing`.
"""
function project(
    invprojfn,
    m::HealpixMap{T,O,AA},
    bmpwidth,
    bmpheight,
    projparams = Dict(),
) where {T<:Number,O,AA}

    center = get(projparams, :center, (0, 0))
    unseen = get(projparams, :unseen, UNSEEN)
    desttype = get(projparams, :desttype, Float32)

    img = Matrix{desttype}(undef, bmpheight, bmpwidth)
    masked = zeros(Bool, bmpheight, bmpwidth)

    anymasked = false
    for j = 1:bmpheight
        y = 2 * (j - 1) / (bmpheight - 1) - 1
        for i = 1:bmpwidth
            x = 2 * (i - 1) / (bmpwidth - 1) - 1
            visible, lat, long = invprojfn(x, y)
            if visible
                value = m.pixels[Healpix.ang2pix(m, lat2colat(lat), long)]
                if ismissing(value) ||
                   isnan(value) ||
                   (!ismissing(unseen) && unseen == value)
                    img[j, i] = NaN
                    masked[j, i] = true
                    anymasked = true
                else
                    img[j, i] = value
                end
            else
                img[j, i] = NaN
            end
        end
    end

    img, masked, anymasked
end

################################################################################

lat2colat(x) = π / 2 - x
colat2lat(x) = π / 2 - x

@doc raw"""
    lat2colat(x)
    colat2lat(x)

Convert colatitude into latitude and vice versa. Both `x` and the
result are expressed in radians.
"""
lat2colat, colat2lat

################################################################################

"""
    equiproj(lat, lon)

Equirectangular projection. Given the latitude `lat` (in radians) and the
longitude (in radians), return a tuple (Bool, Number, Number) where the
first Boolean is a flag telling if the point falls within the projection (true)
or not (false), and the two numbers are the x and y coordinates of the point
on the projection plane (both are in the range [−1, 1]).
"""
function equiproj(lat, lon)
    # We use `rem2pi` because we need angles in the range [-π, +π]
    x, y = (
        rem2pi(lon, RoundNearest) / π,
        2lat / π,
    )
    (true, x, y)
end

"""
    equiprojinv(x, y)

Inverse equirectangular projection. Given a point `(x, y)`
on the plane `[-1, 1] × [-1, 1]`, return a tuple `(Bool, Number, Number)`
where the first Boolean is a flag telling if the point falls
within the projection (`true`) or not (`false`), and the two numbers
are the latitude and longitude in radians.
"""
function equiprojinv(x, y)
    ((-1 ≤ x ≤ 1) && (-1 ≤ y ≤ 1)) || return (false, 0, 0)

    (true, π / 2 * y, π * x)
end

function find_mollweide_theta(ϕ; threshold = 1e-7)
    abs(abs(ϕ) - π/2) < threshold && return ϕ

    θ = ϕ
    while true
        new_θ = θ - (2θ + sin(2θ) - π * sin(ϕ)) / (2 + 2cos(2θ))
        (abs(new_θ - θ) < threshold) && return θ
        θ = new_θ
    end
end

"""
    mollweideproj(lat, lon)

Mollweide projection. Given the latitude `lat` (in radians) and the
longitude (in radians), return a tuple `(Bool, Number, Number)` where the
first Boolean is a flag telling if the point falls within the projection (`true`)
or not (`false`), and the two numbers are the x and y coordinates of the point
on the projection plane (both are in the range [`−1, 1]`).
"""
function mollweideproj(lat, lon)
    θ = find_mollweide_theta(lat)

    (true, -1 / π * lon * cos(θ), sin(θ))
end

"""
    mollweideprojinv(x, y)

Inverse Mollweide projection. Given a point `(x, y)` on the plane,
with `x ∈ [-1, 1]`, `y ∈ [-1, 1]`, return a 3-tuple of type
`(Bool, Number, Number)`. The boolean specifies if `(x, y)` falls within
the map (`true`) or not (`false`), the second and third arguments are
the latitude and longitude in radians.
"""
function mollweideprojinv(x, y)
    # See https://en.wikipedia.org/wiki/Mollweide_projection, we set
    #
    #     R = 1/√2
    #
    # x ∈ [-1, 1], y ∈ [-1, 1]

    x^2 + y^2 ≥ 1 && return (false, 0.0, 0.0)

    sinθ = y
    cosθ = sqrt(1 - sinθ^2)
    θ = asin(sinθ)

    lat = asin((2θ + 2 * sinθ * cosθ) / π)
    long = -2 * π * x / (2cosθ)
    (true, lat, long)

end

"""
    orthoinv(x, y, ϕ1, λ0)

Inverse orthographic projection centered on `(ϕ1, λ0).` Given a
point `(x, y)` on the plane, with `x ∈ [-1, 1]`, `y ∈ [-1, 1]`, return
a 3-tuple of type `(Bool, Number, Number)`. The boolean specifies
if `(x, y)` falls within the map (`true`) or not (`false`), the second
and third arguments are the latitude and longitude in radians.
"""
function orthoinv(x, y, ϕ1, λ0)
    # Assume R = 1/√2. The notation ϕ1, λ0 closely follows
    # the book "Map projections — A working manual" by
    # John P. Snyder (page 145 and ff.)

    R = 1
    ρ = √(x^2 + y^2)
    if ρ > R
        return (false, zero(ϕ1), zero(λ0))
    end

    c = asin(ρ / R)
    sinc, cosc = sin(c), cos(c)
    if cosc < 0
        return (false, zero(ϕ1), zero(λ0))
    end

    if ρ ≈ 0
        return (true, ϕ1, λ0)
    end

    ϕ = asin(cosc * sin(ϕ1) + y * sinc * cos(ϕ1) / ρ)
    if ϕ1 ≈ π / 2
        λ = λ0 + atan(x, -y)
    elseif ϕ1 ≈ -π / 2
        λ = λ0 + atan(x, y)
    else
        λ = λ0 + atan(x * sinc, (ρ * cos(ϕ1) * cosc - y * sin(ϕ1) * sinc))
    end

    (true, ϕ, λ)
end

"""
    ortho2inv(x, y, ϕ1, λ0)

Inverse stereo orthographic projection centered on `(ϕ1, λ0)`. Given
a point `(x, y)` on the plane, with `x ∈ [-1, 1]`, `y ∈ [-1, 1]`, return
a 3-tuple of type `(Bool, Number, Number)`. The boolean specifies
if `(x, y)` falls within the map (`true`) or not (`false`), the second
and third arguments are the latitude and longitude in radians.
"""
function ortho2inv(x, y, ϕ1, λ0)
    x ≤ 0 && return orthoinv(2x + 1, y, ϕ1, λ0)
    orthoinv(2x - 1, y, ϕ1, λ0 + π)
end

"""
    gnominv(x, y, ϕ1, λ0, fov_rad)

Gnomonic projection centered on `(ϕ1, λ0)`, with a field of view
equal to `fov_rad` (in radians).  Given a point `(x, y)` on the plane,
with `x ∈ [-1, 1]`, `y ∈ [-1, 1]`, return a 3-tuple of type `(Bool,
Number, Number)`. The boolean specifies if `(x, y)` falls within
the map (`true`) or not (`false`), the second and third arguments
are the latitude and longitude in radians.
"""
function gnominv(x, y, ϕ1, λ0, ψ0, fov_rad)
    # We fix a Earth radius such that the field of view of the
    # projection is the one expected. Note that the formula
    # diverges if `fov_rad` is 90° (as expected).
    R = 1.0 / tan(fov_rad)

    # We use basic geometry to cast a 3D ray from the center of
    # the sphere to the tangent plane, and then we compute the
    # intersection between the ray and the sphere.
    gamma = 1 / sqrt(1 + (x^2 + y^2) / R^2)
    vecx, vecy, vecz = gamma * R, -gamma * x, -gamma * y

    # We implement the rotations in the following order:
    # 1. Rotation by -ψ0 (orientation) around the x axis
    # 2. Rotation by λ0 (longitude) around the z axis
    # 3. Rotation by ϕ1 (latitude) around the y axis

    sin_ϕ1, cos_ϕ1 = sincos(ϕ1)
    sin_λ0, cos_λ0 = sincos(λ0)
    sin_ψ0, cos_ψ0 = sincos(-ψ0) # Change the sign to match Healpy conventions

    vecxrot = (
        cos_ϕ1 * cos_λ0 * vecx +
        (-cos_ψ0 * cos_ϕ1 * sin_λ0 + sin_ψ0 * sin_ϕ1) * vecy +
        (-sin_ψ0 * cos_ϕ1 * sin_λ0 + sin_ϕ1 * cos_ψ0) * vecz
    )
    vecyrot = (
        (sin_λ0) * vecx +
        (cos_ψ0 * cos_λ0 + sin_ψ0 * sin_ϕ1 * sin_λ0) * vecy +
        (-sin_ψ0 * cos_λ0) * vecz
    )
    veczrot = (
        (-sin_ϕ1 * cos_λ0) * vecx +
        (sin_ψ0 * sin_ϕ1 * sin_λ0 + cos_ϕ1 * sin_ψ0) * vecy +
        (-sin_ψ0 * sin_ϕ1 * sin_λ0 + cos_ϕ1 * cos_ψ0) * vecz
    )

    theta, phi = vec2ang(vecxrot, vecyrot, veczrot)
    (true, colat2lat(theta), phi)
end

################################################################################

"""
    equirectangular(m::HealpixMap{T,O,AA}; kwargs...) where {T <: Number, O, AA}

High-level wrapper around `project` for equirectangular projections.
"""
function equirectangular(m::HealpixMap{T,O,AA}, projparams = Dict()) where {T<:Number,O,AA}
    width = get(projparams, :width, 720)
    height = get(projparams, :height, width)
    project(equiprojinv, m, width, height, projparams)
end

"""
    mollweide(m::HealpixMap{T, O, AA}, projparams = Dict()) where {T <: Number, O, AA}

High-level wrapper around `project` for Mollweide projections.

The following parameters can be set in the `projparams` dictionary:

- `width`: width of the image, in pixels (default: 720 pixels)
- `height`: height of the image, in pixels; if not specified, it will be assumed
  to be equal to `width`
"""
function mollweide(m::HealpixMap{T,O,AA}, projparams = Dict()) where {T<:Number,O,AA}
    width = get(projparams, :width, 720)
    height = get(projparams, :height, width ÷ 2)
    project(mollweideprojinv, m, width, height, projparams)
end

"""
    orthographic(m::HealpixMap{T,O}, projparams = Dict()) where {T <: Number,O <: Order}

High-level wrapper around `project` for orthographic projections.

The following parameters can be set in the `projparams` dictionary:

- `width`: width of the image, in pixels (default: 720 pixels)
- `height`: height of the image, in pixels; if not specified, it will be assumed
  to be equal to `width`
- `center`: position of the pixel in the middle of the left globe (*latitude* and
  longitude).
"""
function orthographic(m::HealpixMap{T,O,AA}, projparams = Dict()) where {T<:Number,O,AA}
    width = get(projparams, :width, 720)
    height = get(projparams, :height, width)
    ϕ0, λ0 = get(projparams, :center, (0, 0))
    project(m, width, height, projparams) do x, y
        orthoinv(x, y, ϕ0, λ0)
    end
end

"""
    orthographic2(m::HealpixMap{T, O, AA}, projparams = Dict()) where {T <: Number, O, AA}

High-level wrapper around `project` for stereo orthographic projections.

The following parameters can be set in the `projparams` dictionary:

- `width`: width of the image, in pixels (default: 720 pixels)
- `height`: height of the image, in pixels; if not specified, it will be assumed
  to be equal to `width`
- `center`: position of the pixel in the middle of the left globe (*latitude* and
  longitude). Default is (0, 0).
"""
function orthographic2(m::HealpixMap{T,O,AA}, projparams = Dict()) where {T<:Number,O,AA}
    width = get(projparams, :width, 720)
    height = get(projparams, :height, width ÷ 2)
    ϕ0, λ0 = get(projparams, :center, (0, 0))
    project(m, width, height, projparams) do x, y
        ortho2inv(x, y, ϕ0, λ0)
    end
end

"""
    gnomonic(m::HealpixMap{T, O, AA}, projparams = Dict()) where {T <: Number, O, AA}

High-level wrapper around `project` for gnomonic projections.

The following parameters can be set in the `projparams` dictionary:

- `width`: width of the image, in pixels (default: 720 pixels)
- `height`: height of the image, in pixels; if not specified, it will be assumed
  to be equal to `width`
- `center`: position and orientation of the pixel in the middle. It is a 3-element
  tuple containing:

  1. The *latitude* of the pixel, in radians
  2. The longitude of the pixel, in radians
  3. The rotation to be applied to the image, in radians

- `fov_rad`: size of the image along the x and y axes, in radians (default: 15°)

# Example

````julia
plot(m, gnomonic, Dict(:fov_rad = deg2rad(1.5), :center = (0, 0, deg2rad(45))))
````
"""
function gnomonic(m::HealpixMap{T,O,AA}, projparams = Dict()) where {T<:Number,O,AA}
    width = get(projparams, :width, 720)
    height = get(projparams, :height, width)
    ϕ0, λ0, ψ0 = get(projparams, :center, (0, 0, 0))
    fov_rad = get(projparams, :fov_rad, deg2rad(15))
    project(m, width, height, projparams) do x, y
        gnominv(x, y, ϕ0, λ0, ψ0, fov_rad)
    end
end

################################################################################

RecipesBase.@recipe function plot(
    m::HealpixMap{T,O,AA},
    projection = mollweide,
    projparams = Dict(),
) where {T<:Number,O,AA}

    img, mask, anymasked = projection(m, projparams)

    if anymasked
        RecipesBase.@series begin
            seriestype --> :shape
            primary --> false
            c --> :grey
            line --> :grey

            width, height = size(mask)
            xm = Float64[]
            ym = Float64[]
            for y = 1:height
                curx = 1
                # Instead of drawing each single masked pixel as a
                # square, squash together long runs of masked pixels
                # into one rectangle
                while curx ≤ width
                    if mask[curx, y]
                        startx = curx
                        while curx ≤ width && mask[curx, y]
                            curx += 1
                        end

                        append!(
                            xm,
                            [startx - 0.5, startx - 0.5, curx + 0.5, curx + 0.5, NaN],
                        )
                        append!(ym, [y - 0.5, y + 0.5, y + 0.5, y - 0.5, NaN])
                    else
                        curx += 1
                    end
                end
            end

            ym, xm
        end
    end

    seriestype --> :heatmap
    aspect_ratio --> 1
    colorbar --> :bottom
    framestyle --> :none

    img
end

@doc raw"""
    plot(m::HealpixMap{T, O, AA}, projection = mollweide, projparams = Dict())

Draw a representation of the map, using some specific projection. The
parameter `projection` must be a function returning the
bitmap. Possible values for `projection` are the following:

- `equirectangular`
- `mollweide`
- `orthographic`
- `orthographic2`
- `gnomonic`

You can define your own projections.

The dictionary `projparams` allows to hack a number of parameters used
in the projection.

# References

See also [`equirectangular`](@ref), [`mollweide`](@ref),
[`orthographic`](@ref), [`orthographic2`](@ref), and
[`gnomonic`](@ref).
"""
plot
