# Map projections

export lat2colat, project, equiprojinv, mollweideprojinv, equirectangular, mollweide

using FileIO
import Cairo
import Images
import ColorSchemes

function drawmapbmp(invprojfn, m::Map{T,O}, bmpwidth, bmpheight;
                    background_color=Images.RGBA{Images.N0f8}(1.0, 1.0, 1.0, 0.0),
                    unseen_color=Images.RGBA{Images.N0f8}(0.5, 0.5, 0.5, 1.0),
                    cs=ColorSchemes.viridis,
                    minval=missing,
                    maxval=missing,
                    unseen=missing) where {T <: Number, O <: Healpix.Order}

    if ismissing(minval)
        minval = convert(Float64, minimum(m.pixels[isfinite.(m.pixels)]))
    end
    if ismissing(maxval)
        maxval = convert(Float64, maximum(m.pixels[isfinite.(m.pixels)]))
    end

    if maxval ≈ minval
        maxval = minval + 1.0
    end

    img = Array{Images.RGBA{Images.N0f8}}(undef, bmpheight, bmpwidth)

    for j in 1:bmpheight
        y = 2 * (j - 1) / (bmpheight - 1) - 1
        for i in 1:bmpwidth
            x = 2 * (i - 1) / (bmpwidth - 1) - 1
            # Flip the sign of the y coordinate because of Cairo's coordinate convention
            visible, lat, long = invprojfn(x, -y)
            if visible
                value = m.pixels[Healpix.ang2pix(m, lat2colat(lat), long)]
                if ismissing(value) || isnan(value) || (
                    !ismissing(unseen) && unseen == value)
                    color = unseen_color
                else
                    color = get(cs, (value - minval) / (maxval - minval))
                end

                img[j, i] = color
            else
                img[j, i] = background_color
            end
        end
    end

    img, minval, maxval
end

"""
   lat2colat(x)

Convert latitude into colatitude. Both `x` and the result are expressed in radians.
"""
lat2colat(x) = π / 2 - x

"""
    project(m::Map{T, O}; kwargs...) where {T <: Number, O <: Order}

Return a 2D bitmap (array) containing a cartographic projection of the map and a
2D bitmap containing a boolean mask. The size of the bitmap is specified by
`figsize`, which must be a 2-tuple. The function `projfn` must be a function
which accepts as input two parameters `x` and `y` (numbers between -1 and 1).

The following keywords can be used in the call:

- `figsize`: 2-tuple specifying the (height, width) of the bitmap in pixels
- `center`: 2-tuple specifying the location (colatitude, longitude) of the sky
  point that is to be placed in the middle of the image (in radians)
- `numfmt`: Number->String function to be used to convert the extrema of the color
  bar into a textual representation. The default is `x -> @sprintf("%g", x)`.
- `returnmask`: Boolean; if true, the function returns a 2-tuple containing
  the image bitmap and a mask bitmap, which is set to true if the pixel falls
  within the carthographic projection, false otherwise. If `returnmask` is
  false, only the image bitmap is returned.
- `show`: Boolean. If true (the default), the map will be displayed. It has no
  effect if `returnmask` is true.  
"""
function project(invprojfn, drawborderfn, m::Map{T,O}, bmpwidth, bmpheight;
                 kwargs...) where {T <: Number, O <: Order}

    args = Dict(kwargs)
    figsize = get(args, :figsize, (400, 800))
    center = get(args, :center, (0, 0))
    returnmask = get(args, :returnmask, false)
    show = get(args, :show, true)
    numfmt = get(args, :numfmt, x -> @sprintf("%g", x))
    colorscheme = get(args, :cs, ColorSchemes.temperaturemap)
    cbarmargin = get(args, :cbarmargin, 10)
    cbarheight = get(args, :cbarheight, 50)
    cbarlblheight = get(args, :cbarlblheight, 20)
    cbarsteps = get(args, :cbarsteps, 80)
    minval = get(args, :minval, missing)
    maxval = get(args, :maxval, missing)

    img, minval, maxval = drawmapbmp(invprojfn, m, bmpwidth, bmpheight,
                                     minval=minval, maxval=maxval,
                                     cs=colorscheme)

    io = IOBuffer()
    Images.save(Stream(format"PNG", io), img)
    seek(io, 0);

    image = Cairo.read_from_png(io)
    imgw = image.width; 
    imgh = image.height;

    width, height = imgw, imgh + cbarheight + cbarlblheight
    c = Cairo.CairoRGBSurface(width, height);
    cr = Cairo.CairoContext(c);

    Cairo.set_source_rgb(cr, 1, 1, 1)
    Cairo.rectangle(cr, 0, 0, width, height)
    Cairo.fill(cr)
    
    Cairo.save(cr);

    Cairo.translate(cr, (width - imgw) / 2, 0.0);
    Cairo.set_source_surface(cr, image, 0, 0);
    Cairo.paint(cr);

    Cairo.set_source_rgb(cr, 0.0, 0.0, 0.0)
    drawborderfn(cr, imgw, imgh)

    Cairo.restore(cr);

    Cairo.save(cr)

    stepsize = (width - 2cbarmargin) / cbarsteps
    Cairo.translate(cr, 0 + cbarmargin, imgh + cbarmargin)
    Cairo.scale(cr, width - 2cbarmargin, cbarheight - 2cbarmargin)

    Cairo.save(cr)
    Cairo.rectangle(cr, 0.0, 0.0, 1.0, 1.0)
    Cairo.clip(cr)
    for i = range(0.0, stop=(1.0 - 1.0 / cbarsteps), length=cbarsteps)
        color = get(colorscheme, i)
        Cairo.set_source_rgb(cr, color.r, color.g, color.b)
        Cairo.rectangle(cr, i, 0, 1.1 / cbarsteps, 1.0)
        Cairo.fill(cr)
    end
    Cairo.restore(cr)

    Cairo.set_source_rgb(cr, 0, 0, 0)
    Cairo.rectangle(cr, 0, 0, 1, 1)
    Cairo.stroke(cr)

    Cairo.restore(cr)

    Cairo.set_source_rgb(cr, 0, 0, 0)
    Cairo.set_font_size(cr, cbarlblheight)

    text_baseline = imgh + cbarheight + cbarlblheight - cbarmargin
    Cairo.move_to(cr, cbarmargin, text_baseline)
    Cairo.show_text(cr, numfmt(minval))

    hilabel = numfmt(maxval)
    Cairo.move_to(cr,
        width - cbarmargin - Cairo.textwidth(cr, hilabel),
        text_baseline)
    Cairo.show_text(cr, hilabel)

    c
end

function equiprojborder(cr, w, h)
    Cairo.rectangle(cr, 0, 0, w, h)
    Cairo.stroke(cr)
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

function mollweideborder(cr, w, h)
    Cairo.save(cr)

    Cairo.translate(cr, w / 2, h / 2)
    Cairo.scale(cr, 1.0, h / w)
    Cairo.translate(cr, -w / 2, -h / 2)
    Cairo.arc(cr, w / 2, h / 2, w / 2, 0, 2π)

    Cairo.restore(cr)

    Cairo.stroke(cr)
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

    if x^2 + y^2 ≥ 1
        return (false, 0, 0)
    end

    sinθ = y
    cosθ = sqrt(1 - sinθ^2)
    θ = asin(sinθ)

    lat = asin((2θ + 2 * sinθ * cosθ) / π)
    long = -2 * π * x / (2cosθ)
    (true, lat, long)

end

function orthoborder(cr, w, h)
    Cairo.arc(cr, w / 2, h / 2, w / 2, 0, 2π)
    Cairo.stroke(cr)
end

"""
    function orthoinv(x, y, ϕ1, λ0)

Inverse orthographic projection centered on (ϕ1, λ0). Given a
point (x, y) on the plane, with x ∈ [-1, 1], y ∈ [-1, 1], return
a 3-tuple of type (Bool, Number, Number). The boolean specifies
if (x, y) falls within the map (true) or not (false), the second
and third arguments are the latitude and longitude in radians.
"""
function orthoinv(x, y, ϕ1, λ0; kwargs...)
    # Assume R = 1/√2. The notation ϕ1, λ0 closely follows
    # the book "Map projections — A working manual" by
    # John P. Snyder (page 145 and ff.)
    
    R = 1
    ρ = √(x^2 + y^2)
    if ρ > R
        return (false, 0, 0)
    end
    
    c = asin(ρ / R)
    sinc, cosc = sin(c), cos(c)
    if cosc < 0
        return (false, 0, 0)
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
    equirectangular(m::Map{T,O}; kwargs...) where {T <: Number, O <: Order}

High-level wrapper around `project` for equirectangular projections.
"""
function equirectangular(m::Map{T,O}; kwargs...) where {T <: Number, O <: Order}
    project(equiprojinv, equiprojborder, m, 720, 720; kwargs...)
end

"""
    mollweide(m::Map{T,O}; kwargs...) where {T <: Number, O <: Order}

High-level wrapper around `project` for Mollweide projections.
"""
function mollweide(m::Map{T,O}; kwargs...) where {T <: Number, O <: Order}
    project(mollweideprojinv, mollweideborder, m, 720, 720 ÷ 2; kwargs...)
end

"""
    orthographic(m::Map{T,O}, ϕ0, λ0; kwargs...) where {T <: Number, O <: Order}

High-level wrapper around `project` for orthographic projections centered around the point (ϕ0, λ0).
"""
function orthographic(m::Map{T,O}, ϕ0, λ0; kwargs...) where {T <: Number, O <: Order}
    project(orthoborder, m, 540, 540; kwargs...) do x, y
        orthoinv(x, y, ϕ0, λ0)
    end
end
