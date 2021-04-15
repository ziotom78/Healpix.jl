################################################################################

ang2pix(map::Map{T,RingOrder,AA}, theta, phi) where {T,AA} =
    ang2pixRing(map.resolution, theta, phi)
ang2pix(map::Map{T,NestedOrder,AA}, theta, phi) where {T,AA} =
    ang2pixNest(map.resolution, theta, phi)
ang2pix(map::PolarizedMap{T,RingOrder,AA}, theta, phi) where {T,AA} =
    ang2pixRing(map.i.resolution, theta, phi)
ang2pix(map::PolarizedMap{T,NestedOrder,AA}, theta, phi) where {T,AA} =
    ang2pixNest(map.i.resolution, theta, phi)

@doc raw"""
    ang2pix{T, O, AA}(map::Map{T, O}, theta, phi)
    ang2pix{T, O, AA}(map::PolarizedMap{T, O}, theta, phi)

Convert the direction specified by the colatitude `theta` (∈ [0, π])
and the longitude `phi` (∈ [0, 2π]) into the index of the pixel in the
Healpix map `map`.
"""
ang2pix

################################################################################

pix2ang(map::Map{T,RingOrder,AA}, ipix) where {T,AA} = pix2angRing(map.resolution, ipix)
pix2ang(map::Map{T,NestedOrder,AA}, ipix) where {T,AA} = pix2angNest(map.resolution, ipix)
pix2ang(map::PolarizedMap{T,RingOrder,AA}, ipix) where {T,AA} =
    pix2angRing(map.i.resolution, ipix)
pix2ang(map::PolarizedMap{T,NestedOrder,AA}, ipix) where {T,AA} =
    pix2angNest(map.i.resolution, ipix)

@doc raw"""
    pix2ang{T, O <: Order}(map::Map{T, O}, ipix) -> (Float64, Float64)
    pix2ang{T, O <: Order}(map::PolarizedMap{T, O}, ipix) -> (Float64, Float64)

Return the pair (`theta`, `phi`), where `theta` is the colatitude and
`phi` the longitude of the direction of the pixel center with index
`ipix`.
"""
pix2ang

################################################################################
# Interpolation

function interpolate(m::Map{T,RingOrder,AA}, θ, ϕ, pixbuf, weightbuf) where {T,AA}
    getinterpolRing(m.resolution, θ, ϕ, pixbuf, weightbuf)

    result = zero(weightbuf[1])
    for i = 1:4
        result += m[pixbuf[i]] * weightbuf[i]
    end

    result
end

function interpolate(m::Map{T,RingOrder,AA}, θ, ϕ) where {T,AA}
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)

    interpolate(m, θ, ϕ, pixbuf, weightbuf)
end

"""
    interpolate(m::Map{T, RingOrder, AA}, θ, ϕ) -> Value
    interpolate(m::Map{T, RingOrder, AA}, θ, ϕ, pixbuf, weightbuf) -> Value

Return an interpolated value of the map along the specified direction.

When provided, the parameters `pixbuf` and `weightbuf` must be
4-element arrays of integer and floating-point values,
respectively. They can be reused across multiple calls of
`interpolate!`, to save heap allocations:

```
pixbuf = Array{Int}(undef, 4)
weightbuf = Array{Float64}(undef, 4)

m = Map{Float64, RingOrder}(1)
for (θ, ϕ) in [(0., 0.), (π/2, π/2)]
    println(interpolate!(m, θ, ϕ, pixbuf, weightbuf))
end
```
"""
interpolate


@doc raw"""
    nest2ring(m_nest::Map{T, NestedOrder, AA}) where {T, AA}

Convert a map from nested to ring order. This version allocates a new array of the same
array type as the input.

# Arguments:
- `m_nest::Map{T, NestedOrder, AA}`: map of type `NestedOrder`

# Returns: 
- `Map{T, RingOrder, AA}`: the input map converted to `RingOrder`

# Examples
```julia-repl
julia> m_nest = Map{Float64,NestedOrder}(rand(nside2npix(64)));

julia> nest2ring(m_nest)
49152-element Map{Float64, RingOrder, Vector{Float64}}:
 0.4703834205807309
 ⋮
 0.3945848051663148
```
"""
function nest2ring(m_nest::Map{T, NestedOrder, AA}) where {T, AA}
    m_ring = Map{T, RingOrder, AA}(AA(undef, size(m_nest.pixels)))
    nest2ring!(m_ring, m_nest)
    return m_ring
end


@doc raw"""
    nest2ring!(m_ring_dst::Map{T, RingOrder, AAR}, 
               m_nest_src::Map{T, NestedOrder, AAN}) where {T, AAN, AAR}

Convert a map from nested to ring order. This version takes a nested map in the 
second argument and writes it to the nested map provided in the first argument,
following the standard Julia `func!(dst, src)` convention.

# Arguments:
- `m_ring_dst::Map{T, NestedOrder, AA}`: map of type `NestedOrder`
- `m_nest_src::Map{T, NestedOrder, AAN}`: map of type `RingOrder`

# Returns: 
- `Map{T, RingOrder, AA}`: the input map converted to `RingOrder`

# Examples
```julia-repl
julia> m_nest = Map{Float64,NestedOrder}(rand(nside2npix(64)));

julia> m_ring = Map{Float64,RingOrder}(64);

julia> nest2ring!(m_ring, m_nest)
49152-element Map{Float64, RingOrder, Vector{Float64}}:
 0.33681791815569895
 ⋮
 0.9092457003948482
```
"""
function nest2ring!(m_ring_dst::Map{T, RingOrder, AAR}, m_nest_src::Map{T, NestedOrder, AAN}) where {T, AAN, AAR}
    res_nest = m_nest_src.resolution
    for i_nest in eachindex(m_nest_src.pixels)
        i_ring = nest2ring(res_nest, i_nest)
        m_ring_dst.pixels[i_ring] = m_nest_src.pixels[i_nest]
    end
    return m_ring_dst
end


@doc raw"""
    ring2nest(m_ring::Map{T, RingOrder, AA}) where {T, AA}

Convert a map from ring to nested order. This version allocates a new array of the same
array type as the input.

# Arguments:
- `m_ring::Map{T, RingOrder, AA}`: map of type `RingOrder`

# Returns: 
- `Map{T, NestedOrder, AA}`: the input map converted to `NestedOrder`

# Examples
```julia-repl
julia> m_ring = Map{Float64,RingOrder}(rand(nside2npix(64)));

julia> ring2nest(m_ring)
49152-element Map{Float64, NestedOrder, Vector{Float64}}:
 0.0673134062168923
 ⋮
 0.703460503535335
```
"""
function ring2nest(m_ring::Map{T, RingOrder, AA}) where {T, AA}
    m_nest = Map{T, NestedOrder, AA}(AA(undef, size(m_ring.pixels)))
    ring2nest!(m_nest, m_ring)
    return m_nest
end


@doc raw"""
    ring2nest!(m_nest_dst::Map{T, NestedOrder, AAN}, 
               m_ring_src::Map{T, RingOrder, AAR}) where {T, AAR, AAN}

Convert a map from ring to nested order. This version takes a nested map in the 
second argument and writes it to the nested map provided in the first argument,
following the standard Julia `func!(dst, src)` convention.

# Arguments:
- `m_nest_dst::Map{T, NestedOrder, AAN}`: map of type `RingOrder`
- `m_ring_src::Map{T, RingOrder, AA}`: map of type `RingOrder`

# Returns: 
- `Map{T, NestedOrder, AA}`: the input map converted to `NestedOrder`

# Examples
```julia-repl
julia> m_ring = Map{Float64,RingOrder}(rand(nside2npix(64)));

julia> m_nest = Map{Float64,RingOrder}(64);

julia> ring2nest!(m_nest, m_ring)
49152-element Map{Float64, NestedOrder, Vector{Float64}}:
 0.0673134062168923
 ⋮
 0.703460503535335
```
"""
function ring2nest!(m_nest_dst::Map{T, NestedOrder, AAN}, m_ring_src::Map{T, RingOrder, AAR}) where {T, AAR, AAN}
    res_ring = m_ring_src.resolution
    for i_ring in eachindex(m_ring_src.pixels)
        i_nest = ring2nest(res_ring, i_ring)
        m_nest_dst.pixels[i_nest] = m_ring_src.pixels[i_ring]
    end
end



"""
    udgrade(input_map::Map{T,O,AA}, output_nside) where {T,O,AA} -> Map{T,O,AA}

Upgrades or downgrades a map to a target nside. Always makes a copy. This is very fast 
for nested orderings, but slow for ring because one needs to transform to nested ordering 
first.

# Arguments:
- `input_map::Map{T,O,AA}`: the map to upgrade/downgrade
- `output_nside`: desired nside

# Returns: 
- `Map{T,O,AA}`: upgraded/downgraded map

# Examples
```julia-repl
julia> A = Map{Float64, NestedOrder}(ones(nside2npix(4)))
192-element Map{Float64, RingOrder, Vector{Float64}}:
 1.0
 ⋮
 1.0

julia> Healpix.udgrade(A, 2)
48-element Map{Float64, NestedOrder, Vector{Float64}}:
 1.0
 ⋮
 1.0
```
"""
function udgrade(map_in::Map{T,O,AA}, nside_out; 
                 threshold=abs(1e-6UNSEEN), pess=false) where {T,O<:NestedOrder,AA}
    nside_in = map_in.resolution.nside
    npix_out = nside2npix(nside_out)
    map_out = Map{T,O,AA}(nside_out)

    if nside_in == nside_out
        map_out.pixels .= input_map.pixels
    elseif nside_out < nside_in  # degrade. loop over input and add them up 
        npratio = nside2npix(nside_in) ÷ nside2npix(nside_out)
        for id = 0:(npix_out-1)
            nobs = 0
            total = 0.0
            for ip = 0:(npratio-1)
                value = map_in[id*npratio + ip + 1]
                if (abs(value - UNSEEN) > threshold) 
                    nobs  = nobs  + 1
                    total = total + value
                end
            end
            map_out[id+1] = UNSEEN
            if pess
                if nobs == npratio
                    map_out[id+1] = total/nobs
                end
            else
                if nobs > 0
                    map_out[id+1] = total/nobs
                end
            end
        end
    else  # we are upgrading. loop over output map pixels, and set them to parent pixel
        npratio = nside2npix(nside_out) ÷ nside2npix(nside_in)
        for iu = 0:(npix_out - 1)
            ip = iu ÷ npratio
            map_out[iu+1] = map_in[ip+1]
        end
    end

    return map_out
end
# just convert to nest and udgrade if we get a ring
function udgrade(map_in::Map{T,O,AA}, nside_out; 
                 threshold=abs(1e-6UNSEEN), pess=true) where {T,O<:RingOrder,AA}
    map_nest = ring2nest(map_in)
    map_out = udgrade(map_nest, nside_out; threshold=threshold, pess=pess)
    return nest2ring(map_out)
end
