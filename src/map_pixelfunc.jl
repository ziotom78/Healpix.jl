################################################################################

ang2pix(map::HealpixMap{T,RingOrder,AA}, theta, phi) where {T,AA} =
    ang2pixRing(map.resolution, theta, phi)
ang2pix(map::HealpixMap{T,NestedOrder,AA}, theta, phi) where {T,AA} =
    ang2pixNest(map.resolution, theta, phi)
ang2pix(map::PolarizedHealpixMap{T,RingOrder,AA}, theta, phi) where {T,AA} =
    ang2pixRing(map.i.resolution, theta, phi)
ang2pix(map::PolarizedHealpixMap{T,NestedOrder,AA}, theta, phi) where {T,AA} =
    ang2pixNest(map.i.resolution, theta, phi)

@doc raw"""
    ang2pix{T, O, AA}(map::HealpixMap{T, O}, theta, phi)
    ang2pix{T, O, AA}(map::PolarizedHealpixMap{T, O}, theta, phi)

Convert the direction specified by the colatitude `theta` (∈ [0, π])
and the longitude `phi` (∈ [0, 2π]) into the index of the pixel in the
Healpix map `map`.
"""
ang2pix

################################################################################

pix2ang(map::HealpixMap{T,RingOrder,AA}, ipix) where {T,AA} = pix2angRing(map.resolution, ipix)
pix2ang(map::HealpixMap{T,NestedOrder,AA}, ipix) where {T,AA} = pix2angNest(map.resolution, ipix)
pix2ang(map::PolarizedHealpixMap{T,RingOrder,AA}, ipix) where {T,AA} =
    pix2angRing(map.i.resolution, ipix)
pix2ang(map::PolarizedHealpixMap{T,NestedOrder,AA}, ipix) where {T,AA} =
    pix2angNest(map.i.resolution, ipix)

@doc raw"""
    pix2ang{T, O <: Order}(map::HealpixMap{T, O}, ipix) -> (Float64, Float64)
    pix2ang{T, O <: Order}(map::PolarizedHealpixMap{T, O}, ipix) -> (Float64, Float64)

Return the pair (`theta`, `phi`), where `theta` is the colatitude and
`phi` the longitude of the direction of the pixel center with index
`ipix`.
"""
pix2ang

################################################################################
# Interpolation

function interpolate(m::HealpixMap{T,RingOrder,AA}, θ, ϕ, pixbuf, weightbuf) where {T,AA}
    getinterpolRing(m.resolution, θ, ϕ, pixbuf, weightbuf)

    result = zero(weightbuf[1])
    for i = 1:4
        result += m[pixbuf[i]] * weightbuf[i]
    end

    result
end

function interpolate(m::HealpixMap{T,RingOrder,AA}, θ, ϕ) where {T,AA}
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)

    interpolate(m, θ, ϕ, pixbuf, weightbuf)
end

"""
    interpolate(m::HealpixMap{T, RingOrder, AA}, θ, ϕ) -> Value
    interpolate(m::HealpixMap{T, RingOrder, AA}, θ, ϕ, pixbuf, weightbuf) -> Value

Return an interpolated value of the map along the specified direction.

When provided, the parameters `pixbuf` and `weightbuf` must be
4-element arrays of integer and floating-point values,
respectively. They can be reused across multiple calls of
`interpolate`, to save heap allocations, and they do not need to be
initialized, as they are used internally:

```
pixbuf = Array{Int}(undef, 4)
weightbuf = Array{Float64}(undef, 4)

m = HealpixMap{Float64, RingOrder}(1)
for (θ, ϕ) in [(0., 0.), (π/2, π/2)]
    println(interpolate(m, θ, ϕ, pixbuf, weightbuf))
end
```

Passing `pixbuf` and `weightbuf` saves some time, as this simple benchmark
shows:

```
julia> @benchmark interpolate(m, rand(), rand(), pixbuf, weightbuf)
BenchmarkTools.Trial: 
  memory estimate:  618 bytes
  allocs estimate:  9
  --------------
  minimum time:     283.184 ns (0.00% GC)
  median time:      296.140 ns (0.00% GC)
  mean time:        348.132 ns (9.55% GC)
  maximum time:     10.627 μs (95.88% GC)
  --------------
  samples:          10000
  evals/sample:     282

julia> @benchmark interpolate(m, rand(), rand())
BenchmarkTools.Trial: 
  memory estimate:  837 bytes
  allocs estimate:  11
  --------------
  minimum time:     329.825 ns (0.00% GC)
  median time:      345.504 ns (0.00% GC)
  mean time:        417.004 ns (11.04% GC)
  maximum time:     13.733 μs (96.21% GC)
  --------------
  samples:          10000
  evals/sample:     223
```

"""
interpolate


@doc raw"""
    nest2ring(m_nest::HealpixMap{T, NestedOrder, AA}) where {T, AA}

Convert a map from nested to ring order. This version allocates a new array of the same
array type as the input.

# Arguments:
- `m_nest::HealpixMap{T, NestedOrder, AA}`: map of type `NestedOrder`

# Returns: 
- `HealpixMap{T, RingOrder, AA}`: the input map converted to `RingOrder`

# Examples
```julia-repl
julia> m_nest = HealpixMap{Float64,NestedOrder}(rand(nside2npix(64)));

julia> nest2ring(m_nest)
49152-element HealpixMap{Float64, RingOrder, Vector{Float64}}:
 0.4703834205807309
 ⋮
 0.3945848051663148
```
"""
function nest2ring(m_nest::HealpixMap{T, NestedOrder, AA}) where {T, AA}
    m_ring = HealpixMap{T, RingOrder, AA}(AA(undef, size(m_nest.pixels)))
    nest2ring!(m_ring, m_nest)
    return m_ring
end


@doc raw"""
    nest2ring!(m_ring_dst::HealpixMap{T, RingOrder, AAR}, 
               m_nest_src::HealpixMap{T, NestedOrder, AAN}) where {T, AAN, AAR}

Convert a map from nested to ring order. This version takes a nested map in the 
second argument and writes it to the nested map provided in the first argument,
following the standard Julia `func!(dst, src)` convention.

# Arguments:
- `m_ring_dst::HealpixMap{T, NestedOrder, AA}`: map of type `NestedOrder`
- `m_nest_src::HealpixMap{T, NestedOrder, AAN}`: map of type `RingOrder`

# Returns: 
- `HealpixMap{T, RingOrder, AA}`: the input map converted to `RingOrder`

# Examples
```julia-repl
julia> m_nest = HealpixMap{Float64,NestedOrder}(rand(nside2npix(64)));

julia> m_ring = HealpixMap{Float64,RingOrder}(64);

julia> nest2ring!(m_ring, m_nest)
49152-element HealpixMap{Float64, RingOrder, Vector{Float64}}:
 0.33681791815569895
 ⋮
 0.9092457003948482
```
"""
function nest2ring!(m_ring_dst::HealpixMap{T, RingOrder, AAR}, m_nest_src::HealpixMap{T, NestedOrder, AAN}) where {T, AAN, AAR}
    res_nest = m_nest_src.resolution
    for i_nest in eachindex(m_nest_src.pixels)
        i_ring = nest2ring(res_nest, i_nest)
        m_ring_dst.pixels[i_ring] = m_nest_src.pixels[i_nest]
    end
    return m_ring_dst
end


@doc raw"""
    ring2nest(m_ring::HealpixMap{T, RingOrder, AA}) where {T, AA}

Convert a map from ring to nested order. This version allocates a new array of the same
array type as the input.

# Arguments:
- `m_ring::HealpixMap{T, RingOrder, AA}`: map of type `RingOrder`

# Returns: 
- `HealpixMap{T, NestedOrder, AA}`: the input map converted to `NestedOrder`

# Examples
```julia-repl
julia> m_ring = HealpixMap{Float64,RingOrder}(rand(nside2npix(64)));

julia> ring2nest(m_ring)
49152-element HealpixMap{Float64, NestedOrder, Vector{Float64}}:
 0.0673134062168923
 ⋮
 0.703460503535335
```
"""
function ring2nest(m_ring::HealpixMap{T, RingOrder, AA}) where {T, AA}
    m_nest = HealpixMap{T, NestedOrder, AA}(AA(undef, size(m_ring.pixels)))
    ring2nest!(m_nest, m_ring)
    return m_nest
end


@doc raw"""
    ring2nest!(m_nest_dst::HealpixMap{T, NestedOrder, AAN}, 
               m_ring_src::HealpixMap{T, RingOrder, AAR}) where {T, AAR, AAN}

Convert a map from ring to nested order. This version takes a nested map in the 
second argument and writes it to the nested map provided in the first argument,
following the standard Julia `func!(dst, src)` convention.

# Arguments:
- `m_nest_dst::HealpixMap{T, NestedOrder, AAN}`: map of type `RingOrder`
- `m_ring_src::HealpixMap{T, RingOrder, AA}`: map of type `RingOrder`

# Returns: 
- `HealpixMap{T, NestedOrder, AA}`: the input map converted to `NestedOrder`

# Examples
```julia-repl
julia> m_ring = HealpixMap{Float64,RingOrder}(rand(nside2npix(64)));

julia> m_nest = HealpixMap{Float64,RingOrder}(64);

julia> ring2nest!(m_nest, m_ring)
49152-element HealpixMap{Float64, NestedOrder, Vector{Float64}}:
 0.0673134062168923
 ⋮
 0.703460503535335
```
"""
function ring2nest!(m_nest_dst::HealpixMap{T, NestedOrder, AAN}, m_ring_src::HealpixMap{T, RingOrder, AAR}) where {T, AAR, AAN}
    res_ring = m_ring_src.resolution
    for i_ring in eachindex(m_ring_src.pixels)
        i_nest = ring2nest(res_ring, i_ring)
        m_nest_dst.pixels[i_nest] = m_ring_src.pixels[i_ring]
    end
end



"""
    udgrade(input_map::HealpixMap{T,O,AA}, output_nside; kw...) where {T,O,AA} -> HealpixMap{T,O,AA}

Upgrades or downgrades a map to a target nside. Always makes a copy. This is very fast 
for nested orderings, but slow for ring because one needs to transform to nested ordering 
first.

# Arguments:
- `input_map::HealpixMap{T,O,AA}`: the map to upgrade/downgrade
- `output_nside`: desired nside

# Keywords:
- `threshold=abs(1e-6UNSEEN)`: absolute tolerance for identifying a bad pixel vs UNSEEN
- `pess=false`: if false, estimate pixels from remaining good pixels when downgrading. 
    if true, the entire downgraded pixel is set to UNSEEN.

# Returns: 
- `HealpixMap{T,O,AA}`: upgraded/downgraded map in the same ordering as the input

# Examples
```julia-repl
julia> A = HealpixMap{Float64, NestedOrder}(ones(nside2npix(4)))
192-element HealpixMap{Float64, RingOrder, Vector{Float64}}:
 1.0
 ⋮
 1.0

julia> Healpix.udgrade(A, 2)
48-element HealpixMap{Float64, NestedOrder, Vector{Float64}}:
 1.0
 ⋮
 1.0
```
"""
function udgrade(map_in::HealpixMap{T,O,AA}, nside_out; 
                 threshold=abs(1e-6UNSEEN), pess=false) where {T,O<:NestedOrder,AA}
    nside_in = map_in.resolution.nside
    npix_out = nside2npix(nside_out)
    map_out = HealpixMap{T,O,AA}(nside_out)

    if nside_in == nside_out
        map_out.pixels .= map_in.pixels
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
function udgrade(map_in::HealpixMap{T,O,AA}, nside_out; 
                 threshold=abs(1e-6UNSEEN), pess=false) where {T,O<:RingOrder,AA}
    map_nest = ring2nest(map_in)
    map_out = udgrade(map_nest, nside_out; threshold=threshold, pess=pess)
    return nest2ring(map_out)
end
