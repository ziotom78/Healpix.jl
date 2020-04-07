# this file contains functions related to spherical harmonic transforms

"""
    iterate_map2alm!(
        maps::Array{Array{Float64,1},1}, alms::Array{Array{Complex{Float64},1},1}, 
        geom_info::Libsharp.GeomInfo, alm_info::Libsharp.AlmInfo, 
        niter::Integer, spin::Integer)

This is an internal function for implementing iterative map2alm, used when the parameter 
`niter` is greater than 0. It synthesizes a map from the alms, subtracts it from the map to 
form a residual map, and then adds the harmonic coefficients of the residual map to the 
alms. It repeats this `niter` times. It performs this in-place on arrays of Float64 and
Complex{Float64}.

# Arguments
- `maps::Array{Array{Float64,1},1}`: an array where each element is a 1D Healpix map array.
- `alms::Array{Array{Complex{Float64},1},1}`: an array where each element is a array of alms
- `geom_info::Libsharp.GeomInfo`: contains information about the pixelization
- `alm_info::Libsharp.AlmInfo`: contains information about the SHT coefficient ordering
- `niter::Integer`: number of iterations to perform
- `spin::Integer`: spin of the field, 0 or 2
"""
function iterate_map2alm!(
        maps::Array{Array{Float64,1},1}, alms::Array{Array{Complex{Float64},1},1}, 
        geom_info::Libsharp.GeomInfo, alm_info::Libsharp.AlmInfo, 
        niter::Integer, spin::Integer)

    ncomp = size(maps, 1)
    residual_maps = [zero(m) for m in maps] # set up buffers
    residual_alms = [zero(a) for a in alms]
    for iter_index in 1:niter
        # synthesize map from the alms
        Libsharp.sharp_execute!(
            Libsharp.SHARP_ALM2MAP, spin, alms, residual_maps,
            geom_info, alm_info, Libsharp.SHARP_DP)
        for i_comp in 1:ncomp
            # subtract synthesized map from original
            residual_maps[i_comp] .= maps[i_comp] .- residual_maps[i_comp]
        end
        # synthesize alms from the residual maps
        Libsharp.sharp_execute!(
            Libsharp.SHARP_MAP2ALM, spin, residual_alms, residual_maps,
            geom_info, alm_info, Libsharp.SHARP_DP)
        for i_comp in 1:ncomp
            # add them to the alms
            alms[i_comp] .= alms[i_comp] .+ residual_alms[i_comp]
        end
    end
end


"""
    map2alm!(map::Map{Float64, RingOrder}, alm::Alm{ComplexF64}; niter::Integer=3)
    map2alm!(map::PolarizedMap{Float64, RingOrder}, alm::Array{Alm{ComplexF64},1}; 
        niter::Integer=3)

This function performs a spherical harmonic transform on the map and places the results
in the passed `alm` object. This function requires types derived from Float64, since it is
done in-place.

# Arguments
- `map::Union{Map{Float64, RingOrder}, PolarizedMap{Float64, RingOrder}`: the map that to 
    perform the spherical harmonic transform on.
- `alm::Alm{ComplexF64}`: the spherical harmonic coefficients to be written to.

# Keywords
- `niter::Integer`: number of iterations of SHTs to perform, to enhance accuracy
"""
function map2alm! end

function map2alm!(map::Map{Float64, RingOrder}, alm::Alm{ComplexF64}; niter::Integer=3)
    geom_info = Libsharp.make_healpix_geom_info(map.resolution.nside, 1)
    alm_info = Libsharp.make_triangular_alm_info(alm.lmax, alm.mmax, 1)
    Libsharp.sharp_execute!(
        Libsharp.SHARP_MAP2ALM, 0, [alm.alm], [map.pixels],
        geom_info, alm_info, Libsharp.SHARP_DP)
    if niter > 0
        iterate_map2alm!([map.pixels], [alm.alm], geom_info, alm_info, niter, 0)
    end
end

function map2alm!(map::PolarizedMap{Float64, RingOrder}, alm::Array{Alm{ComplexF64},1}; 
        niter::Integer=3)
    geom_info = Libsharp.make_healpix_geom_info(map.i.resolution.nside, 1)
    alm_info = Libsharp.make_triangular_alm_info(alm[1].lmax, alm[1].mmax, 1)

    maps_0 = [map.i.pixels]
    maps_2 = [map.q.pixels, map.u.pixels]
    alms_0 = [alm[1].alm]
    alms_2 = [alm[2].alm, alm[3].alm]
    Libsharp.sharp_execute!(
        Libsharp.SHARP_MAP2ALM, 0, alms_0, [map.i.pixels],
        geom_info, alm_info, Libsharp.SHARP_DP)
    Libsharp.sharp_execute!(
        Libsharp.SHARP_MAP2ALM, 2, alms_2, maps_2,
        geom_info, alm_info, Libsharp.SHARP_DP)

    if niter > 0
        iterate_map2alm!(maps_0, alms_0, geom_info, alm_info, niter, 0)
        iterate_map2alm!(maps_2, alms_2, geom_info, alm_info, niter, 2)
    end
end


"""
    map2alm(map::Map{Float64, RingOrder}; 
        lmax=nothing, mmax=nothing, niter::Integer=3)
    map2alm(m::Map{T, RingOrder}; lmax=nothing, mmax=nothing, 
        niter::Integer=3) where T <: Real

Compute the spherical harmonic coefficients of a map. To enhance precision, more iterations
of the transforms can be performed by passing a nonzero `niter`. The underlying SHT library
libsharp performs all calculations in Cdouble, so all inputs are converted to types based on
Float64.

# Arguments
- `map::Union{Map{T, RingOrder}, PolarizedMap{T, RingOrder}`: the map that to 
    perform the spherical harmonic transform on.

# Keywords
- `lmax::Integer`: the maximum ℓ coefficient, will default to 3*nside-1 if not specified.
- `mmax::Integer`: the maximum m coefficient
- `niter::Integer`: number of SHT iterations, to enhance precision. Defaults to 3

# Returns
- `Alm{ComplexF64}`: the spherical harmonic coefficients corresponding to the map
"""
function map2alm end

function map2alm(map::Map{Float64, RingOrder}; 
        lmax=nothing, mmax=nothing, niter::Integer=3)

    nside = map.resolution.nside
    lmax = isnothing(lmax) ? 3 * nside - 1 : lmax
    mmax = isnothing(mmax) ? lmax : mmax
    nalms = numberOfAlms(lmax, mmax)
    alm = Alm{ComplexF64}(lmax, mmax, zeros(ComplexF64, nalms))

    map2alm!(map, alm; niter=niter)
    return alm
end

function map2alm(map::PolarizedMap{Float64, RingOrder};
        lmax=nothing, mmax=nothing, niter::Integer=3)

    nside = map.i.resolution.nside
    lmax = isnothing(lmax) ? 3 * nside - 1 : lmax
    mmax = isnothing(mmax) ? lmax : mmax
    nalms = numberOfAlms(lmax, mmax)
    alms = [Alm{ComplexF64}(lmax, mmax, zeros(ComplexF64, nalms)),
            Alm{ComplexF64}(lmax, mmax, zeros(ComplexF64, nalms)),
            Alm{ComplexF64}(lmax, mmax, zeros(ComplexF64, nalms))]

    map2alm!(map, alms; niter=niter)
    return alms
end

# convert maps to Float64
function map2alm(map::Map{T, RingOrder}; lmax=nothing, mmax=nothing, 
        niter::Integer=3) where T <: Real
    map_float = Map{Float64, RingOrder}(convert(Array{Float64,1}, map.pixels))
    return map2alm(map_float, lmax=lmax, mmax=mmax, niter=niter)
end

# convert PolarizedMap to Float64
function map2alm(map::PolarizedMap{T, RingOrder}; lmax=nothing, mmax=nothing, 
                 niter::Integer=3) where T <: Real
    m_i = convert(Array{Float64,1}, map.i)
    m_q = convert(Array{Float64,1}, map.q)
    m_u = convert(Array{Float64,1}, map.u)
    pol_map_float = PolarizedMap{Float64, RingOrder}(m_i, m_q, m_u)
    return map2alm(pol_map_float, lmax=lmax, mmax=mmax, niter=niter)
end


"""
    alm2map!(alm::Alm{ComplexF64}, map::Map{Float64, RingOrder})
    alm2map!(alm::Array{Alm{ComplexF64},1}, map::PolarizedMap{Float64, RingOrder})

This function performs a spherical harmonic transform on the map and places the results
in the passed `alm` object. This function requires types derived from Float64, since it is
done in-place.

# Arguments
- `alm::Alm{ComplexF64}`: the spherical harmonic coefficients to perform the spherical 
    harmonic transform on.
- `map::Union{Map{Float64, RingOrder}, PolarizedMap{Float64, RingOrder}`: the map to be 
    written to.
"""
function alm2map! end

# in-place alm2map for spin-0
function alm2map!(alm::Alm{ComplexF64}, map::Map{Float64, RingOrder})
    geom_info = Libsharp.make_healpix_geom_info(map.resolution.nside, 1) 
    alm_info = Libsharp.make_triangular_alm_info(alm.lmax, alm.mmax, 1)
    Libsharp.sharp_execute!(
        Libsharp.SHARP_ALM2MAP, 0, [alm.alm], [map.pixels],
        geom_info, alm_info, Libsharp.SHARP_DP)
end

# in-place alm2map for TEB to IQU
function alm2map!(alm::Array{Alm{ComplexF64},1}, map::PolarizedMap{Float64, RingOrder})
    geom_info = Libsharp.make_healpix_geom_info(map.i.resolution.nside, 1) 
    alm_info = Libsharp.make_triangular_alm_info(alm[1].lmax, alm[1].mmax, 1)

    Libsharp.sharp_execute!(
        Libsharp.SHARP_ALM2MAP, 0, 
        [alm[1].alm], [map.i.pixels],
        geom_info, alm_info, Libsharp.SHARP_DP)

    Libsharp.sharp_execute!(
        Libsharp.SHARP_ALM2MAP, 2, 
        [alm[2].alm, alm[3].alm],  [map.q.pixels, map.u.pixels],
        geom_info, alm_info, Libsharp.SHARP_DP)
end


"""
    alm2map(alm::Alm{ComplexF64}, nside::Integer)
    alm2map(alm::Alm{T}, nside::Integer) where T
    alm2map(alm::Array{Alm{ComplexF64},1}, nside::Integer)
    alm2map(alms::Array{Alm{T},1}, nside::Integer) where T

Compute a map from spherical harmonic coefficients. The underlying SHT library libsharp 
performs all calculations in Cdouble, so all inputs are converted to types based on Float64.

# Arguments
- `alm`: the spherical harmonic coefficients to transform. If of type `Alm{T}`, we assume a 
    spin-0 spherical harmonic transform. If an array of `Alm` is passed, we assume that
    the components correspond to T, E, and B coefficients.

# Keywords
- `nside::Integer`: Healpix resolution parameter

# Returns
- `Map{Float64, RingOrder}` or `PolarizedMap{Float64, RingOrder}` depending on if the input
    alm is of type `Alm{T}` or `Array{Alm{T}}` respectively.
"""
function alm2map end

# create a new set of spin-0 maps and project the coefficients to the map
function alm2map(alm::Alm{ComplexF64}, nside::Integer)
    npix = nside2npix(nside)
    map = Map{Float64, RingOrder}(zeros(Float64, npix))
    alm2map!(alm, map)
    return map
end

# create a new set of IQU maps and project the coefficients to the map
function alm2map(alm::Array{Alm{ComplexF64},1}, nside::Integer)
    npix = nside2npix(nside)
    map = PolarizedMap{Float64, RingOrder}(
        zeros(Float64, npix), zeros(Float64, npix), zeros(Float64, npix))
    alm2map!(alm, map)
    return map
end

# convert to ComplexF64 Alm for spin-0 if passed some other type
function alm2map(alm::Alm{T}, nside::Integer) where T
    alm_float = Alm{ComplexF64}(
        alm.lmax, alm.mmax, convert(Array{ComplexF64,1}, alm.alm))
    return alm2map(alm_float, nside)
end

# convert to ComplexF64 Alm for TEB if passed some other type
function alm2map(alms::Array{Alm{Complex{T}},1}, nside::Integer) where T <: Real
    lmax = alms[1].lmax
    mmax = alms[1].mmax
    alm_t = Alm{ComplexF64}(lmax, mmax, 
        convert(Array{ComplexF64,1}, alms[1].alm))
    alm_e = Alm{ComplexF64}(lmax, mmax, 
        convert(Array{ComplexF64,1}, alms[2].alm))
    alm_b = Alm{ComplexF64}(lmax, mmax, 
        convert(Array{ComplexF64,1}, alms[3].alm))
    return alm2map([alm_t, alm_e, alm_b], nside)
end

raw"""
    alm2cl(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}

Compute ``C_{\ell}`` from the spherical harmonic coefficients of two fields.

# Arguments
- `alm₁::Alm{Complex{T}}`: the spherical harmonic coefficients of the first field
- `alm₂::Alm{Complex{T}}`: the spherical harmonic coefficients of the second field

# Returns
- `Array{T}` containing C_{\ell}, with the first element referring to ℓ=0.
"""
function alm2cl(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}
    (alm₁.lmax != alm₂.lmax) && throw(ArgumentError("Alm lmax do not match."))
    (alm₁.mmax != alm₂.mmax) && throw(ArgumentError("Alm mmax do not match."))
    (alm₁.mmax < alm₂.lmax) && throw(ArgumentError("Alm mmax < lmax."))

    lmax = alm₁.lmax
    cl = zeros(T, lmax+1)
    for l in 0:lmax
        for m in 1:l 
            index = almIndex(alm₁, l, m)
            cl[l+1] += 2 * alm₁.alm[index] * conj(alm₂.alm[index])
        end
        index0 = almIndex(alm₁, l, 0)
        cl[l+1] += alm₁.alm[index0] * conj(alm₂.alm[index0])
        cl[l+1] = cl[l+1] / (2 * l + 1)
    end
    return cl
end
alm2cl(alm::Alm{Complex{T}}) where {T <: Number} = alm2cl(alm, alm)

