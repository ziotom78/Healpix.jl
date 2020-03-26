
# internal function for iterating map2alm
function iterate_map2alm!(
        maps::Array{Array{Float64,1},1}, 
        alms::Array{Array{Complex{Float64},1},1}, 
        geom_info::Libsharp.GeomInfo, alm_info::Libsharp.AlmInfo, 
        niter::Integer, spin::Integer)

    ncomp = size(maps, 1)
    residual_maps = [zero(m) for m in maps] # set up buffers
    residual_alms = [zero(a) for a in alms]
    for iter_index in 1:niter
        Libsharp.sharp_execute!(
            Libsharp.SHARP_ALM2MAP, spin, alms, residual_maps,
            geom_info, alm_info, Libsharp.SHARP_DP)
        for i_comp in 1:ncomp
            residual_maps[i_comp] .= maps[i_comp] .- residual_maps[i_comp]
        end
        Libsharp.sharp_execute!(
            Libsharp.SHARP_MAP2ALM, spin, residual_alms, residual_maps,
            geom_info, alm_info, Libsharp.SHARP_DP)
        for i_comp in 1:ncomp
            alms[i_comp] .= alms[i_comp] .+ residual_alms[i_comp]
        end
    end
end

"""
map2alm(...)

Compute the spherical harmonic coefficients of a map. This will
convert your map to Float64 if it isn't already, and will output a ComplexF64 Alm.
"""
function map2alm end

"""
map2alm!(...)

Compute in-place the spherical harmonic coefficients of a map.
"""
function map2alm! end

function map2alm!(map::Map{Float64, RingOrder}, alm::Alm{ComplexF64}; 
        niter::Integer=3)
    geom_info = Libsharp.make_healpix_geom_info(map.resolution.nside, 1)
    alm_info = Libsharp.make_triangular_alm_info(alm.lmax, alm.mmax, 1)
    Libsharp.sharp_execute!(
        Libsharp.SHARP_MAP2ALM, 0, [alm.alm], [map.pixels],
        geom_info, alm_info, Libsharp.SHARP_DP)
    if niter > 0
        iterate_map2alm!([map.pixels], [alm.alm], geom_info, alm_info, niter, 0)
    end
end

function map2alm(map::Map{Float64, RingOrder}; 
        lmax::Integer=nothing, mmax::Integer=nothing, niter::Integer=3)

    nside = map.resolution.nside
    lmax = isnothing(lmax) ? 3 * nside - 1 : lmax
    mmax = isnothing(mmax) ? lmax : mmax
    nalms = numberOfAlms(lmax, mmax)
    alm = Alm{ComplexF64}(lmax, mmax, zeros(ComplexF64, nalms))

    map2alm!(map, alm; niter=niter)
    return alm
end

function map2alm(m::Map{T, RingOrder}; 
        lmax::Integer=nothing, mmax::Integer=nothing, niter::Integer=3) where T
    # convert Map to Float64
    map_float = Map{Float64, RingOrder}(convert(Array{Float64,1}, m.pixels))
    return map2alm(map_float, lmax=lmax, mmax=mmax, niter=niter)
end

function map2alm!(
        map::PolarizedMap{Float64, RingOrder}, 
        alm::Array{Alm{ComplexF64},1}; 
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

function map2alm(map::PolarizedMap{Float64, RingOrder};
        lmax::Integer=nothing, mmax::Integer=nothing, niter::Integer=3)

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

function map2alm(map::PolarizedMap{T, RingOrder};
        lmax::Integer=nothing, mmax::Integer=nothing, niter::Integer=3) where T
    # convert PolarizedMap to Float64
    m_i = convert(Array{Float64,1}, map.i)
    m_q = convert(Array{Float64,1}, map.q)
    m_u = convert(Array{Float64,1}, map.u)
    pol_map_float = PolarizedMap{Float64, RingOrder}(m_i, m_q, m_u)
    return map2alm(pol_map_float, lmax=lmax, mmax=mmax, niter=niter)
end

"""
    alm2map(m::Alm{T}, lmax::Integer, mmax::Integer)

Compute a map from spherical harmonic coefficients.
"""
function alm2map end

"""
    alm2map!(...)

Compute in-place a map from spherical harmonic coefficients.
"""
function alm2map! end

function alm2map!(alm::Alm{ComplexF64}, map::Map{Float64, RingOrder})
     # unweighted healpix map
    geom_info = Libsharp.make_healpix_geom_info(map.resolution.nside, 1) 
    alm_info = Libsharp.make_triangular_alm_info(alm.lmax, alm.mmax, 1)
    Libsharp.sharp_execute!(
        Libsharp.SHARP_ALM2MAP, 0, [alm.alm], [map.pixels],
        geom_info, alm_info, Libsharp.SHARP_DP)
end

function alm2map(alm::Alm{ComplexF64}, nside::Integer)
    npix = nside2npix(nside)
    map = Map{Float64, RingOrder}(zeros(Float64, npix))
    alm2map!(alm, map)
    return map
end

function alm2map(alm::Alm{T}, nside::Integer) where T
    alm_float = Alm{ComplexF64}(
        alm.lmax, alm.mmax, convert(Array{ComplexF64,1}, alm.alm))
    return alm2map(alm_float, nside)
end

function alm2map!(
        alm::Array{Alm{ComplexF64},1}, 
        map::PolarizedMap{Float64, RingOrder})
    # unweighted healpix map
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

function alm2map(alm::Array{Alm{ComplexF64},1}, nside::Integer)
    npix = nside2npix(nside)
    map = PolarizedMap{Float64, RingOrder}(
        zeros(Float64, npix), zeros(Float64, npix), zeros(Float64, npix))
    alm2map!(alm, map)
    return map
end

function alm2map(alms::Array{Alm{T},1}, nside::Integer) where T
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

