
"""
    map2alm(...)

Compute the spherical harmonic coefficients of a map. This will
convert your map to Float64 if it isn't already, and will output a ComplexF64 Alm.
"""
function map2alm end


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

function map2alm(m::Map{Float64, RingOrder}; 
        lmax::Integer=nothing, mmax::Integer=nothing, niter::Integer=3)

    nside = m.resolution.nside
    lmax = isnothing(lmax) ? 3 * nside - 1 : lmax
    mmax = isnothing(mmax) ? lmax : mmax

    geom_info = Libsharp.make_healpix_geom_info(nside, 1)  # unweighted healpix map
    alm_info = Libsharp.make_triangular_alm_info(lmax, mmax, 1)

    npix = Libsharp.map_size(geom_info)
    nalms = Libsharp.alm_count(alm_info)

    alms = [zeros(ComplexF64, nalms)]
    maps = [m.pixels]

    Libsharp.sharp_execute!(
        Libsharp.SHARP_MAP2ALM, 0, alms, maps,
        geom_info, alm_info, Libsharp.SHARP_DP)

    if niter > 0
        iterate_map2alm!(maps, alms, geom_info, alm_info, niter, 0)
    end

    return Alm{ComplexF64}(lmax, mmax, alms[1])
end

function map2alm(m::Map{T, RingOrder}; 
        lmax::Integer=nothing, mmax::Integer=nothing, niter::Integer=3) where T
    # convert Map to Float64
    map_float = Map{Float64, RingOrder}(convert(Array{Float64,1}, m.pixels))
    return map2alm(map_float, lmax=lmax, mmax=mmax, niter=niter)
end

function map2alm(m::PolarizedMap{Float64, RingOrder};
        lmax::Integer=nothing, mmax::Integer=nothing, niter::Integer=3)

    nside = m.i.resolution.nside
    lmax = isnothing(lmax) ? 3 * nside - 1 : lmax
    mmax = isnothing(mmax) ? lmax : mmax

    geom_info = Libsharp.make_healpix_geom_info(nside, 1)  # unweighted healpix map
    alm_info = Libsharp.make_triangular_alm_info(lmax, mmax, 1)
    npix = Libsharp.map_size(geom_info)
    nalms = Libsharp.alm_count(alm_info)

    maps_0 = [m.i.pixels]
    alms_0 = [zeros(ComplexF64, nalms)]
    maps_2 = [m.q.pixels, m.u.pixels]
    alms_2 = [zeros(ComplexF64, nalms), zeros(ComplexF64, nalms)]

    Libsharp.sharp_execute!(
        Libsharp.SHARP_MAP2ALM, 0, alms_0, maps_0,
        geom_info, alm_info, Libsharp.SHARP_DP)

    Libsharp.sharp_execute!(
        Libsharp.SHARP_MAP2ALM, 2, alms_2, maps_2,
        geom_info, alm_info, Libsharp.SHARP_DP)

    if niter > 0
        iterate_map2alm!(maps_0, alms_0, geom_info, alm_info, niter, 0)
        iterate_map2alm!(maps_2, alms_2, geom_info, alm_info, niter, 2)
    end

    return [
        Alm{ComplexF64}(lmax, mmax, alms_0[1]),
        Alm{ComplexF64}(lmax, mmax, alms_2[1]),
        Alm{ComplexF64}(lmax, mmax, alms_2[2])
    ]
end

function map2alm(m::PolarizedMap{T, RingOrder};
        lmax::Integer=nothing, mmax::Integer=nothing, niter::Integer=3) where T
    # convert PolarizedMap to Float64
    m_i = convert(Array{Float64,1}, m.i)
    m_q = convert(Array{Float64,1}, m.q)
    m_u = convert(Array{Float64,1}, m.u)
    pol_map_float = PolarizedMap{Float64, RingOrder}(m_i, m_q, m_u)
    return map2alm(pol_map_float, lmax=lmax, mmax=mmax, niter=niter)
end

"""
    alm2map(m::Alm{T}, lmax::Integer, mmax::Integer)

Compute a map from spherical harmonic coefficients.
"""
function alm2map end

function alm2map(alm::Alm{ComplexF64}, nside::Integer)
    lmax = alm.lmax
    mmax = alm.mmax

    geom_info = Libsharp.make_healpix_geom_info(nside, 1)  # unweighted healpix map
    alm_info = Libsharp.make_triangular_alm_info(lmax, mmax, 1)

    npix = Libsharp.map_size(geom_info)
    nalms = Libsharp.alm_count(alm_info)

    alm_ComplexF64 = convert(Array{ComplexF64,1}, alm.alm)
    alms = [alm_ComplexF64]
    maps = [zeros(Float64, npix)]

    Libsharp.sharp_execute!(
        Libsharp.SHARP_ALM2MAP, 0, alms, maps,
        geom_info, alm_info, Libsharp.SHARP_DP)

    return Map{Float64, RingOrder}(maps[1])
end

function alm2map(alm::Alm{T}, nside::Integer) where T
    alm_float = Alm{ComplexF64}(
        alm.lmax, alm.mmax, convert(Array{ComplexF64,1}, alm.alm))
    return alm2map(alm_float, nside)
end


function alm2map(alms::Array{Alm{ComplexF64},1}, nside::Integer)
    lmax = alms[1].lmax
    mmax = alms[1].mmax

    # unweighted healpix map
    geom_info = Libsharp.make_healpix_geom_info(nside, 1) 
    alm_info = Libsharp.make_triangular_alm_info(lmax, mmax, 1)
    npix = Libsharp.map_size(geom_info)
    nalms = Libsharp.alm_count(alm_info)

    maps_0 = [zeros(Float64, npix)]
    maps_2 = [zeros(Float64, npix), zeros(Float64, npix)]
    alms_0 = [alms[1].alm]
    alms_2 = [alms[2].alm, alms[3].alm]

    Libsharp.sharp_execute!(
        Libsharp.SHARP_ALM2MAP, 0, alms_0, maps_0,
        geom_info, alm_info, Libsharp.SHARP_DP)

    Libsharp.sharp_execute!(
        Libsharp.SHARP_ALM2MAP, 2, alms_2, maps_2,
        geom_info, alm_info, Libsharp.SHARP_DP)

    return PolarizedMap{Float64, RingOrder}(
        maps_0[1], maps_2[1], maps_2[2])
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

