
"""
    map2alm(...)

Compute the spherical harmonic coefficients of a map. This will
convert your map to Float64 if it isn't already, and will output a ComplexF64 Alm.
"""
function map2alm end


# internal function for iterating map2alm
function iterate_map2alm!(maps::Array{Array{T,1},1}, alms::Array{Array{Complex{T},1},1}, 
        geom_info::Libsharp.GeomInfo, alm_info::Libsharp.AlmInfo, 
        niter::Integer, spin::Integer) where T

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

function map2alm(m::Map{T, RingOrder}; 
        lmax::Integer=nothing, mmax::Integer=nothing, niter::Integer=3) where T

    nside = m.resolution.nside
    lmax = isnothing(lmax) ? 3 * nside - 1 : lmax
    mmax = isnothing(mmax) ? lmax : mmax

    geom_info = Libsharp.make_healpix_geom_info(nside, 1)  # unweighted healpix map
    alm_info = Libsharp.make_triangular_alm_info(lmax, mmax, 1)

    npix = Libsharp.map_size(geom_info)
    nalms = Libsharp.alm_count(alm_info)

    alms = [zeros(ComplexF64, nalms)]
    maps = [convert(Array{Float64,1}, m.pixels)]

    Libsharp.sharp_execute!(
        Libsharp.SHARP_MAP2ALM, 0, alms, maps,
        geom_info, alm_info, Libsharp.SHARP_DP)

    if niter > 0
        iterate_map2alm!(maps, alms, geom_info, alm_info, niter, 0)
    end

    return Alm{ComplexF64}(lmax, mmax, alms[1])
end


function map2alm(m::PolarizedMap{T, RingOrder};
        lmax::Integer=nothing, mmax::Integer=nothing, niter::Integer=3) where T

    nside = m.i.resolution.nside
    lmax = isnothing(lmax) ? 3 * nside - 1 : lmax
    mmax = isnothing(mmax) ? lmax : mmax

    geom_info = Libsharp.make_healpix_geom_info(nside, 1)  # unweighted healpix map
    alm_info = Libsharp.make_triangular_alm_info(lmax, mmax, 1)
    npix = Libsharp.map_size(geom_info)
    nalms = Libsharp.alm_count(alm_info)

    m_i = convert(Array{Float64,1}, m.i)
    m_q = convert(Array{Float64,1}, m.q)
    m_u = convert(Array{Float64,1}, m.u)

    maps_0 = [m_i]
    alms_0 = [zeros(ComplexF64, nalms)]
    maps_2 = [m_q, m_u]
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


"""
    alm2map(m::Alm{T}, lmax::Integer, mmax::Integer)

Compute a map from spherical harmonic coefficients.
"""
function alm2map end

function alm2map(alm::Alm{T}, nside::Integer) where T
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


function alm2map(alms::Array{Alm{T},1}, nside::Integer) where T
    lmax = alms[1].lmax
    mmax = alms[1].mmax

    # unweighted healpix map
    geom_info = Libsharp.make_healpix_geom_info(nside, 1) 
    alm_info = Libsharp.make_triangular_alm_info(lmax, mmax, 1)

    npix = Libsharp.map_size(geom_info)
    nalms = Libsharp.alm_count(alm_info)

    alm_t = convert(Array{ComplexF64,1}, alms[1].alm)
    alm_e = convert(Array{ComplexF64,1}, alms[2].alm)
    alm_b = convert(Array{ComplexF64,1}, alms[3].alm)

    maps_0 = [zeros(Float64, npix)]
    maps_2 = [zeros(Float64, npix), zeros(Float64, npix)]
    alms_0 = [alm_t]
    alms_2 = [alm_e, alm_b]

    Libsharp.sharp_execute!(
        Libsharp.SHARP_ALM2MAP, 0, alms_0, maps_0,
        geom_info, alm_info, Libsharp.SHARP_DP)

    Libsharp.sharp_execute!(
        Libsharp.SHARP_ALM2MAP, 2, alms_2, maps_2,
        geom_info, alm_info, Libsharp.SHARP_DP)

    return PolarizedMap{Float64, RingOrder}(
        maps_0[1], maps_2[1], maps_2[2])
end

