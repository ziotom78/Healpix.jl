
"""
    map2alm(...)

Compute the spherical harmonic coefficients of a map. This will
convert your map to Float64 if it isn't already, and will output a ComplexF64 Alm.
"""
function map2alm end

function map2alm(m::Map{T, RingOrder}, 
        lmax::Integer, mmax::Integer) where T
    nside = m.resolution.nside

    geom_info = Libsharp.make_healpix_geom_info(nside, 1)  # unweighted healpix map
    alm_info = Libsharp.make_triangular_alm_info(lmax, mmax, 1)

    npix = Libsharp.map_size(geom_info)
    nalms = Libsharp.alm_count(alm_info)

    alms = [zeros(ComplexF64, nalms)]
    maps = [convert(Array{Float64,1}, m.pixels)]

    Libsharp.sharp_execute!(
        Libsharp.SHARP_MAP2ALM, 0, alms, maps,
        geom_info, alm_info, Libsharp.SHARP_DP)

    return Alm{ComplexF64}(lmax, mmax, alms[1])
end

map2alm(m::Map{T, RingOrder}, lmax::Integer) where T = map2alm(m, lmax, lmax)
map2alm(m::Map{T, RingOrder}) where T = map2alm(m, 3 * m.resolution.nside - 1)


function map2alm(m::PolarizedMap{T, RingOrder}, 
        lmax::Integer, mmax::Integer) where T

    nside = m.i.resolution.nside
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

    return [
        Alm{ComplexF64}(lmax, mmax, alms_0[1]),
        Alm{ComplexF64}(lmax, mmax, alms_2[1]),
        Alm{ComplexF64}(lmax, mmax, alms_2[2])
    ]
end

map2alm(m::PolarizedMap{T, RingOrder}, lmax::Integer) where T = map2alm(m, lmax, lmax)
map2alm(m::PolarizedMap{T, RingOrder}) where T = map2alm(m, 3 * m.resolution.nside - 1)


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

