

"""
    map2alm(...)

Compute the spherical harmonic coefficients of a real scalar map. This will
convert your map to Float64 if it isn't already, and will output a ComplexF64 Alm.
"""
function map2alm end

function map2alm(m::Map{T, O}, lmax::Integer, mmax::Integer) where {T, O}
    nside = m.resolution.nside

    geom_info = Libsharp.make_healpix_geom_info(nside, 1)  # unweighted healpix map
    alm_info = Libsharp.make_triangular_alm_info(lmax, mmax, 1)

    npix = Libsharp.map_size(geom_info)
    nalms = Libsharp.alm_count(alm_info)

    alms = [ones(ComplexF64, nalms)]
    maps = [convert(Array{Float64,1}, m.pixels)]

    spin = 0
    Libsharp.sharp_execute!(
        Libsharp.SHARP_MAP2ALM, spin, alms, maps,
        geom_info, alm_info, Libsharp.SHARP_DP
    )

    return Alm{ComplexF64}(lmax, mmax, alms[1])
end

map2alm(m::Map{T, O}, lmax::Integer) where {T, O} = map2alm(m, lmax, lmax)
map2alm(m::Map{T, O}) where {T, O} = map2alm(m, 3 * m.resolution.nside - 1)


"""
    alm2map(m::Map{T, O}, lmax::Integer, mmax::Integer) where {T, O}

Compute a real scalar map from spherical harmonic coefficients.
"""
function alm2map(alm::Alm{T}, nside::Integer) where T
    lmax = alm.lmax
    mmax = alm.mmax

    geom_info = Libsharp.make_healpix_geom_info(nside, 1)  # unweighted healpix map
    alm_info = Libsharp.make_triangular_alm_info(lmax, mmax, 1)

    npix = Libsharp.map_size(geom_info)
    alms = Libsharp.alm_count(alm_info)

    alm_ComplexF64 = convert(Array{ComplexF64,1}, alm.alm)
    alms = [alm_ComplexF64]
    maps = [zeros(Float64, npix)]

    spin = 0
    Libsharp.sharp_execute!(
        Libsharp.SHARP_ALM2MAP, spin, alms, maps,
        geom_info, alm_info, Libsharp.SHARP_DP
    )

    return Map{Float64, RingOrder}(maps[1])
end