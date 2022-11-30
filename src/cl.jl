"""
    readClFromFITS{T <: Real}(f::CFITSIO.FITSFile, t::Type{T}; col_num = 2) -> Vector{T}
    readClFromFITS{T <: Real}(fileName::String, t::Type{T}; col_num = 2) -> Vector{T}
Read a set of C_ℓ coefficients from a FITS file.
"""
function readClFromFITS(f::CFITSIO.FITSFile, t::Type{T}; col_num = 2) where {T <: Real}
    numOfRows = CFITSIO.fits_get_num_rows(f)

    Cl = Vector{T}(undef, numOfRows)

    CFITSIO.fits_read_col(f, col_num, 1, 1, Cl)

    return Cl
end

function readClFromFITS(fileName, t::Type{T}; col_num = 2) where {T <: Real}
    f = CFITSIO.fits_open_table(fileName)
    try
        result = readClFromFITS(f, t; col_num = 2)
        return result
    finally
        CFITSIO.fits_close_file(f)
end
end


##############################################################################

"""
    writeClToFITS(f::CFITSIO.FITSFile, Cl::Vector{T}) where {T <: Real}
    writeClToFITS(fileName, Cl::Vector{T}; overwrite = true) where {T <: Real}
Write a set of C_ℓ coefficients to a FITS file.
"""
function writeClToFITS(f::CFITSIO.FITSFile, Cl::Vector{T}) where {T <: Real}

    idx = Vector{Int64}(1:length(Cl))

    CFITSIO.fits_write_col(f, 1, 1, 1, idx)
    CFITSIO.fits_write_col(f, 2, 1, 1, Cl)

end

function writeClToFITS(fileName, Cl::Vector{T}; overwrite = true) where {T <: Real}
    if overwrite
        f = CFITSIO.fits_clobber_file(fileName)
    else
        f = CFITSIO.fits_create_file(fileName)
    end
    try
        CFITSIO.fits_create_binary_tbl(f, 0, [("Index", "1I", ""),("C_l", "1D", "")], "Cl")
        writeClToFITS(f, Cl)
    finally
    end
    CFITSIO.fits_close_file(f)
end


#########################################################################

"""
    dl2cl(dl::AbstractVector{T}, lmin::Integer) where {T <: Real}

Convert a set of ``D_{\\ell}`` to ``C_{\\ell}`` power spectrum, where
``C_{\\ell} = 2\\pi D_{\\ell} / \\ell (\\ell + 1)``. The first components are
set to zero if not present. The monopole component is set to zero in any case to avoid Inf values.

# Arguments:
- `dl::AbstractVector{T}` : Array of D_ℓ components
- `lmin::Integer` : minimum l in the representation of the Dℓ power spectrum

# Returns:
- `Vector{T}` : Array of C_ℓ power spectrum components
"""
function dl2cl(dl::AbstractVector{T}, lmin::Integer) where {T <: Real}
    (lmin >= 0) || throw(DomainError(lmin, "`lmin` is not positive or zero"))
    lmax = length(dl)+lmin-1
    l_s = Vector{Int}(0:lmax)
    #fill the missing initial components (monopole and/or dipole, ...) with zeros
    head = zeros(lmin)
    dl = append!(head, dl)
    cl =  dl .* 2π ./ (l_s .* (l_s .+ 1))
    cl[1] = 0
    return cl
end


"""
    cl2dl(cl::AbstractVector{T}, lmin::Integer) where {T <: Real}

Convert a set of ``C_{\\ell}`` to ``D_{\\ell}`` power spectrum, where
``D_{\\ell} = \\ell (\\ell + 1) C_{\\ell} / 2\\pi``.
The first components are set to zero if not present.

# Arguments:
- `cl::AbstractVector{T}` : Array of C_ℓ components
- `lmin::Integer` : minimum l in the representation of the C_ℓ power spectrum

# Returns:
- `Vector{T}` : Array of D_ℓ power spectrum components
"""
function cl2dl(cl::AbstractVector{T}, lmin::Integer) where {T <: Real}
    (lmin >= 0) || throw(DomainError(lmin, "`lmin` is not positive or zero"))
    lmax = length(cl)+lmin-1
    l_s = Vector{Int}(0:lmax)
    head = zeros(lmin)
    #fill the missing initial components (monopole and/or dipole) with zeros
    cl = append!(head, cl)
    dl =  cl ./ 2π .* (l_s .* (l_s .+ 1))
    return dl
end

##########################################################################

"""
    synalm!(cl::Vector{T}, alm::Alm{ComplexF64, Vector{ComplexF64}}, rng::AbstractRNG) where {T <: Real}
    synalm!(cl::Vector{T}, alm::Alm{ComplexF64, Vector{ComplexF64}}) where {T <: Real}

Generate a set of ``a_{\\ell m}`` from a given power spectra ``C_{\\ell}``.
The output is written into the `Alm` object passed in input.

<<<<<<< HEAD
# Arguments:
=======
# ARGUMENTS
>>>>>>> 873b35eff1611e0a5f8921640a037354a4f82bb7
- `cl::AbstractVector{T}`: The array representing the power spectrum components ``C_{\\ell}``,
starting from `` \\ell = 0 ``.
- `alm::Alm{Complex{T}}`: The array representing the spherical harmonics coefficients ``a_{\\ell m}``
we want to write the result into.
- `rng::AbstractRNG` : (optional) the RNG to be used for generating the ``a_{\\ell m}``. It allows
to set the seed beforehand guaranteeing the reproducibility of the process.
"""
function synalm!(cl::Vector{T}, alm::Alm{ComplexF64, Vector{ComplexF64}}, rng::AbstractRNG) where {T <: Real}
    cl_size = length(cl)
    lmax = alm.lmax
    mmax = alm.mmax
    (cl_size - 1 >= lmax) || throw(DomainError(cl_size, "not enough C_l's to generate Alm"))

    for l = 0:lmax
        if l <= mmax
            maxm = l
        else
            maxm = mmax
        end
        for m = 0:maxm
            i = almIndex(alm, l, m)
            #for m=0 the alm must be real, since alm^R_l,0 = alm^C_l,0, if the field is real!
            alm.alm[i] = randn(rng, ifelse(m > 0, ComplexF64, Float64))*sqrt(cl[l+1]) #sqrt bc it's the variance
        end
    end
end

synalm!(cl::Vector{T}, alm::Alm{ComplexF64, Vector{ComplexF64}}) where {T <: Real} =
    synalm!(cl, alm, Random.GLOBAL_RNG)

"""
    synalm(cl::Vector{T}, lmax::Integer, mmax::Integer, rng::AbstractRNG) where {T <: Real}
    synalm(cl::Vector{T}, lmax::Integer, mmax::Integer) where {T <: Real}
    synalm(cl::Vector{T}, lmax::Integer, rng::AbstractRNG) where {T <: Real}
    synalm(cl::Vector{T}, lmax::Integer) where {T <: Real}
    synalm(cl::Vector{T}, rng::AbstractRNG) where {T <: Real}
    synalm(cl::Vector{T}) where {T <: Real}

Generate a set of ``a_{\\ell m}`` from a given power spectra ``C_{\\ell}``.
The output is written into a new `Alm` object of given lmax.

<<<<<<< HEAD
# Arguments:
=======
# ARGUMENTS
>>>>>>> 873b35eff1611e0a5f8921640a037354a4f82bb7
- `cl::AbstractVector{T}`: The array representing the power spectrum components ``C_{\\ell}``,
starting from `` \\ell = 0 ``.
- `lmax::Integer`: the maximum ``ℓ`` coefficient, will default to `length(cl)-1` if not specified.
- `mmax::Integer`: the maximum ``m`` coefficient, will default to `lmax` if not specified.
- `rng::AbstractRNG` : (optional) the RNG to be used for generating the ``a_{\\ell m}``. It allows
to set the seed beforehand guaranteeing the reproducibility of the process.
"""
function synalm(cl::Vector{T}, lmax::Integer, mmax::Integer, rng::AbstractRNG) where {T <: Real}
    cl_size = length(cl)
    (cl_size - 1 >= lmax) || throw(DomainError(cl_size, "not enough C_l's to generate Alm"))
    alm = Alm{ComplexF64, Vector{ComplexF64}}(lmax, mmax)
    synalm!(cl, alm, rng)
    alm
end

synalm(cl::Vector{T}, lmax::Integer, mmax::Integer) where {T <: Real} =
    synalm(cl, lmax, mmax, Random.GLOBAL_RNG)

synalm(cl::Vector{T}, lmax::Integer, rng::AbstractRNG) where {T <: Real} =
    synalm(cl, lmax, lmax, rng)

synalm(cl::Vector{T}, lmax::Integer) where {T <: Real} =
    synalm(cl, lmax, lmax, Random.GLOBAL_RNG)

synalm(cl::Vector{T}, rng::AbstractRNG) where {T <: Real} =
    synalm(cl, length(cl) - 1, length(cl) - 1, rng)

synalm(cl::Vector{T}) where {T <: Real} =
    synalm(cl, length(cl) - 1, length(cl) - 1, Random.GLOBAL_RNG)

#########################################################################

"""
    synfast!(cl::Vector{T}, map::HealpixMap{T, RingOrder}, lmax::Integer, rng::AbstractRNG) where {T <: Real}
    synfast!(cl::Vector{T}, map::HealpixMap{T, RingOrder}, lmax::Integer) where {T <: Real}
    synfast!(cl::Vector{T}, map::HealpixMap{T, RingOrder}, rng::AbstractRNG) where {T <: Real}
    synfast!(cl::Vector{T}, map::HealpixMap{T, RingOrder}) where {T <: Real}

Generate a map from a given power spectra ``C_{\\ell}``. The result is saved into
the `HealpixMap` passed in input.

<<<<<<< HEAD
# Arguments:
=======
# ARGUMENTS
>>>>>>> 873b35eff1611e0a5f8921640a037354a4f82bb7
- `cl::AbstractVector{T}`: The array representing the power spectrum components ``C_{\\ell}``.
- `map::HealpixMap{T, RingOrder}`: the map that will contain the result.
- `lmax::Integer`: the maximum ``ℓ`` coefficient, will default to `length(cl)-1` if not specified.
- `rng::AbstractRNG` : (optional) the RNG to be used for generating the ``a_{\\ell m}``. It allows
to set the seed beforehand guaranteeing the reproducibility of the process.
"""
function synfast!(cl::Vector{T}, map::HealpixMap{T, RingOrder}, lmax::Integer, rng::AbstractRNG) where {T <: Real}
    cl_size = length(cl)
    (cl_size - 1 >= lmax) || throw(DomainError(cl_size, "not enough C_l's to generate Alm"))
    alm = Alm{ComplexF64, Vector{ComplexF64}}(lmax, lmax)
    synalm!(cl, alm)
    alm2map!(alm, map)
end

synfast!(cl::Vector{T}, map::HealpixMap{T, RingOrder}, lmax::Integer) where {T <: Real} =
    synfast!(cl, map, lmax, Random.GLOBAL_RNG)

synfast!(cl::Vector{T}, map::HealpixMap{T, RingOrder}, rng::AbstractRNG) where {T <: Real} =
    synfast!(cl, map, length(cl) - 1, rng)

synfast!(cl::Vector{T}, map::HealpixMap{T, RingOrder}) where {T <: Real} =
    synfast!(cl, map, length(cl) - 1, Random.GLOBAL_RNG)

#########################################################################

"""
    synfast(cl::Vector{T}, nside::Integer, lmax::Integer, rng::AbstractRNG) where {T <: Real}
    synfast(cl::Vector{T}, nside::Integer, lmax::Integer) where {T <: Real}
    synfast(cl::Vector{T}, nside::Integer, rng::AbstractRNG) where {T <: Real}
    synfast(cl::Vector{T}, nside::Integer) where {T <: Real}

Generate a `HealpixMap` with given Nside, from a given power spectra ``C_{\\ell}``.

<<<<<<< HEAD
# Arguments:
- `cl::AbstractVector{T}`: The array representing the power spectrum components ``C_{\\ell}``.
- `nside::Integer`: nside of the map that will contain the result.
- `lmax::Integer`: the maximum ``ℓ`` coefficient, will default to `length(cl)`-1 if not specified.
=======
# ARGUMENTS
- `cl::AbstractVector{T}`: The array representing the power spectrum components ``C_{\\ell}``.
- `nside::Integer`: nside of the map that will contain the result.
- `lmax::Integer`: the maximum ``ℓ`` coefficient, will default to `length(cl)-1` if not specified.
>>>>>>> 873b35eff1611e0a5f8921640a037354a4f82bb7
- `rng::AbstractRNG` : (optional) the RNG to be used for generating the ``a_{\\ell m}``. It allows
to set the seed beforehand guaranteeing the reproducibility of the process.
"""
function synfast(cl::Vector{T}, nside::Integer, lmax::Integer, rng::AbstractRNG) where {T <: Real}
    map=HealpixMap{T, RingOrder}(nside)
    synfast!(cl, map, lmax, rng)
    map
end

synfast(cl::Vector{T}, nside::Integer, lmax::Integer) where {T <: Real} =
    synfast(cl, nside, lmax, Random.GLOBAL_RNG)

synfast(cl::Vector{T}, nside::Integer, rng::AbstractRNG) where {T <: Real} =
    synfast(cl, nside, length(cl) - 1, rng)

synfast(cl::Vector{T}, nside::Integer) where {T <: Real} =
    synfast(cl, nside, length(cl) - 1, Random.GLOBAL_RNG)

#########################################################################

"""
    anafast(map::HealpixMap{Float64, RingOrder, AA}; lmax=nothing, mmax=nothing, niter::Integer = 3) where {T <: Real,AA <: AbstractArray{T,1}} -> Vector{Float64}
    anafast(map₁::HealpixMap{Float64, RingOrder, AA}, map₂::HealpixMap{Float64, RingOrder, AA}; lmax=nothing, mmax=nothing, niter::Integer = 3) where {T <: Real,AA <: AbstractArray{T,1}} -> Vector{Float64}

Computes the power spectrum of a Healpix map, or the cross-spectrum between two maps if `map2` is given.
No removal of monopole or dipole is performed. The input maps must be in ring-ordering.

<<<<<<< HEAD
# Arguments:
- `map₁::HealpixMap{Float64, RingOrder, AA}`: the spherical harmonic coefficients of the first field
- `map₂::HealpixMap{Float64, RingOrder, AA}`: the spherical harmonic coefficients of the second field

# Returns:
=======
# Arguments
- `map₁::HealpixMap{Float64, RingOrder, AA}`: the spherical harmonic coefficients of the first field
- `map₂::HealpixMap{Float64, RingOrder, AA}`: the spherical harmonic coefficients of the second field

# Returns
>>>>>>> 873b35eff1611e0a5f8921640a037354a4f82bb7
- `Array{T}` containing ``C_{\\ell}``, with the first element referring to ℓ=0.
"""
function anafast(map::HealpixMap{Float64, RingOrder, AA}; lmax=nothing, mmax=nothing, niter::Integer = 3) where {T <: Real,AA <: AbstractArray{T,1}}
    alm2cl(map2alm(map; lmax = lmax, mmax = mmax, niter=niter))
end

function anafast(
    map₁::HealpixMap{Float64, RingOrder, AA},
    map₂::HealpixMap{Float64, RingOrder, AA};
    lmax=nothing,
    mmax=nothing,
    niter::Integer = 3
    ) where {T <: Real,AA <: AbstractArray{T,1}}

    alm2cl(map2alm(map₁; lmax = lmax, mmax = mmax, niter=niter), map2alm(map₂; lmax = lmax, mmax = mmax, niter=niter))
end
