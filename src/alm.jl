# Definition of the composite type Alm

"""An array of harmonic coefficients (a_ℓm).

The type `T` is used for the value of each harmonic coefficient, and
it must be a `Number` (one should however only use complex types for
this). The type `AA` is used to store the array of coefficients; a
typical choice is `Vector`.

A `Alm` type contains the following fields:

- `alm`: the array of harmonic coefficients
- `lmax`: the maximum value for ``ℓ``
- `mmax`: the maximum value for ``m``
- `tval`: maximum number of ``m`` coefficients for the maximum ``ℓ``

The ``a_{\\ell m}`` are stored by ``m``: if ``\\ell_{max}`` is 16, the first 16 elements
are ``m=0``, ``\\ell=0-16``, then the following 15 elements are ``m=1``, ``\\ell=1-16``,
then ``m=2``, ``\\ell=2-16`` and so on until the last element, the 153th, is ``m=16``, ``\\ell=16``.
"""
mutable struct Alm{T <: Number,AA <: AbstractArray{T,1}}
    alm::AA
    lmax::Int
    mmax::Int
    tval::Int

    Alm{T,AA}(lmax, mmax) where {T <: Number,AA <: AbstractArray{T,1}} =
        new{T,AA}(zeros(T, numberOfAlms(lmax, mmax)), lmax, mmax, 2lmax + 1)

    function Alm{T,AA}(lmax, mmax, arr::AA) where {T <: Number,AA <: AbstractArray{T,1}}
        (numberOfAlms(lmax, mmax) == length(arr)) || throw(DomainError())

        new{T,AA}(arr, lmax, mmax, 2lmax + 1)
    end
end

Alm{T}(lmax, mmax) where {T <: Number} = Alm{T,Array{T,1}}(lmax, mmax)
Alm(lmax, mmax) = Alm{ComplexF64}(lmax, mmax)
Alm(lmax, mmax, arr::AA) where {T <: Number,AA <: AbstractArray{T,1}} =
    Alm{T,AA}(lmax, mmax, arr)

################################################################################

"""
    numberOfAlms(lmax::Integer, mmax::Integer) -> Integer
    numberOfAlms(lmax::Integer) -> Integer

Return the size of the array of complex numbers needed to store the
a_ℓm coefficients in the range of ℓ and m specified by `lmax` and
`mmax`. If `mmax` is not specified, it is assumed to be equal to
`lmax`. If `lmax` and `mmax` are inconsistent or negative, a
`DomainError` exception is thrown.
"""
function numberOfAlms(lmax, mmax)
    (lmax >= 0) || throw(DomainError(lmax, "`lmax` is not positive or zero"))
    (mmax >= 0) || throw(DomainError(mmax, "`mmax` is not positive or zero"))
    (0 ≤ mmax ≤ lmax) ||
        throw(DomainError((lmax, mmax), "`lmax` and `mmax` are inconsistent"))

    div((mmax + 1) * (mmax + 2), 2) + (mmax + 1) * (lmax - mmax)
end

numberOfAlms(lmax) = numberOfAlms(lmax, lmax)

shr(x, y) = x >> y
shr(x::Array{T}, y) where {T} = [a >> y for a in x]

almIndexL0(alm::Alm{T}, m) where {T} = shr((m .* (alm.tval .- m)), 1) .+ 1
almIndex(alm::Alm{T}, l, m) where {T} = almIndexL0(alm, m) .+ l

################################################################################

"""
    readAlmFromFITS{T <: Complex}(f::CFITSIO.FITSFile, t::Type{T}) -> Alm{T}
    readAlmFromFITS{T <: Complex}(fileName::String, t::Type{T}) -> Alm{T}

Read a set of a_ℓm coefficients from a FITS file. If the code fails,
CFITSIO will raise an exception. (Refer to the CFITSIO library for more
information.)
"""
function readAlmFromFITS(f::CFITSIO.FITSFile, t::Type{T}) where {T <: Complex}
    numOfRows = CFITSIO.fits_get_num_rows(f)

    idx = Array{Int64}(undef, numOfRows)
    almReal = Array{Float64}(undef, numOfRows)
    almImag = Array{Float64}(undef, numOfRows)

    CFITSIO.fits_read_col(f, 1, 1, 1, idx)
    CFITSIO.fits_read_col(f, 2, 1, 1, almReal)
    CFITSIO.fits_read_col(f, 3, 1, 1, almImag)

    l = floor.(Int64, sqrt.(idx .- 1))
    m = idx .- l.^2 .- l .- 1
    if count(x -> x < 0, m) > 0
        throw(DomainError())
    end

    result = Alm{T}(maximum(l), maximum(m))
    i = almIndex(result, l, m)
    result.alm = complex.(almReal[i], almImag[i])
    result
end

function readAlmFromFITS(fileName, t::Type{T}) where {T <: Complex}
    f = CFITSIO.fits_open_table(fileName)
    try
        result = readAlmFromFITS(f, t)
        result
    finally
        CFITSIO.fits_close_file(f)
end
end

###############################################################

"""
    almExplicitIndex(lmax) -> Vector{Int}
    almExplicitIndex(lmax, mmax) -> Vector{Int}
    almExplicitIndex(alm::Alm{T}) where {T} -> Vector{Int}

Compute the explicit index scheme, i.e. ``\\mathrm{index} = \\ell^2 + \\ell + m + 1``
up to a certain ``ℓ`` and ``m`` if specified, or taken from the `Alm` passed.
If not passed, `mmax` is defaulted to `lmax`. If `lmax` and `mmax` are inconsistent
or negative, a `DomainError` exception is thrown.

"""

function almExplicitIndex(lmax, mmax)
    lmax >= 0) || throw(DomainError(lmax, "`lmax` is not positive or zero"))
    (mmax >= 0) || throw(DomainError(mmax, "`mmax` is not positive or zero"))
    (0 ≤ mmax ≤ lmax) ||
        throw(DomainError((lmax, mmax), "`lmax` and `mmax` are inconsistent"))

    # Step 1: count the number of elements in the output
    count = 0
    for m = 0:mmax
        (lmax ≥ m) && (count += (lmax - m) + 1)
    end

    # Step 2: allocate the vector for the output in one batch
    idx = Vector{Int}(undef, count)

    # Step 3: fill the output
    i = 1
    for m = 0:mmax
        for l = m:lmax
            idx[i] = l^2 + l + m + 1
            i += 1
        end
    end

    idx
end

almExplicitIndex(lmax) = almExplicitIndex(lmax, lmax)

almExplicitIndex(alm::Alm{T}) where {T} = almExplicitIndex(alm.lmax, alm.mmax)


############################################################################

"""
    writeAlmToFITS{T <: Complex}(f::CFITSIO.FITSFile, alm::Alm{Complex{T}})
    writeAlmToFITS{T <: Complex}(fileName::String, alm::Alm{Complex{T}})

Write a set of a_ℓm coefficients into a FITS file. If the code fails,
CFITSIO will raise an exception. (Refer to the CFITSIO library for more
information.)
In the fits file the alms are written with explicit index scheme,
``\\mathrm{index} = \\ell^2 + \\ell + m + 1``, possibly out of order (check `almExplicitIndex`).
"""
function writeAlmToFITS(f::CFITSIO.FITSFile, alm::Alm{Complex{T}}) where {T <: Number}

    idx = almExplicitIndex(alm)
    almReal = real(alm.alm)
    almImag = imag(alm.alm)

    CFITSIO.fits_write_col(f, 1, 1, 1, idx)
    CFITSIO.fits_write_col(f, 2, 1, 1, almReal)
    CFITSIO.fits_write_col(f, 3, 1, 1, almImag)

end

function writeAlmToFITS(fileName, alm::Alm{Complex{T}}; overwrite = true) where {T <: Number}
    if overwrite
        f = CFITSIO.fits_clobber_file(fileName)
    else
        f = CFITSIO.fits_create_file(fileName)
    end
    try
        CFITSIO.fits_create_binary_tbl(f, 0, [("Index", "1I", ""),("Re[alm]", "1D", ""),("Im[alm]", "1D", "")], "alm")
        writeAlmToFITS(f, alm)
    finally
        CFITSIO.fits_close_file(f)
end
end

########################################################################

"""
    alm2cl(alm::Alm{Complex{T}}) where {T <: Number}
    alm2cl(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}

Compute ``C_{\\ell}`` from the spherical harmonic coefficients of one
or two fields.

# Arguments
- `alm₁::Alm{Complex{T}}`: the spherical harmonic coefficients of the first field
- `alm₂::Alm{Complex{T}}`: the spherical harmonic coefficients of the second field

# Returns
- `Array{T}` containing C_{\\ell}, with the first element referring to ℓ=0.
"""
function alm2cl(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}
    (alm₁.lmax != alm₂.lmax) && throw(ArgumentError("Alm lmax do not match."))
    (alm₁.mmax != alm₂.mmax) && throw(ArgumentError("Alm mmax do not match."))
    (alm₁.mmax < alm₂.lmax) && throw(ArgumentError("Alm mmax < lmax."))

    lmax = alm₁.lmax
    cl = zeros(T, lmax + 1)
    for l = 0:lmax
        for m = 1:l
            index = almIndex(alm₁, l, m)
            cl[l + 1] += 2 * real(alm₁.alm[index] * conj(alm₂.alm[index]))
        end
        index0 = almIndex(alm₁, l, 0)
        cl[l + 1] += real(alm₁.alm[index0] * conj(alm₂.alm[index0]))
        cl[l + 1] = cl[l + 1] / (2 * l + 1)
    end
    cl
end

alm2cl(alm::Alm{Complex{T}}) where {T <: Number} = alm2cl(alm, alm)

# eq 54, 55 of https://arxiv.org/abs/astro-ph/0008228
# these are asymptotic beams for σ² << 1

###########################################################################

"""
    gaussbeam(fwhm::T, lmax::Int; pol=false) where T

Compute the Gaussian beam window function ``B_{\\ell}`` given the FWHM of the beam in radians, where
``C_{\\ell, \\mathrm{measured}} = B_{\\ell}^2 C_{\\ell}``. This beam is valid in the limit of
``\\sigma^2 \\ll 0``, which is the case for all high-resolution CMB experiments.

# Arguments
- `fwhm::T`: FWHM of the Gaussian beam in radians
- `lmax::Int`: maximum multipole ℓ
- `pol=false`: if false, returns the spin-0 beam for i.e. intensity. if true, returns the spin-2 beam

# Returns
- `Array{T,1}` containing ``B_{\\ell}``, with the first element referring to ℓ=0.
"""
function gaussbeam(fwhm::T, lmax::Int; pol=false) where T
    Bl = Array{T,1}(undef, lmax+1)
    fwhm²_to_σ² = 1 / (8log(T(2)))  # constant
    σ² = fwhm²_to_σ² * fwhm^2

    if pol
        for l = 0:lmax
            Bl[l+1] = exp(-(l * (l+1) - 4) * σ² / 2)
        end
    else
        for l = 0:lmax
            Bl[l+1] = exp(-l * (l+1) * σ² / 2)
        end
    end
    return Bl
end

###########################################################################

"""
    almxfl!(alms::Alm{Complex{T}}, fl::AbstractVector{T}) where {T <: Number}

Multiply IN-PLACE an a_ℓm by a vector b_ℓ representing an ℓ-dependent function.

# ARGUMENTS
- `alms::Alm{Complex{T}}`: The array representing the spherical harmonics coefficients
- `fl::AbstractVector{T}`: The array giving the factor f_ℓ by which to multiply a_ℓm

"""

function almxfl!(alms::Alm{Complex{T}}, fl::AbstractVector{T}) where {T <: Number}

    lmax = alms.lmax
    mmax = alms.mmax
    fl_size = length(fl)

    for l = 0:lmax
        if l < fl_size
            f = fl[l + 1]
        else
            f = 0
        end
        if l <= mmax
            maxm = l
        else
            maxm = mmax
        end

        for m = 0:maxm
            i = almIndex(alms, l, m)
            alms.alm[i] = alms.alm[i]*f
        end
    end
end

"""
    almxfl(alms::Alm{Complex{T}}, fl::AbstractVector{T}) where {T <: Number}

Multiply an a_ℓm by a vector b_ℓ representing an ℓ-dependent function, without changing
the a_ℓm passed in input.

# ARGUMENTS
- `alms::Alm{Complex{T}}`: The array representing the spherical harmonics coefficients
- `fl::AbstractVector{T}`: The array giving the factor f_ℓ by which to multiply a_ℓm

#RETURNS
- `Alm{Complex{T}}`: The result of a_ℓm * f_ℓ.
"""
function almxfl(alm::Alm{Complex{T}}, fl::AbstractVector{T}) where {T <: Number}
    lmax = alm.lmax
    mmax = alm.mmax
    alm_new = deepcopy(alm)
    almxfl!(alm_new, fl)
    alm_new
end
