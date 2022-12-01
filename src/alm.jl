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
#INDEXING

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
    (lmax >= 0) || throw(DomainError(lmax, "`lmax` is not positive or zero"))
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


""" each_ell(alm::Alm{Complex{T}}, m::Integer) where {T <: Number} -> Vector{Int}
    each_ell(alm::Alm{Complex{T}}, ms::AbstractArray{I, 1}) where {T <: Number, I <: Integer} -> Vector{Int}

    Returns an array of all the allowed ℓ values in `alm` for the given `m`.
"""
function each_ell(alm::Alm{Complex{T}}, m::Integer) where {T <: Number}
    (m <= alm.mmax) || throw(DomainError(m, "`m` is greater than mmax"))
    [l for l in m:alm.lmax]
end

function each_ell(alm::Alm{Complex{T}}, ms::AbstractArray{I, 1}) where {T <: Number, I <: Integer}
    reduce(vcat, [each_ell(alm, m) for m in ms])
end

""" each_ell_idx(alm::Alm{Complex{T}}, m::Integer) where {T <: Number} -> Vector{Int}
    each_ell_idx(alm::Alm{Complex{T}}, ms::AbstractArray{I, 1}) where {T <: Number, I <: Integer} -> Vector{Int}

    Returns an array of the indexes of the harmonic coefficients in `alm` corresponding
    to all the ℓ values for the given m value(s).
"""
function each_ell_idx(alm::Alm{Complex{T}}, m::Integer) where {T <: Number}
    (m <= alm.mmax) || throw(DomainError(m, "`m` is greater than mmax"))
    [i for i in almIndex(alm, m, m):almIndex(alm, alm.lmax, m)]
end

function each_ell_idx(alm::Alm{Complex{T}}, ms::AbstractArray{I, 1}) where {T <: Number, I <: Integer}
    reduce(vcat, [each_ell_idx(alm, m) for m in ms])
end

""" each_m(alm::Alm{Complex{T}}, l::Integer) where {T <: Number} -> Vector{Int}
    each_m(alm::Alm{Complex{T}}, ls::AbstractArray{I, 1}) where {T <: Number, I <: Integer} -> Vector{Int}

    Returns an array containing all the allowed m values in `alm` for the given ℓ value(s).
"""
function each_m(alm::Alm{Complex{T}}, l::Integer) where {T <: Number}
    (l <= alm.lmax) || throw(DomainError(l, "`l` is greater than lmax"))
    if l <= alm.mmax
        maxm = l
    else
        maxm = alm.mmax
    end
    [m for m in 0:maxm]
end

function each_m(alm::Alm{Complex{T}}, ls::AbstractArray{I, 1}) where {T <: Number, I <: Integer}
    reduce(vcat, [each_m(alm, l) for l in ls])
end

""" each_m_idx(alm::Alm{Complex{T}}, l::Integer) where {T <: Number} -> Vector{Int}
    each_m_idx(alm::Alm{Complex{T}}, ls::AbstractArray{I, 1}) where {T <: Number, I <: Integer} -> Vector{Int}

    Returns an array of the indexes of the harmonic coefficients in `alm` corresponding
    to all the allowed m values for the given ℓ value(s).
"""
function each_m_idx(alm::Alm{Complex{T}}, l::Integer) where {T <: Number}
    (l <= alm.lmax) || throw(DomainError(l, "`l` is greater than lmax"))
    [almIndex(alm, l, m) for m in each_m(alm, l)]
end

function each_m_idx(alm::Alm{Complex{T}}, ls::AbstractArray{I, 1}) where {T <: Number, I <: Integer}
    reduce(vcat, [each_m_idx(alm, l) for l in ls])
end

""" each_ell_m(alm::Alm{Complex{T}}) where {T <: Number} -> Vector{Int}

    Returns an array of tuples `(l, m)` of all the ℓ and m values of `alm` in
    m-major order (the same order as how the harmonic coefficients are stored in `Alm` objects).
"""
function each_ell_m(alm::Alm{Complex{T}}) where {T <: Number}
    ell_m = Vector(undef, numberOfAlms(alm.lmax, alm.mmax))
    i = 1
    for m in 0:alm.mmax
        for l in each_ell(alm, m)
            ell_m[i] = (l,m)
            i += 1
        end
    end
    ell_m
end

import Base: eachindex

""" eachindex(alm::Alm{Complex{T}}) where {T <: Number}

    Works as `eachindex(alm.alm)`.
"""
function eachindex(alm::Alm{Complex{T}}) where {T <: Number}
    eachindex(alm.alm)
end
############################################################################

"""
    writeAlmToFITS(f::CFITSIO.FITSFile, alm::Alm{Complex{T}}) where {T <: Number}
    writeAlmToFITS(fileName, alm::Alm{Complex{T}}; overwrite = true) where {T <: Number}

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
    alm2cl(alm::Alm{Complex{T}}) where {T <: Number} -> Vector{T}
    alm2cl(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number} -> Vector{T}

Compute ``C_{\\ell}`` from the spherical harmonic coefficients of one
or two fields.

# Arguments
- `alm₁::Alm{Complex{T}}`: the spherical harmonic coefficients of the first field
- `alm₂::Alm{Complex{T}}`: the spherical harmonic coefficients of the second field

# Returns
- `Array{T}` containing ``C_{\\ell}``, with the first element referring to ℓ=0.
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

###########################################################################

# eq 54, 55 of https://arxiv.org/abs/astro-ph/0008228
# these are asymptotic beams for σ² << 1

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
#ALM ALGEBRA

"""
    almxfl!(alm::Alm{Complex{T}}, fl::AbstractVector{T}) where {T <: Number}

Multiply IN-PLACE an a_ℓm by a vector b_ℓ representing an ℓ-dependent function.

# ARGUMENTS
- `alms::Alm{Complex{T}}`: The array representing the spherical harmonics coefficients
- `fl::AbstractVector{T}`: The array giving the factor f_ℓ by which to multiply a_ℓm

"""
function almxfl!(alm::Alm{Complex{T}}, fl::AbstractVector{T}) where {T <: Number}

    lmax = alm.lmax
    mmax = alm.mmax
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
            i = almIndex(alm, l, m)
            alm.alm[i] = alm.alm[i]*f
        end
    end
end

"""
    almxfl(alm::Alm{Complex{T}}, fl::AbstractVector{T}) where {T <: Number} -> Alm{T}

Multiply an a_ℓm by a vector b_ℓ representing an ℓ-dependent function, without changing
the a_ℓm passed in input.

# ARGUMENTS
- `alm::Alm{Complex{T}}`: The array representing the spherical harmonics coefficients
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


import Base: +, -, *, /, \
import LinearAlgebra: dot

""" +(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}

    Perform the element-wise sum in a_ℓm space.
    A new `Alm` object is returned.
"""
function +(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}
    (length(alm₁.alm) == length(alm₂.alm)) || throw(DomainError("Alms sizes not matching"))
    (alm₁.lmax == alm₂.lmax) || throw(DomainError("lmax's not matching"))
    (alm₁.mmax == alm₂.mmax) || throw(DomainError("mmax's not matching"))

    res_alm = Alm(alm₁.lmax, alm₁.mmax, Vector{Complex{T}}(undef, length(alm₁.alm)))

     @inbounds for i in eachindex(alm₁)
        res_alm.alm[i] = alm₁.alm[i] + alm₂.alm[i]
    end
    res_alm
end

""" -(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}

    Perform the element-wise subtraction in a_ℓm space.
    A new `Alm` object is returned.
"""
function -(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}
    (length(alm₁.alm) == length(alm₂.alm)) || throw(DomainError("Alms sizes not matching"))
    (alm₁.lmax == alm₂.lmax) || throw(DomainError("lmax's not matching"))
    (alm₁.mmax == alm₂.mmax) || throw(DomainError("mmax's not matching"))

    res_alm = Alm(alm₁.lmax, alm₁.mmax, Vector{Complex{T}}(undef, length(alm₁.alm)))

     @inbounds for i in eachindex(alm₁)
        res_alm.alm[i] = alm₁.alm[i] - alm₂.alm[i]
    end
    res_alm
end

""" *(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}

    Perform the element-wise product in a_ℓm space.
    A new `Alm` object is returned.
"""
function *(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}
    (length(alm₁.alm) == length(alm₂.alm)) || throw(DomainError("Alms sizes not matching"))
    (alm₁.lmax == alm₂.lmax) || throw(DomainError("lmax's not matching"))
    (alm₁.mmax == alm₂.mmax) || throw(DomainError("mmax's not matching"))
    lmax = alm₁.lmax

    #first part of the alm arrays, where for m=0 whe have alm^R_l,0 = alm^C_l,0,
    #and thus only real values are interesting (imag. should be 0)
    res_alm = Alm(lmax, alm₁.mmax, Vector{Complex{T}}(undef, length(alm₁.alm)))
    @inbounds for i in 1:lmax+1
        res_alm.alm[i] = real(alm₁.alm[i]) * real(alm₂.alm[i])
    end
    #we then compute the rest, where alm^R_l,m = √2 Re{alm^C_l,m}, alm^R_l,-m = √2Im{alm^C_l,m}
    @inbounds for i in lmax+2:length(alm₁.alm)
        res_alm.alm[i] = real(alm₁.alm[i]) * real(alm₂.alm[i]) + im * imag(alm₁.alm[i]) * imag(alm₂.alm[i])
    end
    res_alm
end

""" *(alm₁::Alm{Complex{T}}, fl::AbstractVector{T}) where {T <: Number}

    Perform the product of an `Alm` object by a function of ℓ in a_ℓm space.
    Note the order of the arguments: this consists in a shortcut of [`almxfl`](@ref),
    therefore a new `Alm` object is returned. Swap the arguments for an in-place
    version.
"""
function *(alm₁::Alm{Complex{T}}, fl::AbstractVector{T}) where {T <: Number}
    almxfl(alm₁, fl)
end

""" *(fl::AbstractVector{T}, alm₁::Alm{Complex{T}}) where {T <: Number}

    Perform the in-place product of an `Alm` object by a function of ℓ in a_ℓm space.
    Note the order of the arguments: this consists in a shortcut of [`almxfl!`](@ref),
    therefore the multiplication is performed IN PLACE.
"""
function *(fl::AbstractVector{T}, alm₁::Alm{Complex{T}}) where {T <: Number}
    almxfl!(alm₁, fl)
    alm₁
end

""" *(alm₁::Alm{Complex{T}}, fl::AbstractVector{T}) where {T <: Number}

    Perform the element-wise product of an `Alm` object by a constant in a_ℓm space.
    Note the order of the arguments: in this case a new `Alm` object is returned.
    Swap the arguments for an in-place version.
"""
function *(alm₁::Alm{Complex{T}}, c::Number) where {T <: Number}
    res_alm = Alm(alm₁.lmax, alm₁.mmax, Vector{Complex{T}}(undef, length(alm₁.alm)))

    @inbounds for i in eachindex(alm₁)
        res_alm.alm[i] = alm₁.alm[i] * c
    end
    res_alm
end

""" *(c::Number, alm₁::Alm{Complex{T}}) where {T <: Number}

    Perform the in-place element-wise product of an `Alm` object by a constant in a_ℓm space.
    Note the order of the arguments: in this case the product is performed IN PLACE.
"""
function *(c::Number, alm₁::Alm{Complex{T}}) where {T <: Number}
    @inbounds for i in eachindex(alm₁)
        alm₁.alm[i] *= c
    end
    alm₁
end

""" /(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}

    Perform an element-wise division in a_ℓm space between two `Alm`s.
    A new `Alm` object is returned.
"""
function /(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}
    (length(alm₁.alm) == length(alm₂.alm)) || throw(DomainError("Alms sizes not matching"))
    (alm₁.lmax == alm₂.lmax) || throw(DomainError("lmax's not matching"))
    (alm₁.mmax == alm₂.mmax) || throw(DomainError("mmax's not matching"))
    lmax = alm₁.lmax

    #first part of the alm arrays, where for m=0 whe have alm^R_l,0 = alm^C_l,0,
    #and thus only real values are interesting (imag. should be 0)
    res_alm = Alm(lmax, alm₁.mmax, Vector{Complex{T}}(undef, length(alm₁.alm)))
    @inbounds for i in 1:lmax+1
        res_alm.alm[i] = real(alm₁.alm[i]) / real(alm₂.alm[i])
    end
    #we then compute the rest, where alm^R_l,m = √2 Re{alm^C_l,m}, alm^R_l,-m = √2Im{alm^C_l,m}
    @inbounds for i in lmax+2:length(alm₁.alm)
        res_alm.alm[i] = real(alm₁.alm[i]) / real(alm₂.alm[i]) + im * imag(alm₁.alm[i]) / imag(alm₂.alm[i])
    end
    res_alm
end

""" /(alm₁::Alm{Complex{T}}, fl::AbstractVector{T}) where {T <: Number}

    Perform an element-wise division by a function of ℓ in a_ℓm space.
    A new `Alm` object is returned.
"""
/(alm₁::Alm{Complex{T}}, fl::AbstractVector{T}) where {T <: Number} = almxfl(alm₁, 1. ./ fl)

""" /(alm₁::Alm{Complex{T}}, c::Number) where {T <: Number}

    Perform an element-wise division by a constant in a_ℓm space.
    A new `Alm` object is returned.
"""
function /(alm₁::Alm{Complex{T}}, c::Number) where {T <: Number}
    res_alm = Alm(alm₁.lmax, alm₁.mmax, Vector{Complex{T}}(undef, length(alm₁.alm)))

    @inbounds for i in eachindex(alm₁)
        res_alm.alm[i] = alm₁.alm[i] / c
    end
    res_alm
end

"""
    Perform an IN-PLACE element-wise division by a function of ℓ in a_ℓm space.
"""
\(fl::AbstractVector{T}, alm₁::Alm{Complex{T}}) where {T <: Number} = almxfl!(alm₁, 1. ./ fl)

"""
    Perform an IN-PLACE element-wise division by a constant in a_ℓm space.
"""
function \(c::Number, alm₁::Alm{Complex{T}}) where {T <: Number}

    @inbounds for i in eachindex(alm₁)
        alm₁.alm[i] /= c
    end
    alm₁
end

""" dot(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}

    Implements the dot product in a_ℓm space.
    The two imput alms must have matching size, lmax and mmax.
    A new `Alm` object is returned.
"""
function dot(alm₁::Alm{Complex{T}}, alm₂::Alm{Complex{T}}) where {T <: Number}
    (length(alm₁.alm) == length(alm₂.alm)) || throw(DomainError("Alms sizes not matching"))
    (alm₁.lmax == alm₂.lmax) || throw(DomainError("lmax's not matching"))
    (alm₁.mmax == alm₂.mmax) || throw(DomainError("mmax's not matching"))
    lmax = alm₁.lmax
    res = 0.0

    #we first compute the l≠0, where alm^R_l,m = √2 Re{alm^C_l,m}, alm^R_l,-m = √2Im{alm^C_l,m}
    @inbounds for i in lmax+2:length(alm₁.alm)
        res += real(alm₁.alm[i]) * real(alm₂.alm[i]) + imag(alm₁.alm[i]) * imag(alm₂.alm[i])
    end
    res *= sqrt(2)

    #then the first part of the alm arrays, where for m=0 whe have alm^R_l,0 = alm^C_l,0,
    #and thus only real values are interesting (imag. should be 0)
    @inbounds for i in 1:lmax+1
        res += real(alm₁.alm[i]) * real(alm₂.alm[i])
    end
    res
end
