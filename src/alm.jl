# Definition of the composite type Alm

"An array of a_ℓm numbers."
mutable struct Alm{T <: Number}
    alm::Array{T, 1}
    lmax::Int
    mmax::Int
    tval::Int

    Alm{T}(lmax, mmax) where {T <: Number} = new(zeros(T, numberOfAlms(lmax, mmax)),
                                                 lmax,
                                                 mmax,
                                                 2lmax + 1)

    function Alm{T}(lmax, mmax, arr::Array{T, 1}) where {T <: Number}
        (numberOfAlms(lmax, mmax) == length(arr)) || throw(DomainError())

        new{T}(arr, lmax, mmax, 2lmax + 1)
    end
end
                                            
################################################################################

"""
    numberOfAlms(lmax::Integer, mmax::Integer) -> Integer
    numberOfAlms(lmax::Integer) -> Integer

Return the size of the array of complex numbers needed to store the
a_lm coefficients in the range of ℓ and m specified by `lmax` and
`mmax`. If `mmax` is not specified, it is assumed to be equal to
`lmax`. If `lmax` and `mmax` are inconsistent or negative, a
`DomainError` exception is thrown.
"""
function numberOfAlms(lmax, mmax)
    (lmax >= 0) || throw(DomainError(lmax, "`lmax` is not positive or zero"))
    (mmax >= 0) || throw(DomainError(mmax, "`mmax` is not positive or zero"))
    (0 ≤ mmax ≤ lmax) || throw(DomainError((lmax, mmax), "`lmax` and `mmax` are inconsistent"))

    div((mmax + 1) * (mmax + 2), 2) + (mmax + 1) * (lmax - mmax)
end

numberOfAlms(lmax) = numberOfAlms(lmax, lmax)

shr(x, y) = x >> y
shr(x::Array{T}, y) where {T} = [a >> y for a in x]

almIndexL0(alm::Alm{T}, m) where {T} = shr((m .* (alm.tval .- m)), 1) .+ 1
almIndex(alm::Alm{T}, l, m) where {T} = almIndexL0(alm, m) .+ l

################################################################################

"""
    readAlmFromFITS{T <: Complex}(f::FITSIO.FITSFile, t::Type{T}) -> Alm{T}
    readAlmFromFITS{T <: Complex}(fileName::String, t::Type{T}) -> Alm{T}

Read a set of a_ℓm coefficients from a FITS file. If the code fails,
FITSIO will raise an exception. (Refer to the FITSIO library for more
information.)
"""
function readAlmFromFITS(f::FITSIO.FITSFile,
                         t::Type{T}) where {T <: Complex}
    numOfRows = FITSIO.fits_get_num_rows(f)

    idx = Array{Int64}(undef, numOfRows)
    almReal = Array{Float64}(undef, numOfRows)
    almImag = Array{Float64}(undef, numOfRows)

    FITSIO.fits_read_col(f, 1, 1, 1, idx)
    FITSIO.fits_read_col(f, 2, 1, 1, almReal)
    FITSIO.fits_read_col(f, 3, 1, 1, almImag)

    l = floor.(Int64, sqrt.(idx .- 1))
    m = idx .- l.^2 .- l .- 1
    if count(x -> x < 0, m) > 0
        throw(DomainError())
    end

    result = Alm{T}(maximum(l), maximum(m))
    i = almIndex(result, l, m)
    result.alm = complex.(almReal[i], almImag[i])
end

function readAlmFromFITS(fileName,
                         t::Type{T}) where {T <: Complex}
    f = FITSIO.fits_open_table(fileName)
    try
        result = readAlmFromFITS(f, t)
        return result
    finally
        FITSIO.fits_close_file(f)
    end
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
            cl[l+1] += 2 * real(alm₁.alm[index] * conj(alm₂.alm[index]))
        end
        index0 = almIndex(alm₁, l, 0)
        cl[l+1] += real(alm₁.alm[index0] * conj(alm₂.alm[index0]))
        cl[l+1] = cl[l+1] / (2 * l + 1)
    end
    return cl
end
alm2cl(alm::Alm{Complex{T}}) where {T <: Number} = alm2cl(alm, alm)

