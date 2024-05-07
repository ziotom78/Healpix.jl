function cosdist_zphi(z1, phi1, z2, phi2)
    z1 * z2 + cos(phi1 - phi2) * sqrt((1 - z1^2) * (1 - z2^2))
end

function checkPixelRing(
    b1::Resolution,
    b2::Resolution,
    pix,
    nr,
    ipix1,
    fct,
    cz,
    cphi,
    cosrp2,
    cpix,
)
    (pix >= nr) && (pix -= nr)
    (pix < 0) && (pix += nr)

    pix += ipix1

    (pix == cpix) && return false

    (px, py, pf) = pix2xyfRing(b1, pix)

    for i in 0:(fct - 1)
        ox = fct * px
        oy = fct * py

        for (curx, cury) in (
            (ox + i, oy),
            (ox + fct - 1, oy + i),
            (ox + fct - 1 - i, oy + fct - 1),
            (ox, oy + fct - 1 - i),
        )
            (pz, pphi) = pix2zphiRing(b2, xyf2pixRing(b2, curx, cury, pf))
            (cosdist_zphi(pz, pphi, cz, cphi) > cosrp2) && (return false)
        end
    end

    true
end


@doc raw"""
    queryDiscRing(resol::Resolution, theta, phi, radius; fact=0)

Return a list of the indices of those pixels whose centers are closer
than `radius` to direction `(theta, phi)`. The three angles `radius`,
`theta`, and `phi` must be expressed in radians.

If `fact` is nonzero, it must be a positive integer; it requires to
carry the computation at a resolution `fact * nside`.

"""
function queryDiscRing(
    resol::Resolution,
    theta,
    phi,
    radius;
    fact=0,
)
    inclusive = (fact != 0)
    result = Int[]

    fct = 1

    if inclusive
        @assert ((1 << ORDER_MAX) / resol.nside) >= fact
        fct = fact
    end

    b2 = Resolution(fct * resol.nside)
    (rsmall, rbig) = if fct > 1
        (radius + max_pixrad(b2), radius + max_pixrad(resol))
    else
        value = inclusive ? (radius + max_pixrad(resol)) : radius
        (value, value)
    end

    (rsmall >= π) && (return 1:resol.numOfPixels)

    rbig = min(pi, rbig)
    (cosrsmall, cosrbig) = (cos(rsmall), cos(rbig))

    z0 = cos(theta)
    xa = 1 / sqrt((1 - z0) * (1 + z0))

    cpix = zphi2pixRing(resol, z0, phi)

    rlat1 = theta - rsmall
    zmax = cos(rlat1)
    irmin = ringAbove(resol, zmax) + 1


    if (rlat1 <= 0) && (irmin > 1)
        # The North Pole is within the disk
        ringinfo = getringinfo(resol, irmin - 1)
        append!(
            result,
            1:(ringinfo.firstPixIdx + ringinfo.numOfPixels - 1),
        )
    end

    if (fct > 1) && (rlat1 > 0)
        irmin = max(1, irmin - 1)
    end

    rlat2 = theta + rsmall
    zmin = cos(rlat2)
    irmax = ringabove(resol, zmin)

    if (fct > 1) && (rlat2 < π)
        irmax = min(resol.nsideTimesFour - 1, irmax + 1)
    end

    for iz in irmin:irmax
        z = ring2z(resol, iz)
        x = (cosrbig - z * z0) * xa
        ysq = 1 - z^2 - x^2
        dphi = if ysq < 0
            (fct == 1) ? 0 : π - 1e-15
        else
            atan(sqrt(ysq), x)
        end

        if dphi > 0
            ringinfo = getringinfo(resol, iz)
            ipix1 = ringinfo.firstPixIdx
            nr = ringinfo.numOfPixels
            shift = ringinfo.shifted ? 0.5 : 0.0

            ipix2 = ipix1 + nr - 1

            ip_lo = floor(Int, nr / 2π * (phi - dphi) - shift) + 1
            ip_hi = floor(Int, nr / 2π * (phi + dphi) - shift)

            if fct > 1
                while ((ip_lo <= ip_hi) &&
                       checkPixelRing(resol, b2, ip_lo, nr, ipix1, fct, z0, phi, cosrsmall, cpix))
                    ip_lo += 1
                end
                while ((ip_hi > ip_lo) &&
                       checkPixelRing(resol, b2, ip_hi, nr, ipix1, fct, z0, phi, cosrsmall, cpix))
                    ip_hi -= 1
                end
            end

            if ip_lo <= ip_hi
                if ip_hi >= nr
                    ip_lo -= nr
                    ip_hi -= nr
                end

                if ip_lo < 0
                    append!(result, ipix1:(ipix1 + ip_hi))
                    append!(result, (ipix1 + ip_lo + nr):ipix2)
                else
                    append!(result, (ipix1 + ip_lo):(ipix1 + ip_hi))
                end
            end
        end
    end

    if (rlat2 >= π) && (irmax + 1 < resol.nsideTimesFour)
        ringinfo = getringinfo(resol, irmax + 1)
        append!(result, ringinfo.firstPixIdx:resol.numOfPixels)
    end

    result
end


@doc raw"""
    queryStripRing(resol::Resolution, theta1, theta2; inclusive=true)

Return a range of the indices of pixels that overlap with
the colatitude range `[theta1, theta2]`. If `inclusive` is set to
`false`, only those pixels whose *centers* lie within the colatitude
range are returned.

This function assumes the RING scheme, because in this case, the
indexes of the pixels cover a range without gaps. Therefore, the
function returns a *range* instead of a *list*, as it is quicker
and occupies far less memory.

"""
function queryStripRing(resol::Resolution, theta1, theta2; inclusive=true)

    if theta1 > theta2
        throw(DomainError("Call to queryStripRing with theta1 > theta2 ($theta1 > $theta2)"))
    end

    maxring = numOfRings(resol)
    ring1 = max(1, 1 + ringAbove(resol, cos(theta1)))
    ring2 = min(maxring, ringAbove(resol, cos(theta2)))

    if inclusive
        ring1 = max(1, ring1 - 1)
        ring2 = min(maxring, ring2 + 1)
    end

    ringinfo = getringinfo(resol, ring1, full = false)
    (sp1, rp1) = (ringinfo.firstPixIdx, ringinfo.numOfPixels)
    getringinfo!(resol, ring2, ringinfo, full = false)
    (sp2, rp2) = (ringinfo.firstPixIdx, ringinfo.numOfPixels)

    (pix1, pix2) = (sp1, sp2 + rp2)

    (pix1 <= pix2) ? (pix1:(pix2 - 1)) : Int[]
end
