# A few general-purpose mathematical functions

function ilog2(argument)

    local result = 0
    local shifted_argument = argument

    while shifted_argument > 0x0000FFFF
        result += 16
        shifted_argument >>= 16
    end

    if shifted_argument > 0x000000FF
        result |= 8
        shifted_argument >>= 8
    end

    if shifted_argument > 0x0000000F
        result |= 4
        shifted_argument >>= 4
    end

    if shifted_argument > 0x00000003
        result |= 2
        shifted_argument>>=2
    end

    if shifted_argument > 0x00000001
        result |= 1
    end

    result
end

"""
    normalizeAngle(x)

Return the same angle as the argument, but in the range [0, 2π). Note that this
is slightly different from mod2pi, as the latter returns a value in the range
[0, 2π].
"""
function normalizeAngle(x)
    while x >= 2π
        x -= 2π
    end

    while x < 0
        x += 2π
    end

    x
end
