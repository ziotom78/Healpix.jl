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
        shifted_argument >>= 2
    end

    if shifted_argument > 0x00000001
        result |= 1
    end

    result
end
