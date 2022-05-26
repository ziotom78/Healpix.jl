# -*- encoding: utf-8 -*-

import Healpix
using Test

const eps = 1e-12

@testset "Resolution" begin
    include("test_resolution.jl")
end

@testset "Math functions" begin
    include("test_math.jl")
end

@testset "Resolution-related functions" begin
    include("test_nside.jl")
end

@testset "Pixel functions" begin
    include("test_pixelfunctions.jl")
end

@testset "Polarized maps" begin
    include("test_polarizedmap.jl")
end

@testset "Interpolation functions" begin
    include("test_interp.jl")
end

@testset "XYF representation" begin
    include("test_xyf.jl")
end

@testset "Iterator interface for maps" begin
    include("test_mapiterators.jl")
end

@testset "Conformability" begin
    include("test_conformability.jl")
end

@testset "Mapmaking" begin
    include("test_mapmaking.jl")
end

@testset "Map I/O" begin
    include("test_mapio.jl")
end

@testset "Projections" begin
    include("test_projections.jl")
end

@testset "Alm coefficients" begin
    include("test_alm.jl")
end

@testset "Spherical Harmonic Transforms" begin
    include("test_sphtfunc.jl")
end

@testset "Query disc" begin
    include("test_querydisc.jl")
end

# Alm creation
