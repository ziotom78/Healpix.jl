var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Healpix\nend"
},

{
    "location": "#Healpix.jl:-an-implementation-of-the-Healpix-tessellation-scheme-in-Julia-1",
    "page": "Introduction",
    "title": "Healpix.jl: an implementation of the Healpix tessellation scheme in Julia",
    "category": "section",
    "text": "This is the documentation of the Healpix.jl package, an implementation of the Healpix spherical tessellation scheme written entirely in Julia.This library is a work-in-progress: if you want something with more functionality, have a look at Libhealpix.jl, as it wraps the Healpix C++ library. This package has the main purpose of providing a Julia-only solution, so that it can easily be used on platforms not supported by the Healpix C++ library (e.g., Windows).This library implements algorithms for converting directions into pixel indices and vice versa. It supports both RING and NESTED schemes, and it employs Julia\'s powerful type system to avoid mistaking one scheme in place of the other."
},

{
    "location": "#Documentation-1",
    "page": "Introduction",
    "title": "Documentation",
    "category": "section",
    "text": "The documentation was built using Documenter.jl.println(\"Documentation built $(now()) with Julia $(VERSION).\") # hide"
},

{
    "location": "#Index-1",
    "page": "Introduction",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "resolutions/#",
    "page": "Working with resolutions",
    "title": "Working with resolutions",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Healpix\nend"
},

{
    "location": "resolutions/#Healpix.Resolution",
    "page": "Working with resolutions",
    "title": "Healpix.Resolution",
    "category": "type",
    "text": "struct Resolution\n\nResolution objects are needed to perform a number of pixel-related functions, e.g., convert a direction into a pixel number and vice versa.\n\nThe fields of a Resolution object are the following:\n\nnside: the NSIDE parameter\nnsideTimesTwo: 2 * NSIDE\nnsideTimesFour: 4 * NSIDE\nnumOfPixels: number of pixels in the map\norder: order of the map\npixelsPerFace: number of pixels in each Healpix face\nncap\nfact2\nfact1\n\n\n\n\n\n"
},

{
    "location": "resolutions/#Healpix.Resolution-Tuple{Integer}",
    "page": "Working with resolutions",
    "title": "Healpix.Resolution",
    "category": "method",
    "text": "Resolution(nside) -> Resolution\n\nCreate a Resolution object, given a value for NSIDE.\n\n\n\n\n\n"
},

{
    "location": "resolutions/#Healpix.nsideok-Tuple{Integer}",
    "page": "Working with resolutions",
    "title": "Healpix.nsideok",
    "category": "method",
    "text": "nsideok(nside::Integer) -> Bool\n\nCheck whether nside is a valid NSIDE parameter.\n\n\n\n\n\n"
},

{
    "location": "resolutions/#Healpix.nside2npix-Tuple{Integer}",
    "page": "Working with resolutions",
    "title": "Healpix.nside2npix",
    "category": "method",
    "text": "nside2npix(nside::Integer) -> Integer\n\nReturn the number of pixels for a Healpix map with the specified NSIDE value. If NSIDE is not an integer power of two, the function throws a DomainError exception.\n\n\n\n\n\n"
},

{
    "location": "resolutions/#Healpix.npix2nside-Tuple{Integer}",
    "page": "Working with resolutions",
    "title": "Healpix.npix2nside",
    "category": "method",
    "text": "npix2nside(npix::Integer) -> Integer\n\nGiven the number of pixels in a Healpix map, return the NSIDE resolution parameter. If the number is invalid, throw a DomainError exception.\n\n\n\n\n\n"
},

{
    "location": "resolutions/#Healpix.nside2pixarea-Tuple{Integer}",
    "page": "Working with resolutions",
    "title": "Healpix.nside2pixarea",
    "category": "method",
    "text": "nside2pixarea(nside::Integer) -> Real\n\nReturn the solid angle of a pixel in a map with the specified NSIDE parameter. The result is expressed in steradians.\n\n\n\n\n\n"
},

{
    "location": "resolutions/#Healpix.nside2resol-Tuple{Integer}",
    "page": "Working with resolutions",
    "title": "Healpix.nside2resol",
    "category": "method",
    "text": "nside2resol(nside::Integer) -> Real\n\nReturn the approximate resolution of a map with the specified NSIDE. The resolution is expressed in radians, and it is the square root of the pixel size.\n\n\n\n\n\n"
},

{
    "location": "resolutions/#Working-with-resolutions-1",
    "page": "Working with resolutions",
    "title": "Working with resolutions",
    "category": "section",
    "text": "A Healpix tessellation is parametrized by a number, called NSIDE, which must be a positive power of 2. It is related to the number of pixels N in the maps by the simple equation N = 12 mathrmNSIDE^2, and it is therefore related to the resolution of the pixelization. Any function working on a Healpix tessellation needs to receive the value of NSIDE. Healpix.jl provides a wrapper around this parameter, the Resolution type, which internally keeps a number of precomputed coefficients to accelerate calculations.The following example prints a table containing details about a few Healpix resolutions:using Healpix # hide\nusing Printf\n\n@printf(\"%-6s\\t%-12s\\t%-12s\\t%-12s\\n\",\n        \"NSIDE\",\n        \"#pix\",\n        \"#pix per face\",\n        \"solid angle\")\nfor poweroftwo in [0, 1, 2, 3, 4, 5]\n    res = Resolution(2 ^ poweroftwo)\n    @printf(\"%6d\\t%12d\\t%12d\\t%12.4f\\n\",\n            res.nside,\n            res.numOfPixels,\n            res.pixelsPerFace,\n            4π / res.numOfPixels)\nendResolution\nResolution(nside::Integer)\nnsideok(nside::Integer)\nnside2npix(nside::Integer)\nnpix2nside(npix::Integer)\nnside2pixarea(nside::Integer)\nnside2resol(nside::Integer)"
},

{
    "location": "pixelfunc/#",
    "page": "Pixel functions",
    "title": "Pixel functions",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Healpix\nend"
},

{
    "location": "pixelfunc/#Pixel-functions-1",
    "page": "Pixel functions",
    "title": "Pixel functions",
    "category": "section",
    "text": "In this section we document the functions that convert from a direction in the sky into a pixel index, and vice versa.First of all, Healpix.jl implements the most basic functions to convert between spherical and Cartesian coordinates. Note that Healpix uses co-latitude instead of latitude:using Healpix # hide\nang2vec(0.0, 0.0)\nvec2ang(0.0, 0.0, 1.0)More interesting functions return the index of the pixel on a Healpix-tessellated sphere. For these functions to work, you have to provide a Resolution object:res = Resolution(16)\nang2pixRing(res, π/2, 0)\nang2pixNest(res, π/2, 0)"
},

{
    "location": "pixelfunc/#Healpix.RingInfo",
    "page": "Pixel functions",
    "title": "Healpix.RingInfo",
    "category": "type",
    "text": "RingInfo\n\nInformation about a ring of pixels, i.e., the set of pixels on a Healpix map which have the same colatitude. The type is \"mutable\", so that one object can begin reused many times without further memory allocations.\n\nThe list of fields defined in this structure is the following:\n\nring: an integer index, running from \nfirstPixIdx: index of the first pixel (using the RING scheme) belonging to this ring\nnumOfPixels: number of consecutive pixels within the ring\ncolatitude_rad: value of the colatitude for this ring (in radians)\nshifted: Boolean flag; it is true if the longitude of the first pixel in the ring is not zero.\n\nReferences\n\nSee also getringinfo! and getringinfo.\n\nExample\n\nimport Healpix\nres = Healpix.Resolution(256)\n\n# Show information about ring #10\nprint(getringinfo(res, 10))\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.getringinfo",
    "page": "Pixel functions",
    "title": "Healpix.getringinfo",
    "category": "function",
    "text": "getringinfo(resol::Resolution, ring; kwargs...) :: RingInfo\n\nReturn a RingInfo structure containing information about the specified ring. For the list of accepted keyword arguments, see getringinfo!.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.getringinfo!",
    "page": "Pixel functions",
    "title": "Healpix.getringinfo!",
    "category": "function",
    "text": "getringinfo!(resol::Resolution, ring, ringinfo::RingInfo; full=true) :: RingInfo\n\nFill the RingInfo structure with information about the specified ring. If full is false, the field colatitude_rad (the most expensive in terms of computation) is set to NaN.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.getinterpolRing",
    "page": "Pixel functions",
    "title": "Healpix.getinterpolRing",
    "category": "function",
    "text": "getinterpolRing(resol::Resolution, θ, ϕ) :: (Array{Int,1}, Array{Float64, 1})\n\nReturn the indices and the weights of the four neighbour pixels for the given direction (θ, ϕ) in a map with the specified resolution.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Ring-functions-1",
    "page": "Pixel functions",
    "title": "Ring functions",
    "category": "section",
    "text": "The Healpix projection has the advantage of storing pixels along iso-latitude rings; this allows to implement efficient spherical-transform functions. Healpix.jl provides a number of functions that manage rings. Many of them use the RingInfo structure, which encodes details about a ring.RingInfo\ngetringinfo\ngetringinfo!\ngetinterpolRing"
},

{
    "location": "pixelfunc/#Healpix.ang2vec-Tuple{Any,Any}",
    "page": "Pixel functions",
    "title": "Healpix.ang2vec",
    "category": "method",
    "text": "ang2vec(theta, phi) -> Array{Float64}\n\nGiven a direction in the sky with colatitude theta and longitude phi (in radians), return an array of 3 elements containing the x, y, and z components of the one-length vector pointing to that direction.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.vec2ang-Tuple{Any,Any,Any}",
    "page": "Pixel functions",
    "title": "Healpix.vec2ang",
    "category": "method",
    "text": "vec2ang(x, y, z) -> (Number, Number)\n\nGiven a vector (not necessarily normalized) whose Cartesian components are x, y, and z, return a pair (theta, phi) containing the colatitude theta and the longitude phi (in radians) of the direction in the sky the vector is pointing at.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.ang2pixNest-Tuple{Resolution,Any,Any}",
    "page": "Pixel functions",
    "title": "Healpix.ang2pixNest",
    "category": "method",
    "text": "ang2pixNest(resol::Resolution, theta, phi) -> Integer\n\nReturn the index of the pixel which contains the point with coordinates (theta, the colatitude, and phi, the longitude), in radians, for a Healpix map with pixels in nested order. Note that pixel indexes are 1-based (this is Julia)!\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.ang2pixRing-Tuple{Resolution,Any,Any}",
    "page": "Pixel functions",
    "title": "Healpix.ang2pixRing",
    "category": "method",
    "text": "ang2pixRing(resol::Resolution, theta, phi) -> Integer\n\nReturn the index of the pixel which contains the point with coordinates (theta, the colatitude, and phi, the longitude), in radians, for a Healpix map with pixels in ring order. Note that pixel indexes are 1-based (this is Julia)!\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.pix2angNest-Tuple{Resolution,Any}",
    "page": "Pixel functions",
    "title": "Healpix.pix2angNest",
    "category": "method",
    "text": "pix2angNest(resol::Resolution, pixel) -> (Float64, Float64)\n\nGiven the (1-based) index of a pixel in a Healpix map in nested order, return a pair containing the (colatitude, longitude) angles corresponding to its center, both expressed in radians.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.pix2angRing-Tuple{Resolution,Any}",
    "page": "Pixel functions",
    "title": "Healpix.pix2angRing",
    "category": "method",
    "text": "pix2angRing(resol::Resolution, pixel) -> (Float64, Float64)\n\nGiven the (1-based) index of a pixel in a Healpix map in ring order, return a pair containing the (colatitude, longitude) angles corresponding to its center, both expressed in radians.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.ring2nest-Tuple{Resolution,Any}",
    "page": "Pixel functions",
    "title": "Healpix.ring2nest",
    "category": "method",
    "text": "ring2nest(resol::Resolution, ipix) :: Int\n\nConvert the number of a pixel from RING to NESTED scheme.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.nest2ring-Tuple{Resolution,Any}",
    "page": "Pixel functions",
    "title": "Healpix.nest2ring",
    "category": "method",
    "text": "nest2ring(resol::Resolution, ipix) :: Int\n\nConvert the number of a pixel from NESTED to RING scheme.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.pix2ringpos-Tuple{Resolution,Any}",
    "page": "Pixel functions",
    "title": "Healpix.pix2ringpos",
    "category": "method",
    "text": "pix2ringpos(resol::Resolution, pixel)\n\nGiven the (1-based) index of a pixel in a Healpix map in ring order, return a pair of numbers (n, i, j) whose meaning is the following:\n\nn can be one of the symbols :northcap, :equator, or :southcap, representing the region of the sky\ni is the ring index, from 1 to 4NSIDE - 1\nj is the pixel-in-ring index\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.pix2xyfNest-Tuple{Resolution,Any}",
    "page": "Pixel functions",
    "title": "Healpix.pix2xyfNest",
    "category": "method",
    "text": "pix2xyfNest(resol::Resolution, ipix) :: (Int, Int, Int)\n\nConvert a pixel number into (x, y, face), using NESTED ordering.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.pix2xyfRing-Tuple{Resolution,Any}",
    "page": "Pixel functions",
    "title": "Healpix.pix2xyfRing",
    "category": "method",
    "text": "pix2xyfRing(resol::Resolution, ipix) :: (Int, Int, Int)\n\nConvert a pixel number into (x, y, face), using RING ordering.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.xyf2pixNest-Tuple{Resolution,Any,Any,Any}",
    "page": "Pixel functions",
    "title": "Healpix.xyf2pixNest",
    "category": "method",
    "text": "xyf2pixNest(resol::Resolution, ix, iy, facenum) :: Int\n\nConvert (x, y, face) into a pixel number, using NESTED ordering.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Healpix.xyf2pixRing-Tuple{Resolution,Any,Any,Any}",
    "page": "Pixel functions",
    "title": "Healpix.xyf2pixRing",
    "category": "method",
    "text": "xyf2pixRing(resol::Resolution, ix, iy, facenum) :: Int\n\nConvert (x, y, face) into a pixel number, using RING ordering.\n\n\n\n\n\n"
},

{
    "location": "pixelfunc/#Reference-1",
    "page": "Pixel functions",
    "title": "Reference",
    "category": "section",
    "text": "ang2vec(theta, phi)\nvec2ang(x, y, z)\nang2pixNest(resol::Resolution, theta, phi)\nang2pixRing(resol::Resolution, theta, phi)\npix2angNest(resol::Resolution, pixel)\npix2angRing(resol::Resolution, pixel)\nring2nest(resol::Resolution, ipix)\nnest2ring(resol::Resolution, ipix)\npix2ringpos(resol::Resolution, pixel)\npix2xyfNest(resol::Resolution, ipix)\npix2xyfRing(resol::Resolution, ipix)\nxyf2pixNest(resol::Resolution, ix, iy, facenum)\nxyf2pixRing(resol::Resolution, ix, iy, facenum)"
},

{
    "location": "mapfunc/#",
    "page": "Map functions",
    "title": "Map functions",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Healpix\nend"
},

{
    "location": "mapfunc/#Map-functions-1",
    "page": "Map functions",
    "title": "Map functions",
    "category": "section",
    "text": "Functions like pix2angNest and ang2pixNest fully define the Healpix tessellation scheme. They are however extremely impractical in a number of situations. It happens often that a large fraction of pixels in a map need to be processed together. Healpix.jl introduces the Map{T, O <: Order} type, which acts as a collection of all the pixels on the sphere. A Map type holds the value of all the pixels in its pixels field, and it keeps track of the ordering (either RING or NESTED). Here is an example that shows how to create a map and initialize it:nside = 32\nm = Map{Float64, RingOrder}(nside)\nm.pixels[:] = 1.0  # Set all pixels to 1Healpix.jl defines the basic operations on maps (sum, subtraction, multiplication, division). These operations can either combine two maps or a map and a scalar value:mollweide(m * 2.0)\nmollweide(m * m)The Map{T, O <: Order} is derived from the abstract type GenericMap{T}, which does not encode the ordering. It is useful for functions that can either work on ring/nested-ordered maps but cannot be executed on plain generic arrays:# Return the number of pixels in the map, regardless of its ordering\nmaplength(m::GenericMap{T}) where T = length(m)\n\n# This returns 12\nmaplength(Map{Float64, RingOrder}(1))\n\n# This too returns 12\nmaplength(Map{Float64, NestedOrder}(1))\n\n# This fails\nmaplength(zeros(Float64, 12))Healpix.jl implements the PolarizedMap{T, O <: Order} type as well. This encodes three maps containing the I/Q/U signal: the intensity (I), and the Q and U Stokes parameters. The three maps must have the same resolution.GenericMap\nMap\nPolarizedMap"
},

{
    "location": "mapfunc/#Healpix.Order",
    "page": "Map functions",
    "title": "Healpix.Order",
    "category": "type",
    "text": "Abstract type representing the ordering of pixels in a Healpix map. See also RingOrder and NestedOrder.\n\n\n\n\n\n"
},

{
    "location": "mapfunc/#Healpix.RingOrder",
    "page": "Map functions",
    "title": "Healpix.RingOrder",
    "category": "type",
    "text": "The RingOrder type should be used when creating Map types in order to specify that the pixels in the map are sorted in ring ordering. (See also NestedOrder.)\n\n\n\n\n\n"
},

{
    "location": "mapfunc/#Healpix.NestedOrder",
    "page": "Map functions",
    "title": "Healpix.NestedOrder",
    "category": "type",
    "text": "The NestedOrder type should be used when creating Map types in order to specify that the pixels in the map are sorted in ring ordering. (See also RingOrder.)\n\n\n\n\n\n"
},

{
    "location": "mapfunc/#Encoding-the-order-1",
    "page": "Map functions",
    "title": "Encoding the order",
    "category": "section",
    "text": "Healpix.jl distinguishes between RING and NEST orderings using Julia\'s typesystem. The abstract type Order has two descendeants, RingOrder and NestedOrder, which are used to instantiate objects of type Map.Order\nRingOrder\nNestedOrder"
},

{
    "location": "mapfunc/#Healpix.pix2ang",
    "page": "Map functions",
    "title": "Healpix.pix2ang",
    "category": "function",
    "text": "pix2ang{T, O <: Order}(map::Map{T, O}, ipix) -> (Float64, Float64)\npix2ang{T, O <: Order}(map::PolarizedMap{T, O}, ipix) -> (Float64, Float64)\n\nReturn the pair (theta, phi), where theta is the colatitude and phi the longitude of the direction of the pixel center with index ipix.\n\n\n\n\n\n"
},

{
    "location": "mapfunc/#Healpix.ang2pix",
    "page": "Map functions",
    "title": "Healpix.ang2pix",
    "category": "function",
    "text": "ang2pix{T, O <: Order}(map::Map{T, O}, theta, phi)\nang2pix{T, O <: Order}(map::PolarizedMap{T, O}, theta, phi)\n\nConvert the direction specified by the colatitude theta (∈ [0, π]) and the longitude phi (∈ [0, 2π]) into the index of the pixel in the Healpix map map.\n\n\n\n\n\n"
},

{
    "location": "mapfunc/#Pixel-functions-1",
    "page": "Map functions",
    "title": "Pixel functions",
    "category": "section",
    "text": "When working with maps, it is not needed to pick between ang2pixNest and ang2pixRing because a Map type already encodes the ordering. Functions pix2ang and ang2pix always choose the correct ordering, but they require a Map instead of a Resolution as their first argument.pix2ang\nang2pix"
},

{
    "location": "mapfunc/#Healpix.saveToFITS",
    "page": "Map functions",
    "title": "Healpix.saveToFITS",
    "category": "function",
    "text": "saveToFITS(map::Map{T, O}, filename::AbstractString, typechar=\"D\", unit=\"\", extname=\"MAP\") where {T <: Number, O <: Order}\nsaveToFITS(map::PolarizedMap{T, O}, filename::AbstractString, typechar=\"D\", unit=\"\", extname=\"MAP\") where {T <: Number, O <: Order}\n\nSave a map into a FITS file. The name of the file is specified in filename; if it begins with !, existing files will be overwritten without warning. The parameter typechar specifies the data type to be used in the FITS file: the default (D) will save 64-bit floating-point values. See the CFITSIO documentation for other values. The keyword unit specifies the measure unit used for the pixels in the map. The keyword extname specifies the name of the HDU where the map pixels will be written.\n\n\n\n\n\n"
},

{
    "location": "mapfunc/#Healpix.savePixelsToFITS",
    "page": "Map functions",
    "title": "Healpix.savePixelsToFITS",
    "category": "function",
    "text": "savePixelsToFITS(map::Map{T}, f::FITSIO.FITSFile, column) where {T <: Number}\n\nSave the pixels of map into the column with index/name column in the FITS file, which must have been already opened.\n\n\n\n\n\n"
},

{
    "location": "mapfunc/#Healpix.readMapFromFITS",
    "page": "Map functions",
    "title": "Healpix.readMapFromFITS",
    "category": "function",
    "text": "readMapFromFITS{T <: Number}(f::FITSIO.FITSFILE, column, t::Type{T})\nreadMapFromFITS{T <: Number}(fileName::String, column, t::Type{T})\n\nRead a Healpix map from the specified (1-base indexed) column in a FITS file. The values will be read as numbers of type T. If the code fails, FITSIO will raise an exception. (Refer to the FITSIO library for more information.)\n\n\n\n\n\n"
},

{
    "location": "mapfunc/#Loading-and-saving-maps-1",
    "page": "Map functions",
    "title": "Loading and saving maps",
    "category": "section",
    "text": "Healpix.jl implements a number of functions to save maps in FITS files.saveToFITSFunction savePixelsToFITS is a low-level function. It knows nothing about the ordering schema used for the pixels, so the caller should manually write the ORDERING keyword in the HDU header by itself.savePixelsToFITSTo load a map from a FITS file, you can use readMapFromFITS.readMapFromFITS"
},

{
    "location": "mapfunc/#Healpix.conformables",
    "page": "Map functions",
    "title": "Healpix.conformables",
    "category": "function",
    "text": "conformables{T, S, O1 <: Order, O2 <: Order}(map1::Map{T, O1},\n                                             map2::Map{S, O2}) -> Bool\nconformables{T, S, O1 <: Order, O2 <: Order}(map1::PolarizedMap{T, O1},\n                                             map2::PolarizedMap{S, O2}) -> Bool\n\nDetermine if two Healpix maps are \"conformables\", i.e., if their shape and ordering are the same.\n\n\n\n\n\n"
},

{
    "location": "mapfunc/#Testing-for-conformability-1",
    "page": "Map functions",
    "title": "Testing for conformability",
    "category": "section",
    "text": "It often happens that two Healpix maps need to be combined together: for instance, pixels on a sky map might need to be masked using a sky mask, or one map might need to be subtracted from another one. «Conformability» means that the operation between the two maps can be done directly on the pixels, without oordering or resolution conversions. The function conformables checks this.using Healpix # hide\nm1 = Map{Float64, RingOrder}(1)\nm2 = Map{Float64, RingOrder}(1)\nm3 = Map{Float64, NestedOrder}(1)\nm4 = Map{Float64, NestedOrder}(2)\nconformables(m1, m2)\nconformables(m1, m3)\nconformables(m1, m4)conformables"
},

{
    "location": "mapfunc/#Healpix.tod2map",
    "page": "Map functions",
    "title": "Healpix.tod2map",
    "category": "function",
    "text": "tod2map{T,O}(pixidx, tod::Array{T}; nside=128) :: (map, hits)\n\nCreate a binned map for a TOD and return a tuple containing the map itself and the hit map.\n\n\n\n\n\n"
},

{
    "location": "mapfunc/#Healpix.combinemaps!",
    "page": "Map functions",
    "title": "Healpix.combinemaps!",
    "category": "function",
    "text": "combinemaps{T, O, H}(destmap::Map{T, O}, desthitmap::Map{H, O}, othermap::Map{T, O}, otherhitmap::Map{H, O})\n\nSum \"othermap\" to \"destmap\", assuming that both maps have been produced by binning TODs. The parameters desthitmap and otherhitmap are the two hit maps. At the end of the call, destmap and desthitmap are updated.\n\n\n\n\n\n"
},

{
    "location": "mapfunc/#Map-making-1",
    "page": "Map functions",
    "title": "Map-making",
    "category": "section",
    "text": "Map-making is the process of converting a time series of measurements into a sky map. The most basic form of map-making is the so-called \"binning\", where samples in the time stream falling within the same sky pixel are averaged. This map-making algorithm is strictly accurate only if the noise in the time stream is white.Healpix.jl implements two functions to perform binning, tod2map and combinemaps!.tod2map\ncombinemaps!"
},

{
    "location": "alm/#",
    "page": "Spherical harmonics",
    "title": "Spherical harmonics",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Healpix\nend"
},

{
    "location": "alm/#Healpix.Alm",
    "page": "Spherical harmonics",
    "title": "Healpix.Alm",
    "category": "type",
    "text": "An array of a_ℓm numbers.\n\n\n\n\n\n"
},

{
    "location": "alm/#Healpix.numberOfAlms",
    "page": "Spherical harmonics",
    "title": "Healpix.numberOfAlms",
    "category": "function",
    "text": "numberOfAlms(lmax::Integer, mmax::Integer) -> Integer\nnumberOfAlms(lmax::Integer) -> Integer\n\nReturn the size of the array of complex numbers needed to store the a_lm coefficients in the range of ℓ and m specified by lmax and mmax. If mmax is not specified, it is assumed to be equal to lmax. If lmax and mmax are inconsistent or negative, a DomainError exception is thrown.\n\n\n\n\n\n"
},

{
    "location": "alm/#Spherical-harmonics-1",
    "page": "Spherical harmonics",
    "title": "Spherical harmonics",
    "category": "section",
    "text": "The support for spherical harmonics in Healpix.jl is still woefully inadequate. Only a few functions to load and store harmonic coefficients are available. Everything revolves around the Alm type:AlmThe number of coefficients in a spherical harmonic expansion is infinite. For obvious reasons, Healpix.jl only allows to store band-limited expansions. The function numberOfAlms returns the number of floating-point numbers used to store the expansion, as a function of the maximum value for ell and m.numberOfAlms"
},

{
    "location": "alm/#Healpix.readAlmFromFITS",
    "page": "Spherical harmonics",
    "title": "Healpix.readAlmFromFITS",
    "category": "function",
    "text": "readAlmFromFITS{T <: Complex}(f::FITSIO.FITSFile, t::Type{T}) -> Alm{T}\nreadAlmFromFITS{T <: Complex}(fileName::String, t::Type{T}) -> Alm{T}\n\nRead a set of a_ℓm coefficients from a FITS file. If the code fails, FITSIO will raise an exception. (Refer to the FITSIO library for more information.)\n\n\n\n\n\n"
},

{
    "location": "alm/#Loading-and-saving-harmonic-coefficients-1",
    "page": "Spherical harmonics",
    "title": "Loading and saving harmonic coefficients",
    "category": "section",
    "text": "readAlmFromFITS"
},

{
    "location": "visualization/#",
    "page": "Visualization",
    "title": "Visualization",
    "category": "page",
    "text": ""
},

{
    "location": "visualization/#Visualization-functions-1",
    "page": "Visualization",
    "title": "Visualization functions",
    "category": "section",
    "text": "Healpix.jl uses RecipesBase to display maps. You need to import Plots in order to display maps, using the plot functions.  Maps are internally treated as heatmaps, so your backend should support this kind of visualization: at the moment, this is true for GR, PlotLy and PyPlot.using Healpix\nusing Plots\ngr()  # Use the GR backend\n\nnside = 8\nm = Map{Float64, RingOrder}(nside)\nm.pixels[:] = 1:length(m.pixels)\nplot(m)\nsavefig(joinpath(\"images\", \"mollweide.png\")) # hide(Image: )A call to plot can provide two additional arguments:A carthographic projection (see below).\nA dictionary containing parameters to be used by the carthographic projection.The following example displays the same map in orthographic coordinates:plot(m, orthographic)\nsavefig(joinpath(\"images\", \"orthographic.png\")) # hide(Image: )"
},

{
    "location": "visualization/#Healpix.project",
    "page": "Visualization",
    "title": "Healpix.project",
    "category": "function",
    "text": "project(invprojfn, m::Map{T, O}, bmpwidth, bmpheight; kwargs...) where {T <: Number, O <: Order}\n\nReturn a 2D bitmap (array) containing a cartographic projection of the map and a 2D bitmap containing a boolean mask. The size of the bitmap is bmpwidth×bmpheight pixels. The function projfn must be a function which accepts as input two parameters x and y (numbers between -1 and 1).\n\nThe following keywords can be used in the call:\n\ncenter: 2-tuple specifying the location (colatitude, longitude) of the sky point that is to be placed in the middle of the image (in radians)\nunseen: by default, Healpix maps use the value -1.6375e+30 to mark unseen pixels. You can specify a different value using this keyword. This should not be used in common applications.\n\nReturn a Array{Union{Missing, Float32}} containing the intensity of each pixel. Pixels falling outside the projection are marked as NaN, and unseen pixels are marked as missing.\n\n\n\n\n\n"
},

{
    "location": "visualization/#Cartographic-projections-1",
    "page": "Visualization",
    "title": "Cartographic projections",
    "category": "section",
    "text": "Plotting is based on project, which takes a map as input and produces a 2-D bitmap containing a representation of the map suitable to be shown using Plots.Although the easiest way to plot a map is to use plot, project might be suitable in those cases where you are just interested in a 2D bitmap. It requires a inverse projection function (mapping the 2D plane to a point on the sphere) and the size of the bitmap, and it returns three values:A 2-D bitmap containing the color level of each pixel. Unseen pixels (e.g., those falling outside the ellipse in a Mollweide projection) are marked as NaN, as well as unseen pixels;\nA 2-D bitmap of Bool values, telling which pixels in the map are masked, i.e., they are marked as UNSEEN, NaN or missing in the Healpix map;\nA Bool flag telling if there is any masked value in the mask (2nd return value, see above). This parameter is returned to optimize calls to plot, but it is obviously redundant.Consider this example:using Healpix\n\nm = Map{Float64, RingOrder}(1)\n# Plot the map on a 20×20 bitmap using an\n# equirectangular projection\nimage, mask, maskflag = project(equiprojinv, m, 20, 20)A number of parameters can be passed to project, in order to taylor the representation. They must not be passed as keyword arguments, because this would clash with the way plot recipes work; instead, you must use a dictionary:# Return a 2-D bitmap of 16-bit floating-point values\nimage, _, _ = project(equiprojinv, m, 20, 20,\n                      Dict(:desttype => Float16))The following dictionary keys are available::desttype: type used for the pixels in the 2-D bitmap returned by project. Default is Float32;\n:unseen: the value marking pixels as unseen, i.e., masked. The default is -1.6375e+30, to preserve compatibility with other Healpix libraries.\n:center: currently this is used only with orthographic projections. It specifies the coordinates of the center of the image (colatitude and longitude, both in radians).project"
},

{
    "location": "visualization/#Healpix.mollweide",
    "page": "Visualization",
    "title": "Healpix.mollweide",
    "category": "function",
    "text": "mollweide(m::Map{T,O}, projparams = Dict()) where {T <: Number,O <: Order}\n\nHigh-level wrapper around project for Mollweide projections.\n\nThe following parameters can be set in the projparams dictionary:\n\nwidth: width of the image, in pixels (default: 720 pixels)\nheight: height of the image, in pixels; if not specified, it will be assumed to be equal to width\n\n\n\n\n\n"
},

{
    "location": "visualization/#Healpix.equirectangular",
    "page": "Visualization",
    "title": "Healpix.equirectangular",
    "category": "function",
    "text": "equirectangular(m::Map{T,O}; kwargs...) where {T <: Number, O <: Order}\n\nHigh-level wrapper around project for equirectangular projections.\n\n\n\n\n\n"
},

{
    "location": "visualization/#Healpix.orthographic",
    "page": "Visualization",
    "title": "Healpix.orthographic",
    "category": "function",
    "text": "orthographic(m::Map{T,O}, projparams = Dict()) where {T <: Number,O <: Order}\n\nHigh-level wrapper around project for orthographic projections.\n\nThe following parameters can be set in the projparams dictionary:\n\nwidth: width of the image, in pixels (default: 720 pixels)\nheight: height of the image, in pixels; if not specified, it will be assumed to be equal to width\ncenter: position of the pixel in the middle of the left globe (latitude and longitude).\n\n\n\n\n\n"
},

{
    "location": "visualization/#Healpix.orthographic2",
    "page": "Visualization",
    "title": "Healpix.orthographic2",
    "category": "function",
    "text": "orthographic2(m::Map{T,O}, projparams = Dict()) where {T <: Number,O <: Order}\n\nHigh-level wrapper around project for stereo orthographic projections.\n\nThe following parameters can be set in the projparams dictionary:\n\nwidth: width of the image, in pixels (default: 720 pixels)\nheight: height of the image, in pixels; if not specified, it will be assumed to be equal to width\ncenter: position of the pixel in the middle of the left globe (latitude and longitude). Default is (0, 0).\n\n\n\n\n\n"
},

{
    "location": "visualization/#Healpix.gnomonic",
    "page": "Visualization",
    "title": "Healpix.gnomonic",
    "category": "function",
    "text": "gnomonic(m::Map{T,O}, projparams = Dict()) where {T <: Number,O <: Order}\n\nHigh-level wrapper around project for gnomonic projections.\n\nThe following parameters can be set in the projparams dictionary:\n\nwidth: width of the image, in pixels (default: 720 pixels)\nheight: height of the image, in pixels; if not specified, it will be assumed to be equal to width\ncenter: position and orientation of the pixel in the middle. It is a 3-element tuple containing:\nThe latitude of the pixel, in radians\nThe longitude of the pixel, in radians\nThe rotation to be applied to the image, in radians\nfov_rad: size of the image along the x and y axes, in radians (default: 15°)\n\nExample\n\nplot(m, gnomonic, Dict(:fov_rad = deg2rad(1.5), :center = (0, 0, deg2rad(45))))\n\n\n\n\n\n"
},

{
    "location": "visualization/#Projection-functions-1",
    "page": "Visualization",
    "title": "Projection functions",
    "category": "section",
    "text": "Functions mollweide, equirectangular, and orthographic can be passed as parameters to plot.mollweide\nequirectangular\northographic\northographic2\ngnomonicThey are based on inverse projection functions, i.e., functions that take a mollweideprojinv\nequiprojinv\northoinv\northo2inv\ngnominv"
},

{
    "location": "misc/#",
    "page": "Miscellanea",
    "title": "Miscellanea",
    "category": "page",
    "text": "DocTestSetup = quote\n    using Healpix\nend"
},

{
    "location": "misc/#Healpix.lat2colat",
    "page": "Miscellanea",
    "title": "Healpix.lat2colat",
    "category": "function",
    "text": "lat2colat(x)    colat2lat(x)\n\nConvert colatitude into latitude and vice versa. Both x and the result are expressed in radians.\n\n\n\n\n\n"
},

{
    "location": "misc/#Healpix.colat2lat",
    "page": "Miscellanea",
    "title": "Healpix.colat2lat",
    "category": "function",
    "text": "lat2colat(x)    colat2lat(x)\n\nConvert colatitude into latitude and vice versa. Both x and the result are expressed in radians.\n\n\n\n\n\n"
},

{
    "location": "misc/#General-purpose-functions-1",
    "page": "Miscellanea",
    "title": "General-purpose functions",
    "category": "section",
    "text": "Healpix.jl implements a few generic functions that can be helpful when doing calculations on the sphere.lat2colat\ncolat2lat"
},

]}
