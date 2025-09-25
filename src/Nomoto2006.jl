# Code for Nomoto2006 SN-II yields
# The mass cut is the mass coordinate (usually in solar masses) that separates the ejecta from the remnant
# This code allows you to select out a subset of isotopes when constructing Nomoto2006SN;
# planning to do away with this in the future

module Nomoto2006

using ..AstrochemicalYields: AbstractYield, α_elements
import ..AstrochemicalYields: isotopes, remnant_mass, ejecta_mass, ejecta_metal_mass, ejecta_alpha_mass
using Interpolations: interpolate, Linear, Gridded, extrapolate, Flat, Throw
# using DataInterpolationsND: NDInterpolation, LinearInterpolationDimension, BSplineInterpolationDimension
using Printf: @sprintf
using StaticArrays: SVector

export Nomoto2006SN

const _Nomoto2006_Zs = (0.0, 0.001, 0.004, 0.02)
const _Nomoto2006_SN_M = (13.0, 15.0, 18.0, 20.0, 25.0, 30.0, 40.0)

struct Nomoto2006Entry{A, B <: AbstractArray{A}, C, D}
    Z::A
    M::B
    E::B
    table::C
    isotopes::D
end

function Nomoto2006Entry(Z::Number)
    if !any(map(Base.Fix1(isapprox, Z), _Nomoto2006_Zs))
        throw(ArgumentError("Argument `Z` must be approximately equal to one of $(_Nomoto2006_Zs)."))
    end
    fname = joinpath(@__DIR__, "data", "Nomoto2006", "SN_Z_"*@sprintf("%03d", round(Int, Z*1000))*".csv")
    lines = readlines(fname)
    Z = parse(Float64, split(lines[1])[end])
    M = parse.(Float64, split(lines[2])[2:end])
    E = parse.(Float64, split(lines[3])[2:end])
    Mcut = parse.(Float64, split(lines[4])[2:end])
    isotopes = Vector{String}(undef, length(lines)-4+1) # +1 for Mcut
    table = Matrix{Float64}(undef, length(isotopes), length(M))
    for (i, line) in enumerate(lines[5:end])
        parts = split(line)
        isotopes[i] = parts[1]
        table[i, :] = parse.(Float64, replace.(parts[2:end], '−' => '-'))
    end
    table[end, :] .= Mcut
    isotopes[end] = "Mcut"
    return Nomoto2006Entry(Z, M, E, table, isotopes)
end

struct Nomoto2006SN{I, B} <: AbstractYield
    itp::B
    # isotopes::NTuple{N, Symbol} # This is now part of the type signature (I) for performance
end
# Base.show(io::IO, ::Nomoto2006SN) = print(io, "Grid of core-collapse SN yields from Nomoto+2006.")
isotopes(::Nomoto2006SN{I}) where I = I
_nt(::Nomoto2006SN{I}) where I = NamedTuple{I, NTuple{length(I), Float64}}
function (x::Nomoto2006SN)(Z, M)
    return _nt(x)(Tuple(x.itp(Z, M)))
end
remnant_mass(x::Nomoto2006SN, Z, M) = x(Z, M).Mcut
ejecta_mass(x::Nomoto2006SN, Z, M) = M - remnant_mass(x, Z, M)
# Filter to include only metals
function filter_metals(nt)
    exclude_keys = (:p, :d, Symbol("3He"), Symbol("4He"), :Mcut) # Symbol("6Li"), Symbol("7Li")
    return (; (k => nt[k] for k in keys(nt) if k ∉ exclude_keys)...)
end
function ejecta_metal_mass(x::Nomoto2006SN{I}, Z, M) where I
    if length(I) == 81 # Fast path if no elements were excluded
        return sum(values(x(Z, M))[5:end-1])
    else
        return sum(values(filter_metals(x(Z, M))))
    end
end
# Filter to include only alpha elements
function filter_alpha(nt)
    include_keys = (Symbol("16O"), Symbol("17O"), Symbol("18O"), Symbol("20Ne"), Symbol("21Ne"), Symbol("22Ne"), Symbol("24Mg"), Symbol("25Mg"), Symbol("26Mg"), Symbol("28Si"), Symbol("29Si"), Symbol("30Si"), Symbol("32S"), Symbol("33S"), Symbol("34S"), Symbol("36S"), Symbol("36Ar"), Symbol("38Ar"), Symbol("40Ar"), Symbol("40Ca"), Symbol("42Ca"), Symbol("43Ca"), Symbol("44Ca"), Symbol("46Ca"), Symbol("48Ca"), Symbol("46Ti"), Symbol("47Ti"), Symbol("48Ti"), Symbol("49Ti"), Symbol("50Ti"))
    # This line will get the indices for the fast path
    # return [i for i in eachindex(values(nt)) if keys(nt)[i] ∈ include_keys]
    if length(nt) == 81 # Fast path if no elements were excluded
        return NamedTuple{include_keys, NTuple{length(include_keys), Float64}}(values(nt)[[14, 15, 16, 18, 19, 20, 22, 23, 24, 27, 28, 29, 31, 32, 33, 34, 37, 38, 39, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54]])
    # ks = keys(nt)
    # good = [any(map(x -> occursin(x, k), α_elements)) for k in string.(ks)] # Also grabs Sc, which we don't want
    # return (; (k => nt[k] for k in ks if k ∈ ks[good])...)
    else
        return NamedTuple(k => nt[k] for k in keys(nt) if k ∈ include_keys)
    end
end
ejecta_alpha_mass(x::Nomoto2006SN, Z, M) = sum(values(filter_alpha(x(Z, M))))

"""
    Nomoto2006SN()
Load the Nomoto+2006 core-collapse supernova yield table. The yield table can be interpolated by calling it with the metal mass fraction `Z` and stellar mass `M` (in solar masses) of the progenitor.

```jldoctest
julia> n = Nomoto2006SN();

julia> n(0.002, 13.5) isa NamedTuple
true
"""
function Nomoto2006SN()
    entries = Nomoto2006Entry.(_Nomoto2006_Zs)
    return Nomoto2006SN(entries, eachindex(entries[1].isotopes))
end
function Nomoto2006SN(isotopes)
    isotopes = string.(isotopes)
    entries = Nomoto2006Entry.(_Nomoto2006_Zs)
    good = [findfirst(x -> i == x, entries[1].isotopes) for i in isotopes]
    return Nomoto2006SN(entries, good)
end

function Nomoto2006SN(entries, good)
    N = length(good)
    iso_mat = [SVector{N}(i.table[good, j]) for i=entries, j=eachindex(entries[1].M)]
    iso_itp = interpolate((SVector(_Nomoto2006_Zs), SVector(_Nomoto2006_SN_M)), iso_mat, Gridded(Linear()))
    iso_itp = extrapolate(iso_itp, Throw())
    isotopes = Tuple(Symbol.(entries[1].isotopes[good]))
    return Nomoto2006SN{isotopes, typeof(iso_itp)}(iso_itp)
end

end # module