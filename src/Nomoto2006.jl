# Code for Nomoto2006 SN-II yields
# The mass cut is the mass coordinate (usually in solar masses) that separates the ejecta from the remnant

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
        m = match(r"^(\d+)([A-Za-z]+)$", parts[1])
        if isnothing(m)
            isotopes[i] = parts[1]
        else
            isotopes[i] = m.captures[2] * m.captures[1]
        end
        table[i, :] = parse.(Float64, replace.(parts[2:end], '−' => '-'))
    end
    table[end, :] .= Mcut
    isotopes[end] = "Mcut"
    return Nomoto2006Entry(Z, M, E, table, isotopes)
end

"""
    Nomoto2006SN()
Load the Nomoto+2006 core-collapse supernova yield table. The yield table can be interpolated by calling it with the metal mass fraction `Z` and stellar mass `M` (in solar masses) of the progenitor.

```jldoctest
julia> n = Nomoto2006SN();

julia> n(0.002, 13.5) isa NamedTuple
true
"""
struct Nomoto2006SN{I, B} <: AbstractYield
    itp::B
    # isotopes::NTuple{N, Symbol} # This is now part of the type signature (I) for performance
end
function Nomoto2006SN()
    entries = Nomoto2006Entry.(_Nomoto2006_Zs)
    iso_mat = [SVector{length(entries[1].isotopes)}(i.table[:, j]) for i=entries, j=eachindex(entries[1].M)]
    iso_itp = interpolate((SVector(_Nomoto2006_Zs), SVector(_Nomoto2006_SN_M)), iso_mat, Gridded(Linear()))
    iso_itp = extrapolate(iso_itp, Throw())
    isotopes = Tuple(Symbol.(entries[1].isotopes))
    return Nomoto2006SN{isotopes, typeof(iso_itp)}(iso_itp)
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
# function filter_metals(nt)
#     exclude_keys = (:p, :d, Symbol("3He"), Symbol("4He"), :Mcut) # Symbol("6Li"), Symbol("7Li")
#     exclude_keys = (:p, :d, Symbol("He3"), Symbol("He4"), :Mcut) # Symbol("6Li"), Symbol("7Li")
#     # return (; (k => nt[k] for k in keys(nt) if k ∉ exclude_keys)...)
#     # Get indices into keys(nt) that we want to use
#     # return [i for i in eachindex(values(nt)) if keys(nt)[i] ∉ exclude_keys]
#     # Get entries in keys(nt) that we want to use
#     # return [k for k in keys(nt) if k ∉ exclude_keys]
# end

function filter_metals(nt::NamedTuple)
    # idxs = 5:80
    # return SVector(values(nt))[idxs]
    # return NamedTuple{A[idxs], NTuple{length(idxs), Float64}}(values(nt)[idxs])
    # names = (Symbol("6Li"), Symbol("7Li"), Symbol("9Be"), Symbol("10B"), Symbol("11B"), Symbol("12C"), Symbol("13C"), Symbol("14N"), Symbol("15N"), Symbol("16O"), Symbol("17O"), Symbol("18O"), Symbol("19F"), Symbol("20Ne"), Symbol("21Ne"), Symbol("22Ne"), Symbol("23Na"), Symbol("24Mg"), Symbol("25Mg"), Symbol("26Mg"), Symbol("26Al"), Symbol("27Al"), Symbol("28Si"), Symbol("29Si"), Symbol("30Si"), Symbol("31P"), Symbol("32S"), Symbol("33S"), Symbol("34S"), Symbol("36S"), Symbol("35Cl"), Symbol("37Cl"), Symbol("36Ar"), Symbol("38Ar"), Symbol("40Ar"), Symbol("39K"), Symbol("40K"), Symbol("41K"), Symbol("40Ca"), Symbol("42Ca"), Symbol("43Ca"), Symbol("44Ca"), Symbol("46Ca"), Symbol("48Ca"), Symbol("45Sc"), Symbol("46Ti"), Symbol("47Ti"), Symbol("48Ti"), Symbol("49Ti"), Symbol("50Ti"), Symbol("50V"), Symbol("51V"), Symbol("50Cr"), Symbol("52Cr"), Symbol("53Cr"), Symbol("54Cr"), Symbol("55Mn"), Symbol("54Fe"), Symbol("56Fe"), Symbol("57Fe"), Symbol("58Fe"), Symbol("59Co"), Symbol("58Ni"), Symbol("60Ni"), Symbol("61Ni"), Symbol("62Ni"), Symbol("64Ni"), Symbol("63Cu"), Symbol("65Cu"), Symbol("64Zn"), Symbol("66Zn"), Symbol("67Zn"), Symbol("68Zn"), Symbol("70Zn"), Symbol("69Ga"), Symbol("71Ga"))
    names = (:Li6, :Li7, :Be9, :B10, :B11, :C12, :C13, :N14, :N15, :O16, :O17, :O18, :F19, :Ne20, :Ne21, :Ne22, :Na23, :Mg24, :Mg25, :Mg26, :Al26, :Al27, :Si28, :Si29, :Si30, :P31, :S32, :S33, :S34, :S36, :Cl35, :Cl37, :Ar36, :Ar38, :Ar40, :K39, :K40, :K41, :Ca40, :Ca42, :Ca43, :Ca44, :Ca46, :Ca48, :Sc45, :Ti46, :Ti47, :Ti48, :Ti49, :Ti50, :V50, :V51, :Cr50, :Cr52, :Cr53, :Cr54, :Mn55, :Fe54, :Fe56, :Fe57, :Fe58, :Co59, :Ni58, :Ni60, :Ni61, :Ni62, :Ni64, :Cu63, :Cu65, :Zn64, :Zn66, :Zn67, :Zn68, :Zn70, :Ga69, :Ga71)
    return nt[names]
    # return NamedTuple{names, NTuple{length(names), Float64}}(values(nt)[idxs])
end

# @generated function filter_metals(nt::NamedTuple{A}) where {A}
#     idxs = 5:80
#     keys = ntuple(i -> A[i], length(idxs))
#     vals = [:(nt.$k) for k in keys]
#     return :(NamedTuple{Tuple($(map(QuoteNode, keys)...))}( (; $(vals...)) ))
# end


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

end # module