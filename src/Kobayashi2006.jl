# Code for Kobayashi2006 SN-II yields
# The mass cut is the mass coordinate (usually in solar masses) that separates the ejecta from the remnant

module Kobayashi2006

using ..AstrochemicalYields: AbstractYield, α_elements
import ..AstrochemicalYields: isotopes, preSN_mass, remnant_mass, ejecta_mass, ejecta_metal_mass, ejecta_alpha_mass
using Interpolations: interpolate, Linear, Gridded, extrapolate, Flat, Throw
# using DataInterpolationsND: NDInterpolation, LinearInterpolationDimension, BSplineInterpolationDimension
using Printf: @sprintf
using StaticArrays: SVector

export Kobayashi2006SN

const _Kobayashi2006_Zs = (0.0, 0.001, 0.004, 0.02)
const _Kobayashi2006_SN_M = (13.0, 15.0, 18.0, 20.0, 25.0, 30.0, 40.0)

struct Kobayashi2006Entry{A, C, D}
    Z::A
    table::C
    isotopes::D
end

function Kobayashi2006Entry(Z::Number)
    if !any(map(Base.Fix1(isapprox, Z), _Kobayashi2006_Zs))
        throw(ArgumentError("Argument `Z` must be approximately equal to one of $(_Kobayashi2006_Zs)."))
    end
    fname = joinpath(@__DIR__, "data", "Kobayashi2006", "SN_Z_"*@sprintf("%03d", round(Int, Z*1000))*".csv")
    lines = readlines(fname)
    isotopes = Vector{String}(undef, length(lines))
    table = Matrix{Float64}(undef, length(isotopes), length(_Kobayashi2006_SN_M))
    for (i, line) in enumerate(lines)
        parts = split(line)
        # Replace carats and reverse order of element and isotope number
        parts[2] = replace(parts[2], "^" => "")
        m = match(r"^(\d+)([A-Za-z]+)$", parts[2])
        if isnothing(m)
            isotopes[i] = replace(parts[2], "_" => "")
        else
            isotopes[i] = m.captures[2] * m.captures[1]
        end
        table[i, :] = parse.(Float64, parts[3:end])
    end
    return Kobayashi2006Entry(float(Z), table, isotopes)
end

"""
    Kobayashi2006SN()
Load the Kobayashi+2006 core-collapse supernova yield table (this is mostly the same as Nomoto+2006). The yield table can be interpolated by calling it with the metal mass fraction `Z` and stellar mass `M` (in solar masses) of the progenitor.

```jldoctest
julia> n = Kobayashi2006SN();

julia> n(0.002, 13.5) isa NamedTuple
true
"""
struct Kobayashi2006SN{I, B} <: AbstractYield
    itp::B
end
function Kobayashi2006SN()
    entries = Kobayashi2006Entry.(_Kobayashi2006_Zs)
    iso_mat = [SVector{length(entries[1].isotopes)}(i.table[:, j]) for i=entries, j=eachindex(_Kobayashi2006_SN_M)]
    iso_itp = interpolate((SVector(_Kobayashi2006_Zs), SVector(_Kobayashi2006_SN_M)), iso_mat, Gridded(Linear()))
    iso_itp = extrapolate(iso_itp, Throw())
    isotopes = Tuple(Symbol.(entries[1].isotopes))
    return Kobayashi2006SN{isotopes, typeof(iso_itp)}(iso_itp)
end
# Base.show(io::IO, ::Kobayashi2006SN) = print(io, "Grid of core-collapse SN yields from Kobayashi+2006.")
isotopes(::Kobayashi2006SN{I}) where I = I
_nt(::Kobayashi2006SN{I}) where I = NamedTuple{I, NTuple{length(I), Float64}}
function (x::Kobayashi2006SN)(Z, M)
    return _nt(x)(Tuple(x.itp(Z, M)))
end
preSN_mass(x::Kobayashi2006SN, Z, M) = x(Z, M).Mfinal
remnant_mass(x::Kobayashi2006SN, Z, M) = x(Z, M).Mcut
function ejecta_mass(x::Kobayashi2006SN, Z, M)
    nt = x(Z, M)
    return nt.Mfinal - nt.Mcut
end
# Filter to include only metals
function filter_metals(nt::NamedTuple)
    names = (:Li6, :Li7, :Be9, :B10, :B11, :C12, :C13, :N14, :N15, :O16, :O17, :O18, :F19, :Ne20, :Ne21, :Ne22, :Na23, :Mg24, :Mg25, :Mg26, :Al27, :Si28, :Si29, :Si30, :P31, :S32, :S33, :S34, :S36, :Cl35, :Cl37, :Ar36, :Ar38, :Ar40, :K39, :K40, :K41, :Ca40, :Ca42, :Ca43, :Ca44, :Ca46, :Ca48, :Sc45, :Ti46, :Ti47, :Ti48, :Ti49, :Ti50, :V50, :V51, :Cr50, :Cr52, :Cr53, :Cr54, :Mn55, :Fe54, :Fe56, :Fe57, :Fe58, :Co59, :Ni58, :Ni60, :Ni61, :Ni62, :Ni64, :Cu63, :Cu65, :Zn64, :Zn66, :Zn67, :Zn68, :Zn70, :Ga69, :Ga71, :Ge70, :Ge72, :Ge73, :Ge74)
    return nt[names]
end
ejecta_metal_mass(x::Kobayashi2006SN, Z, M) = sum(values(filter_metals(x(Z, M))))

# Filter to include only alpha elements
function filter_alpha(nt::NamedTuple)
    names = (:O16, :O17, :O18, :Ne20, :Ne21, :Ne22, :Mg24, :Mg25, :Mg26, :Si28, :Si29, :Si30, :S32, :S33, :S34, :S36, :Ar36, :Ar38, :Ar40, :Ca40, :Ca42, :Ca43, :Ca44, :Ca46, :Ca48, :Ti46, :Ti47, :Ti48, :Ti49, :Ti50)
    return nt[names]
end
ejecta_alpha_mass(x::Kobayashi2006SN, Z, M) = sum(values(filter_alpha(x(Z, M))))

end # module