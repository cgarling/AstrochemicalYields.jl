# Code for Nomoto 1997 SN-Ia yields
"""Module implementing the Nomoto+1997 SN-Ia yield tables."""
module Nomoto1997

using ..AstrochemicalYields: AbstractYield, α_elements
import ..AstrochemicalYields: isotopes, remnant_mass, preSN_mass, ejecta_mass, ejecta_metal_mass, ejecta_alpha_mass
# using Interpolations: interpolate, Linear, Gridded, extrapolate, Flat, Throw
using Printf: @sprintf
using StaticArrays: SVector

export Nomoto1997SNIa

const _Nomoto1997_models = (:TypeII, :W70, :W7, :WDD1, :WDD2, :WDD3)

"""
    Nomoto1997SNIa()
Loads the Nomoto+1997 SN-Ia yield table. There is no initial metallicity or mass dependence, so no interpolation is needed. Returned type has fields corresponding to the different models available in the table: `(TypeII, W70, W7, WDD1, WDD2, WDD3)`. Each of these entries is itself a `NamedTuple` with keys corresponding to the relevant isotopes.
"""
struct Nomoto1997SNIa{A} <: AbstractYield
    TypeII::A
    W70::A
    W7::A
    WDD1::A
    WDD2::A
    WDD3::A
end
function Nomoto1997SNIa()
    fname = joinpath(@__DIR__, "data", "Nomoto1997_SNIa.csv")
    lines = readlines(fname)
    isotopes = Vector{String}(undef, length(lines)-1)
    table = Matrix{Float64}(undef, length(isotopes), 6)
    for (i, line) in enumerate(lines[2:end])
        parts = split(line)
        m = match(r"^(\d+)([A-Za-z]+)$", parts[1])
        if isnothing(m)
            isotopes[i] = parts[1]
        else
            isotopes[i] = m.captures[2] * m.captures[1]
        end
        table[i, :] = parse.(Float64, replace.(parts[2:end], '−' => '-'))
    end
    nt = NamedTuple{Tuple(Symbol.(isotopes)), NTuple{length(isotopes), Float64}}
    # return (TypeII = nt(table[:,1]), W70 = nt(table[:,2]), W7 = nt(table[:,3]), WDD1 = nt(table[:,4]), WDD2 = nt(table[:,5]), WDD3 = nt(table[:,6]))
    return Nomoto1997SNIa((nt(table[:,i]) for i in 1:6)...)
end
isotopes(x::Nomoto1997SNIa) = keys(x.TypeII)
preSN_mass(x::Nomoto1997SNIa, Z=0, M=0) = 1.4 # Chandrasekhar mass
remnant_mass(x::Nomoto1997SNIa, Z=0, M=0) = 0.0 # Assuming no compact remnant
ejecta_mass(x::Nomoto1997SNIa, Z=0, M=0) = preSN_mass(x, Z, M)

function filter_metals(nt::NamedTuple)
    names = (:C12, :C13, :N14, :N15, :O16, :O17, :O18, :F19, :Ne20, :Ne21, :Ne22, :Na23, :Mg24, :Mg25, :Mg26, :Al27, :Si28, :Si29, :Si30, :P31, :S32, :S33, :S34, :S36, :Cl35, :Cl37, :Ar36, :Ar38, :Ar40, :K39, :K41, :Ca40, :Ca42, :Ca43, :Ca44, :Ca46, :Ca48, :Sc45, :Ti46, :Ti47, :Ti48, :Ti49, :Ti50, :V50, :V51, :Cr50, :Cr52, :Cr53, :Cr54, :Mn55, :Fe54, :Fe56, :Fe57, :Fe58, :Co59, :Ni58, :Ni60, :Ni61, :Ni62, :Ni64, :Cu63, :Cu65, :Zn64, :Zn66, :Zn67, :Zn68)
    return nt[names]
end
# filter_metals(x::Nomoto1997SNIa) = NamedTuple{_Nomoto1997_models}(filter_metals(getproperty(x, i)) for i in _Nomoto1997_models)
function filter_metals(x::Nomoto1997SNIa)
    names = (:C12, :C13, :N14, :N15, :O16, :O17, :O18, :F19, :Ne20, :Ne21, :Ne22, :Na23, :Mg24, :Mg25, :Mg26, :Al27, :Si28, :Si29, :Si30, :P31, :S32, :S33, :S34, :S36, :Cl35, :Cl37, :Ar36, :Ar38, :Ar40, :K39, :K41, :Ca40, :Ca42, :Ca43, :Ca44, :Ca46, :Ca48, :Sc45, :Ti46, :Ti47, :Ti48, :Ti49, :Ti50, :V50, :V51, :Cr50, :Cr52, :Cr53, :Cr54, :Mn55, :Fe54, :Fe56, :Fe57, :Fe58, :Co59, :Ni58, :Ni60, :Ni61, :Ni62, :Ni64, :Cu63, :Cu65, :Zn64, :Zn66, :Zn67, :Zn68)
    return NamedTuple{_Nomoto1997_models, NTuple{length(_Nomoto1997_models), NamedTuple{names, NTuple{length(names), Float64}}}}(getproperty(x, i)[names] for i in _Nomoto1997_models)
end
function ejecta_metal_mass(x::Nomoto1997SNIa, Z=0, M=0)
    f = filter_metals(x)
    # return (; (i => sum(f[i]) for i in _Nomoto1997_models)...)
    return NamedTuple{_Nomoto1997_models, NTuple{length(_Nomoto1997_models), Float64}}(sum(i) for i in f)
end

function filter_alpha(nt::NamedTuple)
    names = (:O16, :O17, :O18, :Ne20, :Ne21, :Ne22, :Mg24, :Mg25, :Mg26, :Si28, :Si29, :Si30, :S32, :S33, :S34, :S36, :Ar36, :Ar38, :Ar40, :Ca40, :Ca42, :Ca43, :Ca44, :Ca46, :Ca48, :Ti46, :Ti47, :Ti48, :Ti49, :Ti50)
    return nt[names]
end
function filter_alpha(x::Nomoto1997SNIa)
    names = (:O16, :O17, :O18, :Ne20, :Ne21, :Ne22, :Mg24, :Mg25, :Mg26, :Si28, :Si29, :Si30, :S32, :S33, :S34, :S36, :Ar36, :Ar38, :Ar40, :Ca40, :Ca42, :Ca43, :Ca44, :Ca46, :Ca48, :Ti46, :Ti47, :Ti48, :Ti49, :Ti50)
    return NamedTuple{_Nomoto1997_models, NTuple{length(_Nomoto1997_models), NamedTuple{names, NTuple{length(names), Float64}}}}(getproperty(x, i)[names] for i in _Nomoto1997_models)
end
ejecta_alpha_mass(x::Nomoto1997SNIa, Z=0, M=0) = NamedTuple{_Nomoto1997_models, NTuple{length(_Nomoto1997_models), Float64}}(sum(i) for i in filter_alpha(x))

end # module