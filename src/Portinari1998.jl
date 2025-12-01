"""Module implementing the Portinari+1998 core-collapse supernova yield tables."""
module Portinari1998

using ..AstrochemicalYields: AbstractYield, α_elements, extend_bounds
import ..AstrochemicalYields: isotopes, preSN_mass, remnant_mass, ejecta_mass, ejecta_metal_mass, ejecta_alpha_mass
using Interpolations: interpolate, Linear, Gridded, extrapolate, Flat, Throw
using Printf: @sprintf
using StaticArrays: SVector

export Portinari1998SN

const _Portinari1998_Zs = (0.0004, 0.004, 0.008, 0.02, 0.05)
const _Portinari1998_SN_M = (6.0, 7.0, 9.0, 12.0, 15.0, 20.0, 30.0, 40.0, 60.0, 100.0, 120.0, 150.0, 200.0, 300.0, 500.0, 1000.0)

struct Portinari1998Entry{A, C, D}
    Z::A
    table::C
    isotopes::D
end

function Portinari1998Entry(Z::Number)
    if !any(map(Base.Fix1(isapprox, Z), _Portinari1998_Zs))
        throw(ArgumentError("Argument `Z` must be approximately equal to one of $(_Portinari1998_Zs)."))
    end
    fname = joinpath(@__DIR__, "data", "Portinari1998", "SN_Z_"*@sprintf("%04d", round(Int, Z*10000))*".csv")
    lines = readlines(fname)
    if length(lines) != length(_Portinari1998_SN_M) + 1
        error("Wrong number of lines in Portinary 1998 data file $fname.")
    end
    isotopes = split(lines[1])[begin+1:end]
    for i in eachindex(isotopes)
        m = match(r"^(\d+)([A-Za-z]+)$", isotopes[i])
        if !isnothing(m)
            isotopes[i] = m.captures[2] * m.captures[1]
        end
    end
    isotopes[end] = "Mcut"
    table = Matrix{Float64}(undef, length(isotopes), length(_Portinari1998_SN_M))
    for (i, line) in enumerate(lines[2:end])
        parts = split(line)
        table[:, i] = parse.(Float64, replace.(parts[2:end], '−' => '-'))
    end
    return Portinari1998Entry(float(Z), table, isotopes)
end

"""
    Portinari1998SN(; bounds=Interpolations.Throw()) <: AbstractYield
Load the Portinari+1998 core-collapse supernova yield table. The yield table can be interpolated by calling it with the metal mass fraction `Z` and stellar mass `M` (in solar masses) of the progenitor. The keyword argument `bounds` should be a valid `Interpolations.jl` extrapolation specifier that will determine how the interpolation is extrapolated (e.g., `Flat()`).

```jldoctest
julia> n = Portinari1998SN();

julia> n(0.002, 13.5) isa NamedTuple
true
```
"""
struct Portinari1998SN{I, B} <: AbstractYield
    itp::B
end
function Portinari1998SN(; bounds=Throw())
    entries = Portinari1998Entry.(_Portinari1998_Zs)
    iso_mat = [SVector{length(entries[1].isotopes)}(i.table[:, j]) for i=entries, j=eachindex(_Portinari1998_SN_M)]
    iso_itp = interpolate((SVector(_Portinari1998_Zs), SVector(_Portinari1998_SN_M)), iso_mat, Gridded(Linear()))
    bounds = extend_bounds(bounds, length(entries[1].isotopes))
    iso_itp = extrapolate(iso_itp, bounds)
    isotopes = Tuple(Symbol.(entries[1].isotopes))
    return Portinari1998SN{isotopes, typeof(iso_itp)}(iso_itp)
end
# Base.show(io::IO, ::Portinari1998SN) = print(io, "Grid of core-collapse SN yields from Portinari+1998.")
isotopes(::Portinari1998SN{I}) where I = I
_nt(::Portinari1998SN{I}) where I = NamedTuple{I, NTuple{length(I), Float64}}
function (x::Portinari1998SN)(Z, M)
    return _nt(x)(Tuple(x.itp(Z, M)))
end
# Note that the ejecta tables of Portinari1998 are *total* ejecta in the sense that they include the pre-SN wind as well.
# They have a separate wind table (2) but it does not cover the low mass (6, 7) or high mass components so neglecting for now.
preSN_mass(x::Portinari1998SN, Z, M) = sum(x(Z, M)) 
remnant_mass(x::Portinari1998SN, Z, M) = x(Z, M).Mcut
function ejecta_mass(x::Portinari1998SN, Z, M)
    nt = x(Z, M)
    return sum(filter(!isnan, nt)) - nt.Mcut
end

# Filter to include only metals
function filter_metals(nt::NamedTuple)
    names = (:C12, :C13, :N14, :N15, :O16, :O17, :O18, :Ne20, :Ne22, :Mg24, :Si28, :S32, :Ca40, :Fe56)
    return nt[names]
end
ejecta_metal_mass(x::Portinari1998SN, Z, M) = sum(filter(!isnan, filter_metals(x(Z, M))))

# Filter to include only alpha elements
function filter_alpha(nt::NamedTuple)
    names = (:O16, :O17, :Ne20, :Ne22, :Mg24, :Si28, :Ca40)
    return nt[names]
end
ejecta_alpha_mass(x::Portinari1998SN, Z, M) = sum(filter(!isnan, filter_alpha(x(Z, M))))

end # module