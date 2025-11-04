module AstrochemicalYields

using ArgCheck: @argcheck, @check
using Interpolations: interpolate, scale, Linear, BSpline, Cubic, Gridded, extrapolate, Flat, Throw
import InverseFunctions: inverse
using QuadGK: quadgk
using StaticArrays: SVector, @SVector

"""
    extend_bounds(x, N::Integer)
If `x isa Number`, returns an `SVector{N, typeof(x)}(x)`, else returns `x`. Used for adapting constant interpolations for Interpolations.jl into multi-valued `SVectors`.
"""
extend_bounds(x, N::Integer) = x
extend_bounds(x::Number, N::Integer) = @SVector fill(x, N)

"""
`AbstractYield` is the abstract supertype for all yield tables. Yield table subtypes should be made callable with initial metal mass fraction `Z` and mass `M` (in solar masses), returning the yield for all isotopes in units of solar masses. Subtypes should additionally implement the following methods:
 - [`isotopes`](@ref)
 - [`remnant_mass`](@ref)
 - [`ejecta_mass`](@ref)
 - [`ejecta_metal_mass`](@ref)
 - [`ejecta_alpha_mass`](@ref)
"""
abstract type AbstractYield end

"""List of α elements used in [`ejecta_alpha_mass`](@ref) taken from Nomoto+2006. C is sometimes included, but not here."""
const α_elements = ("O", "Ne", "Mg", "Si", "S", "Ar", "Ca", "Ti")

"""
    α_isotopes(x::AbstractYield)
Convenience function that returns the isotope symbols for the given `x::AbstractYield` for the α-elements $(join(α_elements, ", ")).
```jldoctest
julia> AstrochemicalYields.α_isotopes(Nomoto2006SN())
(:O16, :O17, :O18, :Ne20, :Ne21, :Ne22, :Mg24, :Mg25, :Mg26, :Si28, :Si29, :Si30, :S32, :S33, :S34, :S36, :Ar36, :Ar38, :Ar40, :Ca40, :Ca42, :Ca43, :Ca44, :Ca46, :Ca48, :Ti46, :Ti47, :Ti48, :Ti49, :Ti50)
```
"""
function α_isotopes(x::AbstractYield)
    isos = isotopes(x)
    names = Tuple(Iterators.flatten((filter(x -> begin y = String(x); startswith(y, i) && !startswith(y, "Sc") end, isos) for i in α_elements)))
    names = Tuple(unique(names))
    return names
end

"""
    isotopes(table::AbstractYield)
Returns a `NTuple{N, Symbol}` giving identifiers for the isotopes available in the yield table.
"""
function isotopes(::AbstractYield) end
"""
    remnant_mass(table::AbstractYield, Z, M)
Returns the remnant mass of a star with initial metallicity `Z` and initial stellar mass `M` (masses in solar masses).
"""
function remnant_mass(::AbstractYield, Z, M) end
"""
    remnant_mass(table::AbstractYield, Z, M)
Returns the pre-supernova mass of a star with initial metallicity `Z` and initial stellar mass `M` (masses in solar masses). This is not always equal to the stellar initial mass as massive stars can lose mass due to winds before becoming supernovae.
"""
function preSN_mass(::AbstractYield, Z, M) end
"""
    ejecta_mass(table::AbstractYield, Z, M)
Returns the mass of all ejected materials for a star with initial metallicity `Z` and initial stellar mass `M` (masses in solar masses). This is generally `preSN_mass(...) - remnant_mass(...)`. 
"""
function ejecta_mass(table::AbstractYield, Z, M) end
"""
    ejecta_metal_mass(table::AbstractYield, Z, M)
Returns the mass of ejected metals (elements heavier than helium) for a star with initial metallicity `Z` and initial stellar mass `M` (masses in solar masses).
"""
function ejecta_metal_mass(table::AbstractYield, Z, M) end
"""
    ejecta_alpha_mass(table::AbstractYield, Z, M)
Returns the mass of ejected alpha elements (O, Ne, Mg, Si, S, Ar, Ca, and Ti) for a star with initial metallicity `Z` and initial stellar mass `M` (masses in solar masses).
"""
function ejecta_alpha_mass(table::AbstractYield, Z, M) end

include("utilities.jl")
include("Nomoto2006.jl")
using .Nomoto2006
include("Kobayashi2006.jl")
using .Kobayashi2006
include("Portinari1998.jl")
using .Portinari1998
include("Nomoto1997.jl")
using .Nomoto1997

# Load submodules
include(joinpath("lifetimes", "lifetimes.jl"))
using .Lifetimes
include(joinpath("delay_time_distributions", "dtds.jl"))
using .DelayTimeDistributions
include(joinpath("AGB", "agbs.jl"))
using .AGB

export inverse # extended methods from other packages
# API generics
export isotopes, preSN_mass, remnant_mass, ejecta_mass, ejecta_metal_mass, ejecta_alpha_mass 
# Other functions
export N_TO
# Model types
export Nomoto1997SNIa, Nomoto2006SN, Kobayashi2006SN, Portinari1998SN, Vincenzo2016, Portinari1998Lifetimes, Karakas2010AGB 

end
