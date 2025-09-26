module AstrochemicalYields

"""
`AbstractYield` is the abstract supertype for all yield tables. Yield table subtypes should be made callable with initial metal mass fraction `Z` and mass `M` (in solar masses), returning the yield for all isotopes in units of solar masses. Subtypes should additionally implement the following methods:
 - [`isotopes`](@ref)
 - [`remnant_mass`](@ref)
 - [`ejecta_mass`](@ref)
 - [`ejecta_metal_mass`](@ref)
 - [`ejecta_alpha_mass`](@ref)
"""
abstract type AbstractYield end

"""List of α elements used in ejecta_alpha_mass taken from Nomoto+2006. C is sometimes included, but not here."""
const α_elements = ("O", "Ne", "Mg", "Si", "S", "Ar", "Ca", "Ti")

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

include("Nomoto2006.jl")
using .Nomoto2006
include("Kobayashi2006.jl")
using .Kobayashi2006
include("Portinari1998.jl")
using .Portinari1998
include("Nomoto1997.jl")
using .Nomoto1997

export isotopes, preSN_mass, remnant_mass, ejecta_mass, ejecta_metal_mass, ejecta_alpha_mass # Generics
export Nomoto1997SNIa, Nomoto2006SN, Kobayashi2006SN, Portinari1998SN # specific implementations

end
