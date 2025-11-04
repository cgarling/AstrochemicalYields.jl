"""Module containing models for AGB yields."""
module AGB

using ..AstrochemicalYields: AbstractYield, extend_bounds
import ..AstrochemicalYields: isotopes, remnant_mass, ejecta_mass, ejecta_metal_mass, ejecta_alpha_mass

include("Karakas2010_agbs.jl")
using .Karakas2010

export Karakas2010AGB

end # module