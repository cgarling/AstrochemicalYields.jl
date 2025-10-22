
module AGB

using ..AstrochemicalYields: AbstractYield, extend_bounds
import ..AstrochemicalYields: isotopes, remnant_mass, ejecta_mass, ejecta_metal_mass, ejecta_alpha_mass


using Interpolations: interpolate, Throw, Gridded, Linear, extrapolate
using StaticArrays: SVector

include("Karakas2010_agbs.jl")

export Karakas2010AGB

end # module