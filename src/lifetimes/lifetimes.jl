"""Module for models describing stellar lifetimes."""
module Lifetimes

using ..AstrochemicalYields: interp_log

using ArgCheck: @argcheck, @check
using Distributions: pdf, cdf, mean
using Interpolations: interpolate, scale, Linear, BSpline, Cubic, Gridded, extrapolate, Flat, Throw
import InverseFunctions: inverse
using QuadGK: quadgk
using StaticArrays: SVector, @SVector

include("Vincenzo2016.jl")
include("Portinari1998_lifetimes.jl")
include("rates.jl")

export Vincenzo2016, Portinari1998Lifetimes, N_TO

end # module