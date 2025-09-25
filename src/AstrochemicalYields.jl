module AstrochemicalYields

using DataInterpolationsND: NDInterpolation, LinearInterpolationDimension, BSplineInterpolationDimension
using Interpolations: interpolate, Linear, Gridded, extrapolate, Flat, Throw
using Printf: @sprintf
using StaticArrays: SVector

include("Nomoto2006.jl")

end
