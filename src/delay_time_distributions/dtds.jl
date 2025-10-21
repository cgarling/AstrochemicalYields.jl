
"""Module implementing delay time distributions (DTDs)."""
module DelayTimeDistributions
import Distributions: ContinuousUnivariateDistribution, pdf, mean
using QuadGK: quadgk

"""Abstract type for delay time distributions; a subtype of `Distributions.ContinuousUnivariateDistribution`."""
abstract type AbstractDelayTimeDistribution <: ContinuousUnivariateDistribution end

"""
    pdf(::AbstractDelayTimeDistribution, t::Real)
Returns the occurence rate of events per year per solar mass of stars (``dN/dt``) formed at a delay time `t` (in years) after a star formation episode.
"""
pdf(::AbstractDelayTimeDistribution, t::Real)

"""
    mean(::AbstractDelayTimeDistribution, t::Real, m::Real=one(t))
Returns the expectation value for the total number of events occuring for a stellar population of mass `m` (in solar masses) within a time `t` (in years) after a star formation episode.
"""
mean(d::AbstractDelayTimeDistribution, t::Real, m::Real=one(t)) = quadgk(x -> m * pdf(d, x), zero(t), t)[1]

# Load individual DTD implementations
include("fire-2.jl")

end # module