
"""Module implementing delay time distributions (DTDs)."""
module DelayTimeDistributions

using ..AstrochemicalYields: Vincenzo2016, inverse

using ArgCheck: @argcheck, @check
using IrrationalConstants: sqrtπ, sqrt2, sqrt2π
using SpecialFunctions: erf

# """Abstract type for delay time distributions (DTDs) and event rates. DTDs cannot be used to infer event rates (e.g., supernova rates) alone and must combine other information. This type can also be used to parameterize event rate models directly (which are common in simulations, for example). Not all methods will be available for event rate models."""
"""Abstract type for delay time distributions (DTDs)."""
abstract type AbstractDelayTimeDistribution end

"""
    rate(::AbstractDelayTimeDistribution, t)
Returns the occurence rate of events per year per solar mass of stars (``dN/dt``) formed at a delay time `t` (in years) after a star formation episode.
"""
rate(::AbstractDelayTimeDistribution, t)

"""
    expectation(::AbstractDelayTimeDistribution, t1, t2, m)
Returns the expectation value for the total number of events occuring for a stellar population of initial mass `m` (in solar masses) between times `t1` and `t2` (in years) after a star formation episode.
"""
function expectation(d::AbstractDelayTimeDistribution, t1, t2, m)
    @argcheck 0 <= t1 <= t2
    return m * quadgk(Base.Fix1(rate, d), t1, t2)[1]
end

# Load individual DTD implementations
include("fire2_dtd.jl")
include("illustris_dtd.jl")

end # module