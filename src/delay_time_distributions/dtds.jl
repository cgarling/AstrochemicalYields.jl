
"""Module implementing delay time distributions (DTDs)."""
module DelayTimeDistributions

using ..AstrochemicalYields: Vincenzo2016, inverse

using ArgCheck: @argcheck, @check
using Distributions: pdf, cdf, mean
using IrrationalConstants: sqrtπ, sqrt2, sqrt2π
using QuadGK: quadgk
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

"""
    N_SNII(imf, t1, t2, TO_mass_func, mmin, mmax)
Returns the expected number of SN-II per solar mass of stars formed in the time interval between `t1` and `t2`.

# Arguments
 - `imf` is an initial mass function supporting Distributions.jl functions `pdf(imf, m)` or `cdf(imf, m)` for initial stellar mass `m`; see InitialMassFunctions.jl for implementations.
 - `t1` is the starting time of the interval in years.
 - `t2` is the ending time of the inverval in years.
 - `TO_mass_func` is a function that takes a time in years (like `t1`, `t2`) and returns the mass (in solar masses) of stars leaving the main sequence.
 - `mmin` is the lower bound on the initial stellar mass for stars that will end their lives as SN-II (typically ~8 solar masses).
 - `mmax` is the upper bound on the initial stellar mass for stars that will end their lives as SN-II (often ~100 solar masses).

    N_SNII(imf, t1, t2, v::Vincenzo2016, Z, mmin, mmax)
This signature uses the [`Vincenzo2016`](@ref) stellar lifetimes to derive the `TO_mass_func` for provided metallicity `Z`.
"""
function N_SNII(imf, t1, t2, TO_mass_func, mmin, mmax)
    @argcheck t2 > t1
    mmin = max(mmin, TO_mass_func(t2))
    mmax = min(mmax, TO_mass_func(t1))
    if mmin > mmax
        return 0.0 * t1
    end
    try
        return (cdf(imf, mmax) - cdf(imf, mmin)) / mean(imf)
    catch
        return quadgk(x -> pdf(imf, x), mmin, mmax)[1] / quadgk(x -> x * pdf(imf, x), 0.08, mmax)[1]
    end
end
function N_SNII(imf, t1, t2, v::Vincenzo2016, Z, mmin, mmax)
    TO_mass_func(x) = inverse(v)(Z, x / 1e9)
    return N_SNII(imf, t1, t2, TO_mass_func, mmin, mmax)
end
# function dN_SNII(imf, t, v::Vincenzo2016, Z, mmin, mmax)
#     # TO_mass_func(x) = inverse(v)(Z, x / 1e9)
#     dm = derivative(inverse, v, Z)
#     TO_mass = inverse(v)(Z, TO_mass_func(t / 1e9))
#     r = 0.0
#     p = pdf(imf, TO_mass)
#     if TO_mass < mmax
#         r += p * dm(TO_mass)
#     end
#     if TO_mass > mmin
#         r += 0
#     end
#     return r
# end

# Load individual DTD implementations
include("fire-2.jl")
include("illustris.jl")

end # module