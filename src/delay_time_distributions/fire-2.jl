"""Implementation of the FIRE-2 SN-Ia rate; see Appendix A of Hopkins+2018."""
struct FIRE2_SNIa <: AbstractDelayTimeDistribution end

function rate(::FIRE2_SNIa, t::T) where T
    if t < 37.53e6 # 37.53 Myr
        return 0.0 * t
    else
        return 5.3e-14 + 1.6e-11 * exp(-((t - 50e6) / 1e7)^2 / 2)
    end
end

function expectation(::FIRE2_SNIa, t1, t2, m)
    @argcheck 0 <= t1 <= t2
    # Deal with section where rate is 0
    if t2 < 37.53e6
        return 0.0 * t1
    else
        t1 = max(t1, 37.53e6)
    end
    return m * (5.3e-14 * (t2 - t1) + 8e-5 * sqrt2π * (erf((5 - t1/1e7) / sqrt2) - erf((5 - t2/1e7) / sqrt2)))
end

struct FIRE2_SNII <: AbstractDelayTimeDistribution end

function rate(::FIRE2_SNII, t)
    if (t < 3.4e6) || (t > 37.53e6)
        return 0.0 * t
    elseif 3.4e6 <= t <= 10.37e6
        return 5.408e-10
    elseif 10.37e6 <= t <= 37.53e6
        return 2.516e-10
    end
end

function expectation(::FIRE2_SNII, t1, t2, m)
    @argcheck 0 <= t1 <= t2
    tmin = max(t1, 3.4e6)
    tmax = min(t2, 37.53e6)
    if (tmin >= tmax) || (tmax <= 3.4e6)
        return 0.0 * t1
    end
    # Compute overlap with each constant-rate interval
    dt1 = max(0.0, min(tmax, 10.37e6) - max(tmin, 3.4e6))
    dt2 = max(0.0, min(tmax, 37.53e6) - max(tmin, 10.37e6))
    return m * (5.408e-10 * dt1 + 2.516e-10 * dt2)
end