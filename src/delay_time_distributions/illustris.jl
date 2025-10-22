"""Implementation of the Illustris SN-Ia delay time distribution; see section 2.3.3 of Vogelsberger2013. SN-Ia rates based on Maoz+2012 and SN-II rates are found by integrating the Chabrier+2003 IMF"""
struct Illustris_SNIa <: AbstractDelayTimeDistribution end

function rate(::Illustris_SNIa, t)
    if t < 40e6
        return 0.0 * t
    else
        return 1.3e-3 * (t / 40e6)^(-1.12) * (0.12 / 40e6)
    end
end

function expectation(::Illustris_SNIa, t1, t2, m)
    @argcheck 0 <= t1 <= t2
    if t2 < 40e6
        return 0.0 * t1
    else
        t1 = max(t1, 40e6) # Rate is 0 for t < 40 Myr
        p2 = t2 * (t2 / 40e6)^(-1.12)
        p1 = t1 * (t1 / 40e6)^(-1.12)
        return m * 1.3e-3 * (0.12 / 40e6) * (p2 - p1) / (1 - 1.12)
    end
end