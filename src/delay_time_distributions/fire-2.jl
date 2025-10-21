
struct FIRE2_SNIa <: AbstractDelayTimeDistribution end

function pdf(::FIRE2_SNIa, t::Real)
    if t < 37.53e6 # 37.53 Myr
        return zero(t)
    else
        return 5.3e-11 + 1.6e-8 * exp10(-((t - 50e6) / 10)^2 / 2)
    end
end