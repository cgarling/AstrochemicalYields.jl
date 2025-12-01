# Code for interpolating the stellar lifetime fits of Vincenzo+2016
# to the PARSEC stellar models
# DOI 10.1093/mnras/stw1145 

"""
    Vincenzo2016(; bounds=Interpolations.Throw())
Interpolator for the stellar lifetime fits of Vincenzo+2016. These are based on PARSEC stellar models. The lifetimes can be interpolated by calling an instance with a metal mass fraction `Z` and a stellar mass `M` (in solar masses), returning the stellar lifetime in Gyr. The keyword argument `bounds` should be a valid `Interpolations.jl` extrapolation specifier that will determine how the interpolation is extrapolated (e.g., `Flat()`).

```jldoctest vincenzo
julia> using AstrochemicalYields: Vincenzo2016

julia> v = Vincenzo2016();

julia> M, Z = 1.0, 1e-2;

julia> v(Z, M) ≈ 9.876997213870718
true
```

The inverse operation (getting the stellar mass of stars just finishing their lives) can be performed by calling `inverse(v::Vincenzo2016)(Z, t)` where `t` is in Gyr.

```jldoctest vincenzo
julia> using AstrochemicalYields: inverse

julia> inverse(v)(Z, 9.876997213870718) ≈ M
true
```

Derivatives are also available. `AstrochemicalYields.Lifetimes.derivative(v::Vincenzo2016, Z)` returns a function that takes an argument `M` (the initial stellar mass (in solar masses) and returns the derivative of the stellar lifetime (in Gyr) with respect to the initial stellar mass evaluated at `M`.

```jldoctest vincenzo
julia> using AstrochemicalYields.Lifetimes: derivative

julia> d = derivative(v, Z);

julia> isapprox(d(M), -35.05809527556551)
true
```

The derivative of the inverse option (the derivative of the stellar mass of stars ending their lives with respect to time) is also available.

```jldoctest vincenzo
julia> dinv = derivative(inverse, v, Z); # dinv is a function taking a time (in Gyr) argument

julia> isapprox(dinv(9.0), -0.0325150150483193)
true
```
"""
struct Vincenzo2016{A}
    itp::A
end

function (v::Vincenzo2016)(Z, M)
    A, B, C = v.itp(Z)
    # Equation 1 in Vincenzo2016
    return A * exp(B * M^(-C))
end

function derivative(v::Vincenzo2016, Z)
    A, B, C = v.itp(Z)
    return M -> - A * B * C * exp(B * M^(-C)) * M^(-C - 1)
end

function inverse(v::Vincenzo2016)
    return (Z, t) -> begin
        A, B, C = v.itp(Z)
        (log(t / A) / B)^(-1 / C)
    end
end

function derivative(::typeof(inverse), v::Vincenzo2016, Z)
    A, B, C = v.itp(Z)
    return t -> - (log(t / A) / B)^(-1 / C - 1) / B / C / t
end

function Vincenzo2016(; bounds=Throw())
    fname = joinpath(@__DIR__, "..", "data", "supplementary_material_vincenzo2016.dat")
    lines = readlines(fname)
    table = Matrix{Float64}(undef, length(lines)-1, 4)
    for (i, line) in enumerate(lines[2:end])
        parts = split(line)
        table[i,:] .= parse.(Float64, parts[1:4])
    end
    table_vec = [SVector{size(table,2)-1}(table[i,2:end]) for i in 1:size(table,1)]
    Z_data = SVector{size(table, 1)}(table[:,1])
    Z = 1e-4:1e-4:0.0299
    if !isapprox(Z, Z_data)
        error("Vincenzo2016 data malformed.")
    end
    # itp = interpolate((Z,), table_vec, Gridded(Linear()))
    itp = scale(interpolate(table_vec, BSpline(Cubic())), Z)
    itp = extrapolate(itp, bounds)
    return Vincenzo2016(itp)
end