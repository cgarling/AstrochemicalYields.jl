# Code for interpolating the stellar lifetime fits of Vincenzo+2016
# to the PARSEC stellar models
# DOI 10.1093/mnras/stw1145 

"""
    Vincenzo2016(; bounds=Throw())
Interpolator for the stellar lifetime fits of Vincenzo+2016. These are based on PARSEC stellar models. The lifetimes can be interpolated by calling an instance with a metal mass fraction `Z` and a stellar mass `M` (in solar masses), returning the stellar lifetime in Gyr. The keyword argument `bounds` should be a valid `Interpolations.jl` extrapolation specifier that will determine how the interpolation is extrapolated (e.g., `Flat()`).

```jldoctest
julia> using AstrochemicalYields: Vincenzo2016

julia> v = Vincenzo2016();

julia> M, Z = 1.0, 1e-2;

julia> v(Z, M)
9.876997213870718
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

function Vincenzo2016(; bounds=Throw())
    fname = joinpath(@__DIR__, "data", "supplementary_material_vincenzo2016.dat")
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