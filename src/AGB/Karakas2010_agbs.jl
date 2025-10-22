
const _Karakas2010_Zs = (0.0001, 0.004, 0.008, 0.02)
const _Karakas2010_M = (1.0, 1.25, 1.5, 1.75, 1.9, 2.0, 2.1, 2.25, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)

struct Karakas2010AGBEntry{A, B <: AbstractArray{A}, C, D}
    Z::A
    M::B
    Mfinal::B
    table::C
    isotopes::D
end

function Karakas2010AGBEntry(fname)
    lines = readlines(fname)
    parts = split.(lines)
    M = unique([parse(Float64, p[1]) for p in parts])
    Z = parse(Float64, parts[1][2]) # unique([parse(Float64, p[2]) for p in parts])
    M_final = unique([parse(Float64, p[3]) for p in parts])
    isotopes = unique([String(p[4]) for p in parts])
    prepend!(isotopes, ["Mfinal"])
    # Mass of 6.5 is available only for Z=0.02;
    # we will remove this model for uniform M sampling
    mass_mask = M .≈ 6.5
    if any(mass_mask)
        M = M[.~mass_mask]
        M_final = M_final[.~mass_mask]
    end
    # isotopes = Vector{String}(undef, length(lines))
    # M = Vector{Float64}(undef, length(lines))
    # table = Matrix{Float64}(undef, length(parts[1]) - 4, length(M))
    table = Array{Float64}(undef, length(parts[1]) - 4, length(isotopes), length(M))
    # table = fill(Inf, (length(parts[1]) - 4, length(isotopes), length(M)))
    mask = vcat(3, 6:lastindex(parts[1]))
    for (i, p) in enumerate(parts) # enumerate(lines)
        Mi = parse(Float64, p[1])
        if !isapprox(Mi, 6.5)
            i = findfirst(x -> p[4] == x, isotopes)
            j = findfirst(x -> Mi == x, M)
            table[:, i, j] = parse.(Float64, replace.(p[mask], '−' => '-'))
            table[:, 1, j] .= M_final[j]
        end
    end
    return Karakas2010AGBEntry(Z, M, M_final, table, isotopes)
end

"""
    Karakas2010AGB(; bounds=Interpolations.Throw()) <: AbstractYield
Load the Karakas+2010 AGB yield table. The yield table can be interpolated by calling it with the metal mass fraction `Z` and stellar mass `M` (in solar masses) of the progenitor. The keyword argument `bounds` should be a valid `Interpolations.jl` extrapolation specifier that will determine how the interpolation is extrapolated (e.g., `Flat()`).

```jldoctest
julia> k = Karakas2010AGB();

julia> result = k(0.02, 4.5);

julia> result isa NamedTuple
true

julia> isapprox(result.Mfinal, 0.852)
true

julia> isapprox(result.g, 3.5865554e-7)
true

julia> result2 = k(0.01, 4.75);

julia> isapprox(result2.Mfinal, 0.8721666666666666)
true

julia> isapprox(result2.p, -0.2422697775)
true

julia> remnant_mass(k, 0.01, 4.75) == k(0.01, 4.75).Mfinal
true

julia> ejecta_mass(k, 0.01, 4.75) == 4.75 - remnant_mass(k, 0.01, 4.75)
true
```
"""
struct Karakas2010AGB{I, B} <: AbstractYield
    itp::B
end

function Karakas2010AGB(; bounds=Throw())

    # Read tablea1
    
    # Read tables a2-a5
    # a6 two models with partial mixing zones that
    # we don't want to use
    entries = Karakas2010AGBEntry.(joinpath(@__DIR__, "..", "data", "Karakas2010", "tablea"*i*".txt") for i in string.(2:5))
    # return entries
    # The entries[1].table is 3-D with shape (length(properties), length(isotopes), length(masses))
    # The properties are different forms of the yield / wind properties, for now we only care about the total mass lost in the wind, so that's all we'll interpolate

    iso_mat = [SVector{length(entries[1].isotopes)}(i.table[2, :, j]) for i=entries, j=eachindex(entries[1].M)]
    # The number of model stars (M) at each metallicity is the same coming out of Karakas2010Entry, but the values of the initial stellar masses are not the same -- some metallicities have 2.0 models and some have 2.1. We will pre-interpolate the values so that at all metallicities, model stars with the same masses are available.
    iso_mat = [interpolate((entries[i].M,), iso_mat[i, :], Gridded(Linear()))(m) for i=eachindex(entries), m=_Karakas2010_M]
    # Need to sort entries according to Z, not sorted by default
    sidx = sortperm([e.Z for e in entries])
    iso_itp = interpolate((SVector(_Karakas2010_Zs), SVector(_Karakas2010_M)), iso_mat[sidx, :], Gridded(Linear()))
    bounds = extend_bounds(bounds, length(entries[1].isotopes))
    iso_itp = extrapolate(iso_itp, bounds)
    isotopes = Tuple(Symbol.(entries[1].isotopes))
    return Karakas2010AGB{isotopes, typeof(iso_itp)}(iso_itp)
end

isotopes(::Karakas2010AGB{I}) where I = I
_nt(::Karakas2010AGB{I}) where I = NamedTuple{I, NTuple{length(I), Float64}}
function (x::Karakas2010AGB)(Z, M)
    return _nt(x)(Tuple(x.itp(Z, M)))
end
remnant_mass(x::Karakas2010AGB, Z, M) = x(Z, M).Mfinal
ejecta_mass(x::Karakas2010AGB, Z, M) = M - remnant_mass(x, Z, M)

# Don't care to add these right now
# function filter_metals(nt::NamedTuple)
#     names = (:Li6, :Li7, :Be9, :B10, :B11, :C12, :C13, :N14, :N15, :O16, :O17, :O18, :F19, :Ne20, :Ne21, :Ne22, :Na23, :Mg24, :Mg25, :Mg26, :Al26, :Al27, :Si28, :Si29, :Si30, :P31, :S32, :S33, :S34, :S36, :Cl35, :Cl37, :Ar36, :Ar38, :Ar40, :K39, :K40, :K41, :Ca40, :Ca42, :Ca43, :Ca44, :Ca46, :Ca48, :Sc45, :Ti46, :Ti47, :Ti48, :Ti49, :Ti50, :V50, :V51, :Cr50, :Cr52, :Cr53, :Cr54, :Mn55, :Fe54, :Fe56, :Fe57, :Fe58, :Co59, :Ni58, :Ni60, :Ni61, :Ni62, :Ni64, :Cu63, :Cu65, :Zn64, :Zn66, :Zn67, :Zn68, :Zn70, :Ga69, :Ga71)
#     return nt[names]
# end
# ejecta_metal_mass(x::Nomoto2006SN, Z, M) = sum(values(filter_metals(x(Z, M))))

# # Filter to include only alpha elements
# function filter_alpha(nt::NamedTuple)
#     names = (:O16, :O17, :O18, :Ne20, :Ne21, :Ne22, :Mg24, :Mg25, :Mg26, :Si28, :Si29, :Si30, :S32, :S33, :S34, :S36, :Ar36, :Ar38, :Ar40, :Ca40, :Ca42, :Ca43, :Ca44, :Ca46, :Ca48, :Ti46, :Ti47, :Ti48, :Ti49, :Ti50)
#     return nt[names]
# end
# ejecta_alpha_mass(x::Nomoto2006SN, Z, M) = sum(values(filter_alpha(x(Z, M))))