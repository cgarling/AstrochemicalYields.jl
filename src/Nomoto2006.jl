# Code for SN / stellar wind yields

struct Nomoto2006Entry{A, B <: AbstractArray{A}, C, D}
    Z::A
    M::B
    E::B
    Mcut::B
    # p::B
    # d::B
    table::C
    isotopes::D
end
# function Nomoto2006Entry(args...)
#     T = promote_type(typeof(args[1]), map(typeof ∘ eltype, args[2:end])...)
#     return Nomoto2006Entry(convert(T, args[1]), )
# end
const _Nomoto2006_Zs = (0.0, 0.001, 0.004, 0.02)
const _Nomoto2006_SN_M = (13.0, 15.0, 18.0, 20.0, 25.0, 30.0, 40.0)
function Nomoto2006Entry(Z::Number)
    if !any(map(Base.Fix1(isapprox, Z), _Nomoto2006_Zs))
        throw(ArgumentError("Argument `Z` must be approximately equal to one of $(_Nomoto2006_Zs)."))
    end
    fname = joinpath(@__DIR__, "data", "Nomoto2006", "SN_Z_"*@sprintf("%03d", round(Int, Z*1000))*".csv")
    lines = readlines(fname)
    Z = parse(Float64, split(lines[1])[end])
    M = parse.(Float64, split(lines[2])[2:end])
    E = parse.(Float64, split(lines[3])[2:end])
    Mcut = parse.(Float64, split(lines[4])[2:end])
    isotopes = Vector{String}(undef, length(lines)-4)
    table = Matrix{Float64}(undef, length(isotopes), length(M))
    for (i, line) in enumerate(lines[5:end])
        parts = split(line)
        isotopes[i] = parts[1]
        table[i, :] = parse.(Float64, replace.(parts[2:end], '−' => '-'))
    end
    # return isotopes, table
    return Nomoto2006Entry(Z, M, E, Mcut, table, isotopes)
end


struct Nomoto2006SN{N, B, C, F <: Function}
    iso_itp::B
    mcut_itp::C
    isotopes::NTuple{N, Symbol}
    constructor::F
end
isotopes(x::Nomoto2006SN) = x.isotopes
# function (x::Nomoto2006SN{N})(Z, M) where N
#     return NamedTuple{x.isotopes, NTuple{N, Float64}}(x.iso_itp(Z, M))
# end
function (x::Nomoto2006SN{N})(Z, M) where N
    return x.constructor(x.iso_itp(Z, M))
end

# function Nomoto2006SN()
#     entries = Nomoto2006Entry.(_Nomoto2006_Zs)
#     mat = Array{Float64}(undef, length(entries), length(entries[1].M), length(entries[1].isotopes))
#     for i in axes(mat, 1)
#         for j in axes(mat, 2)
#             mat[i,j,:] .= entries[i].table[:,j]
#         end
#     end
#     # return mat
#     return NDInterpolation(mat, (LinearInterpolationDimension(collect(_Nomoto2006_Zs)),
#                            LinearInterpolationDimension(entries[1].M)))
# end

function Nomoto2006SN()
    entries = Nomoto2006Entry.(_Nomoto2006_Zs)
    return Nomoto2006SN(entries, eachindex(entries[1].isotopes))
end
function Nomoto2006SN(isotopes)
    isotopes = string.(isotopes)
    entries = Nomoto2006Entry.(_Nomoto2006_Zs)
    good = [findfirst(x -> i == x, entries[1].isotopes) for i in isotopes]
    return Nomoto2006SN(entries, good)
end

function Nomoto2006SN(entries, good)
    N = length(good)
    iso_mat = [SVector{N}(i.table[good, j]) for i=entries, j=eachindex(entries[1].M)]
    iso_itp = interpolate((SVector(_Nomoto2006_Zs), SVector(_Nomoto2006_SN_M)), iso_mat, Gridded(Linear()))
    iso_itp = extrapolate(iso_itp, Throw())
    mcut_mat = [i.Mcut[j] for i=entries, j=eachindex(entries[1].M)]
    mcut_itp = interpolate((SVector(_Nomoto2006_Zs), SVector(_Nomoto2006_SN_M)), mcut_mat, Gridded(Linear()))
    mcut_itp = extrapolate(mcut_itp, Throw())
    isotopes = Tuple(Symbol.(entries[1].isotopes))
    constructor(x) = NamedTuple{isotopes, NTuple{N, Float64}}(x)
    return Nomoto2006SN(iso_itp, mcut_itp, isotopes, constructor)
end
# Nomoto2006Entry(0.0, Float64[13, 15, 18, 20, 25, 30, 40], Float64[1, 1, 1, 1, 1, 1, 1],
#                 Float64[1.57, 1.48, 1.65, 1.66, 1.92, 2.07, 2.89], 
#                 Float64[6.59E+00, 7.58E+00, 8.43E+00, 8.77E+00, 1.06E+01, 1.17E+01, 1.40E+01],
#                 Float64[1.49E−16, 1.69E−16, 1.28E−16, 8.66E−17, 2.02E−16, 1.34E−16, 3.46E−16],
#                 [])
# function convert_data(fname::AbstractString)
#     lines = readlines(fname)
#     Z = parse(Float64, split(lines[1])[end])
#     M = parse.(Float64, split(lines[2])[2:end])
#     E = parse.(Float64, split(lines[3])[2:end])
#     Mcut = parse.(Float64, split(lines[4])[2:end])
#     isotopes = [split(line)[1] for line in lines[5:end]]
# end