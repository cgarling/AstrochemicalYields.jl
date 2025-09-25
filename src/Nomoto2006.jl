# Code for Nomoto2006 SN-II yields

const _Nomoto2006_Zs = (0.0, 0.001, 0.004, 0.02)
const _Nomoto2006_SN_M = (13.0, 15.0, 18.0, 20.0, 25.0, 30.0, 40.0)
struct Nomoto2006Entry{A, B <: AbstractArray{A}, C, D}
    Z::A
    M::B
    E::B
    table::C
    isotopes::D
end

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
    isotopes = Vector{String}(undef, length(lines)-4+1) # +1 for Mcut
    table = Matrix{Float64}(undef, length(isotopes), length(M))
    for (i, line) in enumerate(lines[5:end])
        parts = split(line)
        isotopes[i] = parts[1]
        table[i, :] = parse.(Float64, replace.(parts[2:end], '−' => '-'))
    end
    table[end, :] .= Mcut
    isotopes[end] = "Mcut"
    return Nomoto2006Entry(Z, M, E, table, isotopes)
end


struct Nomoto2006SN{I, B}
    itp::B
    # isotopes::NTuple{N, Symbol} # This is now part of the type signature (I) for performance
end
isotopes(::Nomoto2006SN{I}) where I = I
_nt(::Nomoto2006SN{I}) where I = NamedTuple{I, NTuple{length(I), Float64}}
function (x::Nomoto2006SN)(Z, M)
    return _nt(x)(Tuple(x.itp(Z, M)))
end

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
    isotopes = Tuple(Symbol.(entries[1].isotopes[good]))
    return Nomoto2006SN{isotopes, typeof(iso_itp)}(iso_itp)
end
