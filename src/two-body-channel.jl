abstract type AbstractChannel end

"""
    TwoBodyChannel(m1, m2; L=0)

Two-body channel representation with masses `m1`, `m2` and angular momentum `L`.

# Fields
- `m1::Complex{Float64}`: Mass of first particle
- `m2::Complex{Float64}`: Mass of second particle  
- `L::Int`: Angular momentum quantum number (currently only L=0 implemented)
"""
struct TwoBodyChannel <: AbstractChannel
    m1::Complex{Float64}
    m2::Complex{Float64}
    L::Int
end

TwoBodyChannel(m1, m2; L::Int=0) = TwoBodyChannel(m1, m2, L)

threshold(ch::TwoBodyChannel) = real(ch.m1 + ch.m2)

"""
    ρ(ch::TwoBodyChannel, m)

Calculate the phase space factor for a two-body channel at mass `m`.

# Arguments
- `ch::TwoBodyChannel`: Two-body channel definition
- `m::Number`: Mass where to evaluate phase space

# Returns
Phase space value ρ(m) for the given channel and mass
"""
function iρ(ch::TwoBodyChannel, m)
    ϕ = -π / 2
    ch.L != 0 && error("not implemented")
    1im *
    sqrt(cis(ϕ) * (m - (ch.m1 + ch.m2))) * cis(-ϕ / 2) *
    sqrt(m + (ch.m1 + ch.m2)) *
    sqrt((m^2 - (ch.m1 - ch.m2)^2)) /
    m^2
end

"""
    ChewMandestam(s, m1sq, m2sq)

Calculate the Chew-Mandelstam function for a two-body channel.

# Arguments
- `s::Number`: Squared center-of-mass energy
- `m1sq::Number`: Squared mass of first particle
- `m2sq::Number`: Squared mass of second particle

# Returns
Chew-Mandelstam function value (complex)
"""
function ChewMandestam(s, m1sq, m2sq)
    # Ensure s is complex for proper branch handling
    s_complex = Complex(s)
    m1, m2 = sqrt(Complex(m1sq)), sqrt(Complex(m2sq))
    #
    sth, spth = (m1 + m2)^2, (m1 - m2)^2
    λh = sqrt(spth - s_complex) * sqrt(sth - s_complex)
    #
    val = 1 / (π) * (
        λh / s_complex * log((m1sq + m2sq - s_complex + λh) / (2 * m1 * m2)) -
        (m1sq - m2sq) * (1.0 / s_complex - 1 / sth) * log(m1 / m2)
    )
    return val
end

"""
    TwoBodyChewMandelstamChannel(m1, m2; L=0)

Two-body channel representation using Chew-Mandelstam function for phase space calculation.

# Fields
- `m1::Complex{Float64}`: Mass of first particle
- `m2::Complex{Float64}`: Mass of second particle  
- `L::Int`: Angular momentum quantum number (currently only L=0 implemented)
"""
struct TwoBodyChewMandelstamChannel <: AbstractChannel
    m1::Complex{Float64}
    m2::Complex{Float64}
    L::Int
end

TwoBodyChewMandelstamChannel(m1, m2; L::Int=0) = TwoBodyChewMandelstamChannel(m1, m2, L)

threshold(ch::TwoBodyChewMandelstamChannel) = real(ch.m1 + ch.m2)

"""
    iρ(ch::TwoBodyChewMandelstamChannel, m)

Calculate the phase space factor for a two-body channel using Chew-Mandelstam function at mass `m`.

# Arguments
- `ch::TwoBodyChewMandelstamChannel`: Two-body channel definition
- `m::Number`: Mass where to evaluate phase space

# Returns
Phase space value ρ(m) for the given channel and mass
"""
function iρ(ch::TwoBodyChewMandelstamChannel, m)
    ch.L != 0 && error("not implemented")
    s = m^2
    m1sq = abs2(ch.m1)
    m2sq = abs2(ch.m2)
    return ChewMandestam(s, m1sq, m2sq)
end
