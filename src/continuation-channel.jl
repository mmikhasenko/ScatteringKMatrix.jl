"""
    ContinuationChannel(channel, where=1)

Wrap a `TwoBodyChewMandelstamChannel` and select its analytic continuation.
The supported values of `where` are:

- `1`: the first-sheet Chew--Mandelstam function;
- `2`: the second-sheet function throughout the complex plane;
- `12`: sheet I for `imag(m) >= 0` and sheet II for `imag(m) < 0`;
- `90`: as for `12`, but only cross the right-hand cut when
  `real(m) > threshold(channel)`.

Mode `12` is useful for evaluating a chosen unphysical sheet around a pole.
Mode `90` displays the physical gluing of the sheets across the right-hand
unitarity cut.

Sheet II is the first sheet plus twice the phase-space factor built from
principal square roots. Modes `12` and `90` switch sheets in the complex
mass plane and assume real particle masses, so the unitarity cut lies on
the real axis. Mode `90` crosses only the positive-mass cut,
`real(m) > threshold(channel)`.
"""
struct ContinuationChannel{C<:TwoBodyChewMandelstamChannel} <: AbstractChannel
    nominal::C
    where::Int

    function ContinuationChannel(
        nominal::C,
        where::Integer=1,
    ) where {C<:TwoBodyChewMandelstamChannel}
        where in (1, 2, 12, 90) || throw(ArgumentError(
            "expected continuation mode 1, 2, 12, or 90; got $where",
        ))
        return new{C}(nominal, Int(where))
    end
end

threshold(ch::ContinuationChannel) = threshold(ch.nominal)

# The discontinuity is expressed with principal square roots.  Keeping this
# separate from `iρ(::TwoBodyChannel, m)`, whose branch prescription is chosen
# for direct physical-axis evaluation, is essential away from the real axis.
function _iρ_minimal(ch::TwoBodyChannel, m)
    return (
        1im *
        sqrt(m - (ch.m1 + ch.m2)) *
        sqrt(m + (ch.m1 + ch.m2)) *
        sqrt(m - (ch.m1 - ch.m2)) *
        sqrt(m + (ch.m1 - ch.m2)) /
        m^2
    )
end

function iρ(ch::ContinuationChannel, m)
    iρ_I = iρ(ch.nominal, m)
    ch.where == 1 && return iρ_I

    minimal = TwoBodyChannel(ch.nominal)
    iρ_II = iρ_I + 2 * _iρ_minimal(minimal, m)
    ch.where == 2 && return iρ_II
    ch.where == 12 && return imag(m) < 0 ? iρ_II : iρ_I
    return real(m) > threshold(ch) && imag(m) < 0 ? iρ_II : iρ_I
end

"""
    continue_channels(channels, m_reference; mode=12)

Wrap Chew--Mandelstam channels and continue those whose threshold is at or
below `m_reference`. Channels above that reference mass remain on sheet I.
This labels the conventional coupled-channel sheet associated with a real
mass interval between adjacent thresholds.
"""
function continue_channels(channels::SVector{N}, m_reference; mode::Integer=12) where {N}
    all(ch -> ch isa TwoBodyChewMandelstamChannel, channels) || throw(ArgumentError(
        "continue_channels requires TwoBodyChewMandelstamChannel entries",
    ))
    continued = map(channels) do ch
        selected_mode = real(m_reference) >= threshold(ch) ? mode : 1
        ContinuationChannel(ch, selected_mode)
    end
    return SVector{N}(continued)
end
