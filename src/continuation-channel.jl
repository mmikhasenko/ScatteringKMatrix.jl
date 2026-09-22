"""
    ContinuationChannel(channel, mode=1)

Wrap a `TwoBodyChannel` or `TwoBodyChewMandelstamChannel` and select its
analytic continuation.
The supported values of `mode` are:

- `1`: the first-sheet phase-space function;
- `2`: the second-sheet function throughout the complex plane;
- `12`: sheet I for `imag(m) >= 0` and sheet II for `imag(m) < 0`;
- `-90`: as for `12`, but only cross the right-hand cut when
  `real(m) > threshold(channel)`. The cut runs downward from threshold.

Mode `12` is useful for evaluating a chosen unphysical sheet around a pole.
Mode `-90` displays the physical gluing of the sheets across the right-hand
unitarity cut.

For square-root phase space, sheet II is the negative of sheet I. For the
Chew--Mandelstam function, it is sheet I plus twice the phase-space factor
built from principal square roots. Modes `12` and `-90` are the cut angles
`-π` and `-π/2`: sheet II fills the wedge from the positive real axis down
to that ray. `AngledCutChannel` takes any other angle, in radians.
"""
struct ContinuationChannel{C<:AbstractTwoBodyChannel} <: AbstractChannel
    nominal::C
    mode::Int

    function ContinuationChannel(
        nominal::C,
        mode::Integer=1,
    ) where {C<:AbstractTwoBodyChannel}
        mode in (1, 2, 12, -90) || throw(ArgumentError(
            "expected continuation mode 1, 2, 12, or -90; got $mode",
        ))
        return new{C}(nominal, Int(mode))
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

# `TwoBodyChannel` historically exposes the square-root phase space with its
# cut running down from threshold.  Undo the sheet switch in that wedge to
# recover the conventional first sheet used by the continuation wrappers.
function _iρ_first_sheet(ch::TwoBodyChannel, m)
    iρ_downward_cut = iρ(ch, m)
    return real(m) > threshold(ch) && imag(m) < 0 ? -iρ_downward_cut : iρ_downward_cut
end

_iρ_first_sheet(ch::TwoBodyChewMandelstamChannel, m) = iρ(ch, m)

function _iρ_second_sheet(ch::TwoBodyChannel, m)
    return -_iρ_first_sheet(ch, m)
end

function _iρ_second_sheet(ch::TwoBodyChewMandelstamChannel, m)
    return iρ(ch, m) + 2 * _iρ_minimal(TwoBodyChannel(ch), m)
end

function iρ(ch::ContinuationChannel, m)
    iρ_I = _iρ_first_sheet(ch.nominal, m)
    ch.mode == 1 && return iρ_I

    iρ_II = _iρ_second_sheet(ch.nominal, m)
    ch.mode == 2 && return iρ_II
    ch.mode == 12 && return imag(m) < 0 ? iρ_II : iρ_I
    # Mode -90: sheet II only below the real axis and above threshold.
    ch.mode == -90 && return real(m) > threshold(ch) && imag(m) < 0 ? iρ_II : iρ_I
    error("unsupported continuation mode $(ch.mode)")
end

"""
    AngledCutChannel(channel, cut_angle)

Two-body phase-space continuation with the unitarity cut rotated to `cut_angle`
radians. The direction is the same convention as Julia's `angle`: `0` lies
on the positive real axis, and negative angles drop the cut into the lower
half-plane. `cut_angle` must lie in `[-π, π]`.

Sheet II fills the wedge between the positive real axis and the ray that
starts at `threshold(channel)` and runs in the direction `cut_angle`.
The wedge is empty at `0`, so that angle stays on sheet I and reproduces
`ContinuationChannel` mode `1`. The values `-π/2` and `-π` reproduce modes
`-90` and `12`.
"""
struct AngledCutChannel{C<:AbstractTwoBodyChannel} <: AbstractChannel
    nominal::C
    cut_angle::Float64

    function AngledCutChannel(
        nominal::C,
        cut_angle::Real,
    ) where {C<:AbstractTwoBodyChannel}
        α = Float64(cut_angle)
        -π <= α <= π || throw(ArgumentError(
            "cut angle must lie in [-π, π] radians; got $cut_angle",
        ))
        return new{C}(nominal, α)
    end
end

threshold(ch::AngledCutChannel) = threshold(ch.nominal)

function iρ(ch::AngledCutChannel, m)
    # Both angles are in radians. `0` is sheet I, and `-π/2` is mode -90.
    ϕ = angle(m - threshold(ch))
    α = ch.cut_angle
    on_II = α < 0 ? (α < ϕ < 0) : (0 < ϕ < α)
    return on_II ? _iρ_second_sheet(ch.nominal, m) : _iρ_first_sheet(ch.nominal, m)
end

"""
    continue_channels(channels, m_reference; mode=12)

Wrap two-body square-root or Chew--Mandelstam channels and continue those
whose threshold is at or below `m_reference`. Channels above that reference
mass remain on sheet I.
This labels the conventional coupled-channel sheet associated with a real
mass interval between adjacent thresholds.

Every entry is a `ContinuationChannel`. For a rotated cut, build an
`SVector` of `AngledCutChannel` directly so the element type stays concrete.
"""
function continue_channels(channels::SVector{N}, m_reference; mode::Integer=12) where {N}
    all(ch -> ch isa AbstractTwoBodyChannel, channels) || throw(ArgumentError(
        "continue_channels requires two-body phase-space channels",
    ))
    continued = map(channels) do ch
        selected_mode = real(m_reference) >= threshold(ch) ? mode : 1
        ContinuationChannel(ch, selected_mode)
    end
    return SVector{N}(continued)
end
