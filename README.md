# ScatteringKMatrix.jl

K-matrix amplitudes for coupled-channel scattering and production.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mmikhasenko/ScatteringKMatrix.jl")
```

## Usage

```julia
using ScatteringKMatrix
using StaticArrays

channels = SVector(
    TwoBodyChannel(1.1, 1.1),
    TwoBodyChannel(1.3, 1.3),
)
K = KMatrix([(M = 4.3, gs = [2.1, 0.0])])
T = TMatrix(K, channels)
A = amplitude(T, 5.0)
```

Two-body phase space:

- `TwoBodyChannel(m1, m2; L=0)` uses the square-root formula.
- `TwoBodyChewMandelstamChannel(m1, m2; L=0)` uses the Chew–Mandelstam function, with a vanishing real part at threshold and an imaginary part that approaches 1 at high energy.

Only `L = 0` is implemented. `ProductionAmplitude` builds the production vector on top of a `TMatrix`. Quasi-two-body channels and a gridded dispersive continuation are available as `QuasiTwoBodyChannel` and `InterpolatedChannel`.

For pole searches, wrap either square-root or Chew–Mandelstam channels in
`ContinuationChannel`. The channel's phase-space function has branch points
in the complex mass plane and can be analytically continued to multiple
sheets. The wrapper's `mode` selects which sheet is evaluated at each mass.
The sheet reached between the first and second thresholds is

```julia
channels_cm = SVector(
    TwoBodyChewMandelstamChannel(1.1, 1.1),
    TwoBodyChewMandelstamChannel(1.3, 1.3),
)
continued = continue_channels(channels_cm, 2.3)
T_unphysical = TMatrix(K, continued)
```

`continue_channels` accepts `TwoBodyChannel` and
`TwoBodyChewMandelstamChannel` entries. Each channel open at the reference
mass is evaluated on sheet I in the upper half-plane and sheet II in the lower
half-plane; channels that are still closed stay on sheet I.
`ContinuationChannel(channel, 2)` keeps the second sheet on both sides of the
real axis. Mode `-90` is the special case of a cut running straight down from
threshold. For `TwoBodyChannel`, that mode is also the branch convention used
by the unwrapped square-root phase space.

`AngledCutChannel(channel, α)` places that cut at an arbitrary angle in
radians, with the same convention as `angle`. `0` stays on sheet I, the same
as mode `1`. Negative angles drop the cut into the lower half-plane:
`-π/6` is 30° below the positive real axis, and `-π/2` reproduces mode
`-90`. Sheet II fills the wedge between the real axis and that ray.
A channel that stays on sheet I uses angle `0`, so the vector keeps one element type.

```julia
rotated = SVector(
    AngledCutChannel(channels_cm[1], -π / 6),
    AngledCutChannel(channels_cm[2], 0),
)
```

See the [continuation heatmap example](docs/README.md) for a side-by-side
visualization of the square-root and Chew–Mandelstam sheets and cuts.

## Tests

From a clone of this repository:

```julia
] test
```

## Notebooks

`notebooks/DD1_pipsi.jl` is a Pluto notebook for πJ/ψ scattering with a D-meson subchannel.

```julia
julia --project=notebooks -e 'using Pkg; Pkg.instantiate(); using Pluto; Pluto.run()'
```

Open `notebooks/DD1_pipsi.jl` from the Pluto file browser.

## License

MIT
