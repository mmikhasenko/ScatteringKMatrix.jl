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

For pole searches, wrap Chew–Mandelstam channels in `ContinuationChannel`.
The sheet reached between the first and second thresholds is

```julia
channels_cm = SVector(
    TwoBodyChewMandelstamChannel(1.1, 1.1),
    TwoBodyChewMandelstamChannel(1.3, 1.3),
)
continued = continue_channels(channels_cm, 2.3)
T_unphysical = TMatrix(K, continued)
```

`continue_channels` accepts only `TwoBodyChewMandelstamChannel` entries. Each
channel open at the reference mass is evaluated on sheet I in the upper
half-plane and sheet II in the lower half-plane; channels that are still
closed stay on sheet I. `ContinuationChannel(channel, 2)` keeps the second
sheet on both sides of the real axis. Mode `90` crosses to sheet II only
through the right-hand cut, where `real(m) > threshold`.

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
