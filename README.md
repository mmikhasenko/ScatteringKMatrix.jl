# ScatteringKMatrix.jl

A Julia package for K-matrix formalism in scattering theory.

## Installation

To install the package, start Julia REPL and run:

```julia
using Pkg
Pkg.add(url="https://github.com/mmikhasenko/ScatteringKMatrix.jl")
```

## Testing

The package comes with a tests>. It is the first step to validate setup.
To run the test, execute in Julia REPL:

```julia
] test
```

## Interactive Examples

This package includes several Pluto notebooks. To run them

- Activate notebooks environment:

```julia
julia> ]
pkg> activate notebooks
pkg> instantiate
# backspace
```

- Start Pluto:
```julia
using Pluto
Pluto.run()
```

- Once Pluto opens in your browser, navigate to the `notebooks` folder in this repository and open any of the following:
   - `example.jl` - Basic introduction to K-matrix formalism
   - `DD1_piJpsi.jl` - Analysis of πJ/ψ scattering with D mesons

## Features

- Implementation of K-matrix formalism for multi-channel scattering
- Support for coupled-channel analysis
- Production amplitudes
- Quasi-two-body decay channels
- Phase space integration tools
- Chew-Mandelstam function for analytic phase space calculation

## Two-Body Channels

The package provides two constructors for two-body decay channels:

- `TwoBodyChannel(m1, m2; L=0)` - Uses direct square-root formula for phase space calculation
- `TwoBodyChewMandelstamChannel(m1, m2; L=0)` - Uses Chew-Mandelstam function for analytic phase space calculation

The Chew-Mandelstam implementation provides an analytic continuation of the phase space factor that is normalized such that the real part is zero at threshold and the imaginary part approaches 1 at high energy. This is useful for dispersion relations and K-matrix calculations.

## License

This package is licensed under the MIT License.
