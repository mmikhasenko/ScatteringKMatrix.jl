module ScatteringKMatrix

using HadronicLineshapes
using Interpolations
using LinearAlgebra
using StaticArrays
using Parameters
using QuadGK

export TwoBodyChannel, iρ
export TwoBodyChewMandelstamChannel
export ChewMandestam
export threshold
include("two-body-channel.jl")

export real_ρ
export nominal_threshold
export QuasiTwoBodyChannel, QuasiTwoBodyChannelBW
include("quasi-two-body-channel.jl")

export InterpolatedChannel
include("interpolation.jl")

export KMatrix, amplitude, npoles, nchannels
export TMatrix, DMatrix, detD, channels
include("t-matrix.jl")

export ProductionAmplitude
export production_pole, production_nonpole
include("production-amplitude.jl")

end # module ScatteringKMatrix
