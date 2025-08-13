module ArraySignalProcessing

using LinearAlgebra
using StatsBase
using Random
using ForwardDiff
using Roots
using Convex
using SCS
using Optimization
using OptimizationOptimJL
using Optim
using Interpolations
using PRIMA
using Peaks
using ProximalAlgorithms
using ProximalOperators
using ProximalCore
using Zygote
using Enzyme
#using DifferentiationInterface: AutoEnzyme

export PhasedArray, IsotropicArray, PhasedArray, NestedArray, steer, IsotropicArrayManifold, SampledArrayManifold,
        AzEl, WaveVec,
        dsb_weights, bartlett, mvdr_weights, mpdr_weights, capon_weights, capon,
        whitenoise, diffnoise, esprit, music, unitary_esprit, lasso, λ_stable, omp, ols, bpdn,
        aic, mdl, wsf, dml, sml, unconditional_signals, find_doas

c_0 = 299792458.0

include("core/core.jl")
include("beamforming/beamforming.jl")
include("direction_finding/direction_finding.jl")
include("detection/detection.jl")
include("utility/utility.jl")
end