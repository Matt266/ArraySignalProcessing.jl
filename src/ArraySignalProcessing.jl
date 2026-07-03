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
using Adapt
using FillArrays

export  IsotropicArray, IsotropicArrayManifold, NestedArray, NestedArrayManifold, SampledArray, SampledArrayManifold,
        TappedDelayLine, TappedDelayLineManifold, AzEl, WaveVec, SlowVec, RAzEl,
        dsb_weights, bartlett, mvdr_weights, mpdr_weights, capon_weights, capon, lcmv_weights, lcmp_weights,
        lcmv_gsc, lcmv, GeneralizedSidelobeCanceler, process, process!, LMS, RLS,
        whitenoise, diffnoise, esprit, music, unitary_esprit, lasso, λ_stable, omp, ols, bpdn,
        aic, mdl, wsf, dml, sml, unconditional_signals, find_doas, classical_crb, classical_crb_conditional_signal, classical_crb_gaussian_signal_unknown_power,
        classical_crb_gaussian_signal_unknown_spectrum, classical_crb_gaussian_signal_known_spectrum

c_0 = 299792458.0

include("core/core.jl")
include("beamforming/beamforming.jl")
include("direction_finding/direction_finding.jl")
include("detection/detection.jl")
include("utility/utility.jl")
end