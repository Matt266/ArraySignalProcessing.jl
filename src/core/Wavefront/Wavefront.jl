abstract type Wavefront end
abstract type PlaneWave <: Wavefront end
abstract type SphericalWave <: Wavefront end

include("AzEl.jl")
include("WaveVec.jl")
include("SlowVec.jl")