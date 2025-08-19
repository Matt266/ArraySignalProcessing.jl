abstract type Wavefront end
abstract type PlaneWave <: Wavefront end
abstract type SphericalWave <: Wavefront end

include("AzEl.jl")
include("WaveVec.jl")
include("SlowVec.jl")
include("RAzEl.jl")

function Base.length(wf::Wavefront)
    # each column is a new coordinate
    # -> length is number of coordinate points
    return size(wf.coords, 2)
end