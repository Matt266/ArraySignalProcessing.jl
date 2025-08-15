abstract type AbstractUpdateMethod end
include("LMS.jl")
include("RLS.jl")

abstract type AbstractBeamformer end
include("gsc.jl")
include("dsb.jl")
include("mvdr.jl")
include("lcmv.jl")