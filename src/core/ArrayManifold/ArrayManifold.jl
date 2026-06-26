include("AbstractArrayManifold.jl")
include("NestedArrayManifold.jl")
include("IsotropicArrayManifold.jl")
include("TappedDelayLineManifold.jl")
include("SampledArrayManifold.jl")

const IsotropicArray = IsotropicArrayManifold
const NestedArray = NestedArrayManifold
const SampledArray = SampledArrayManifold
const TappedDelayLine = TappedDelayLineManifold