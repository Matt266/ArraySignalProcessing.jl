#TODO Refactor Naming and Dispatch
# There are several layers for naming:
#   1. paramater model: classical, bayesian, hybrid
#   2. observation vector: gaussian, ... (currently now non-gaussian use case)
#   3. source signal model: unkown nonrandom (conditional), unkown random (conditional) -> mostly gaussian source signal
# Maybe there is some better way to encapsulate that in a AbstractSignalModel type to dispatch on that?
# Then, extensions like mutual-coupling matrix C (x = CAx+n) could also be handled in a unified structure
include("classical_crb/classical_crb.jl")
include("bayesian_crb/bayesian_crb.jl")
include("hybrid_crb/hybrid_crb.jl")