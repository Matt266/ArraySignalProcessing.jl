"""
wsf(pa::AbstractPhasedArray, Rxx, d, DoAs, f; fs=nothing, c=c_0, optimizer=NelderMead(), maxiters=1e3)

DoA estimation using weighted subspace fitting (WSF) .

arguments:
----------
    pa: AbstractPhasedArray to calculate the wsf estimation for
    Rxx: covariance matrix of the array which is used for estimation
    DoAs: vector/matrix of initial DoAs as starting point for WSF
    f: center/operating frequency
    c: propagation speed of the wave
    optimizer: used optimizer to solve the WSF problem
    maxiters: maximum optimization iterations

References:
-----------
M. Viberg and B. Ottersten, ‘Sensor array processing based on subspace fitting’, IEEE Trans. Signal Process., vol. 39, no. 5, pp. 1110–1121, May 1991.

B. Ottersten and M. Viberg, ‘Analysis of subspace fitting based methods for sensor array processing’, in International Conference on Acoustics, Speech, and Signal Processing, Glasgow, UK, 2003.

H. Krim and M. Viberg, ‘Two decades of array signal processing research: the parametric approach’, IEEE Signal Process. Mag., vol. 13, no. 4, pp. 67–94, Jul. 1996.

M. Pesavento, M. Trinh-Hoang, and M. Viberg, ‘Three More Decades in Array Signal Processing Research: An optimization and structure exploitation perspective’, IEEE Signal Process. Mag., vol. 40, no. 4, pp. 92–106, Jun. 2023.
"""
function wsf(pa::AbstractPhasedArray, Rxx, DoAs, f; c=c_0, coords=:azel, optimizer=:prima, maxiters=1e3)
    p = pa, Rxx, f, c
    wsf_cost = function(angles, p)
        pa, Rxx, f, c = p
        d = size(angles, 2)
        Λ, U = eigen(Rxx, sortby= λ -> -abs(λ))

        Λs = diagm(Λ[1:d])
        Λn = diagm(Λ[d+1:length(Λ)])

        Us = U[:, 1:d]
        Un = U[:, d+1:size(U, 2)]

        # estimate noise variance as mean of noise eigenvalues
        # optimal weights as in Ottersten and Viberg ‘Analysis of subspace fitting based methods for sensor array processing’
        σ² = mean(Λn)
        Λest = Λs - σ²*I
        W = Λest^2*inv(Λs)

        #A = hcat(steerphi.(Ref(pa), f, ϕ; fs=fs, c=c, direction=Incoming)...)
        A = steer(pa, f, angles; c=c, coords=coords)
        #PA = A*inv(A'*A)*A'
        PA = A*pinv(A)

        cost =  tr((I-PA)*Us*W*Us')
        return real(cost)
    end

    if optimizer == :prima
        # prima does not require parameters 
        # but starting points must be a vector
        shape = size(DoAs)
        function obj_func(angles)
            return wsf_cost(reshape(angles, shape), p)
        end
        result, _ = prima(obj_func, vec(DoAs))
        return reshape(result, shape)
    else
        # e.g., optimizer=NelderMead()
        f = OptimizationFunction(wsf_cost, Optimization.AutoForwardDiff())
        p = OptimizationProblem(f, DoAs, p)
        s = solve(p, optimizer; maxiters=maxiters)
        return s.u
    end
end