"""
dml(pa::AbstractPhasedArray, Rxx, DoAs, f; c=c_0, coords=:azel, optimizer=NelderMead(), maxiters=1e3)

DoA estimation using Deterministic Maximum Likelihood (DML).

arguments:
----------
    pa: AbstractPhasedArray to calculate the dml estimator for
    Rxx: covariance matrix of the array which is used for estimation
    DoAs: vector/matrix of initial DoAs as starting point
    f: center/operating frequency
    c: propagation speed of the wave
    optimizer: used optimizer to solve the DML problem
    maxiters: maximum optimization iterations

References:
-----------
H. Krim and M. Viberg, ‘Two decades of array signal processing research: the parametric approach’, IEEE Signal Process. Mag., vol. 13, no. 4, pp. 67–94, Jul. 1996.
"""
function dml(pa::AbstractPhasedArray, Rxx, DoAs, f; c=c_0, coords=:azel, optimizer=:prima, maxiters=1e3)
    p = pa, Rxx, f, c
    dml_cost = function(angles, p)
        pa, Rxx, f, c = p

        A = steer(pa, f, angles; c=c, coords=coords)
        #PA = A*inv(A'*A)*A'
        PA = A*pinv(A)

        cost =  tr((I-PA)*Rxx)
        return real(cost)
    end

    if optimizer == :prima
        # prima does not require parameters 
        # but starting points must be a vector
        shape = size(DoAs)
        function obj_func(angles)
            return dml_cost(reshape(angles, shape), p)
        end
        result, _ = prima(obj_func, vec(DoAs))
        return reshape(result, shape)
    else
        # e.g., optimizer=NelderMead()
        f = OptimizationFunction(dml_cost, Optimization.AutoForwardDiff())
        p = OptimizationProblem(f, DoAs, p)
        s = solve(p, optimizer; maxiters=maxiters)
        return s.u
    end
end