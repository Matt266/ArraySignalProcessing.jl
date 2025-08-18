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
function dml(pa::AbstractPhasedArray, Rxx, DoAs, f, c=c_0; 
            optimizer = Optim.LBFGS(), steer_kwargs=(), problem_kwargs=(), solve_kwargs=())
    p = pa, Rxx, f, c
    dml_cost = function(angles, p)
        pa, Rxx, f, c = p

        A = steer(pa, angles, f, c; steer_kwargs...)
        #PA = A*inv(A'*A)*A'
        PA = A*pinv(A) #TODO: check why the pinv() call here throws an error with CuArrys

        cost =  tr((I-PA)*Rxx)
        return real(cost)
    end

    f = OptimizationFunction(dml_cost, Optimization.AutoZygote())
    p = OptimizationProblem(f, DoAs, p; problem_kwargs...)
    s = solve(p, optimizer; solve_kwargs...)
    return s.u
end