using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: D_nounits as D, t_nounits as t, setu, setp, getu, getp

Ts = 0.1
@mtkmodel Pendulum begin
    @parameters begin
        g = 9.8
        L = 0.4
        K = 1.2
        m = 0.3
        τ = 0.0 # input
    end
    @variables begin
        θ(t) = 0.0 # state
        ω(t) = 0.0 # state
        y(t) # output
    end
    @equations begin
        D(θ)    ~ ω
        D(ω)    ~ -g/L*sin(θ) - K/m*ω + τ/m/L^2
        y       ~ θ * 180 / π
    end
end

@mtkbuild mtk_model = Pendulum()
prob = ODEProblem(mtk_model, nothing, (0.0, Ts); initializealg=CheckInit())
integ = init(prob, Tsit5())
set_x = setu(mtk_model, Initial.(unknowns(mtk_model)))
get_x = getu(mtk_model, unknowns(mtk_model))
set_u = setp(mtk_model, [mtk_model.τ])
get_u = getp(mtk_model, [mtk_model.τ])
get_h = getu(mtk_model, [mtk_model.y])
p = (integ, set_x, set_u, get_h)

function f!(xnext, x, u, _, p)
    (integ, set_x, set_u, _) = p
    set_x(integ, x)
    set_u(integ, u)
    reinit!(integ; reinit_dae=true)
    step!(integ)
    xnext .= integ.u # sol.u is the state, called x in the function
    nothing
end

xnext = zeros(2)

println("Initial")
f!(xnext, zeros(2), 0.0, nothing, p)
@show xnext

println("Just changing state x - doesn't work, xnext stays at zero")
f!(xnext, ones(2), 0.0, nothing, p)
@show xnext

println("Just changing input u - does work, xnext is non-zero")
f!(xnext, zeros(2), 1.0, nothing, p)
@show xnext
nothing
