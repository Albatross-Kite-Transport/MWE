using ModelingToolkit, OrdinaryDiffEq, Plots
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: t_nounits as t, D_nounits as D


@mtkmodel PWM begin
    # @extend v, i = oneport = OnePort()
    @parameters begin
        T = 0.1
        duty = 0.5
        Vcc = 5
    end
    @components begin
        output = RealOutput(u_start=Vcc)
        trigger = Modulo()
    end
    @equations begin
        D(output.u) ~ 0
        connect(t, T, )
    end
    @continuous_events begin
        (t ~ T*duty) => [output.u ~ 0]
        (t ~ T) => [output.u ~ Vcc]
    end
end

@mtkmodel PWM_test begin
    @parameters begin
        R = 1.0
        V = 5.0
    end
    @components begin
        resistor = Resistor(R = R)
        VDD = Voltage()
        pwm = PWM(Vcc = V)
        ground = Ground()
    end
    @equations begin
        connect(pwm.output, VDD.V)
        connect(VDD.p, resistor.p)
        connect(ground.g, VDD.n, resistor.n)
    end

end
@mtkbuild sys = PWM_test()
unknowns(sys)
prob = ODEProblem(sys, [sys.pwm.output.u => 0.0], (0, 1.0); dt=0.001)
sol = solve(prob, Tsit5())

plot(sol, idxs = [sys.resistor.i],
    title = "Circuit Demonstration")