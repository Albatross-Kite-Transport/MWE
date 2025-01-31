# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffness
# for l < l_0), n tether segments, tether drag and reel-in and reel-out. 
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D
using ControlPlots


# Define user input structures (keep the same)
struct Point
    fixed::Bool
    position::Union{Vector{Float64}, Nothing}
    velocity::Union{Vector{Float64}, Nothing}
    mass::Union{Float64, Nothing}
    force::Vector{Float64}
end

struct Tether
    points::Tuple{Int, Int}
    l0::Union{Float64, Nothing}  # Rest length for spring
    stiffness::Float64  # New: spring stiffness
    damping::Float64    # New: damping coefficient
end

struct Pulley
    tethers::Tuple{Int, Int}
    sum_length::Float64
end

# Example user input
points = [
    Point(true,  zeros(3), zeros(3), nothing, zeros(3)),  # Fixed point
    Point(false, [0, 0, -1], zeros(3), 1.0, [0, 0, -9.81]), # Pulley
    Point(true,  [-5.0, 0.0, -1.0], zeros(3), 1.0, [0., 0., -9.81]),  # Fixed point
    Point(false, [5.0, 0.0, -1.0], zeros(3), 1.1, [0., 0., -9.81]),  # Movable point
]

tethers = [
    Tether((1, 2), 1.0, 1000.0, 5.0),
    Tether((2, 3), 5.0, 1000.0, 5.0),
    Tether((2, 4), 5.0, 1000.0, 5.0),
]

pulleys = [
    Pulley((2, 3), 10.0)
]

"""
Acc on point is:
if fixed:
    0
else:
    F/m

Force on point is:
    - sum of minus spring forces on tether-connected points
    - if tether-connected point is pulley: keep direction of pulley, but take minus spring force from the other point connected to the pulley

l0 on pulley:
    - stretched_length_1 / (stretched_length_1 + stretched_length_2) * sum_length
"""


@with_kw mutable struct Settings3 @deftype Float64
    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    v_wind_tether::Vector{Float64} = [2, 0.0, 0.0]
    rho = 1.225
    cd_tether = 0.958
    l0 = 50                                      # initial tether length             [m]
    v_ro = 2                                     # reel-out speed                  [m/s]
    d_tether = 4                                 # tether diameter                  [mm]
    rho_tether = 724                             # density of Dyneema            [kg/m³]
    c_spring = 614600                            # unit spring constant              [N]
    rel_compression_stiffness = 0.01             # relative compression stiffness    [-]
    damping = 473                                # unit damping constant            [Ns]
    segments::Int64 = 2                          # number of tether segments         [-]
    α0 = π/10                                    # initial tether angle            [rad]
    duration = 20.0                                # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
end

function set_tether_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    se.d_tether = d
    se.c_spring = c_spring_4mm * (d/4.0)^2
    se.damping = damping_4mm * (d/4.0)^2
end
                              
function calc_initial_state(points, tethers, pulleys)
    POS0 = zeros(3, length(points))
    VEL0 = zeros(3, length(points))
    L0 = zeros(length(pulleys))
    V0 = zeros(length(pulleys))
    for i in eachindex(points)
        POS0[:, i] .= points[i].position
        VEL0[:, i] .= points[i].velocity
    end
    for i in eachindex(pulleys)
        L0[i] = tethers[pulleys[i].tethers[1]].l0
        V0[i] = 0.0
    end
    @show POS0, VEL0, L0, V0
    POS0, VEL0, L0, V0
end

function model(se)
    @parameters c_spring0=se.c_spring/(se.l0/se.segments) l_seg=se.l0/se.segments
    @parameters rel_compression_stiffness = se.rel_compression_stiffness
    @variables begin 
        pos(t)[1:3, eachindex(points)]
        vel(t)[1:3, eachindex(points)]
        acc(t)[1:3, eachindex(points)]
        force(t)[1:3, eachindex(points)]

        pulley_force(t)[eachindex(pulleys)]
        pulley_l0(t)[eachindex(pulleys)] # first tether length in pulley
        pulley_vel(t)[eachindex(pulleys)]
        pulley_acc(t)[eachindex(pulleys)]

        segment(t)[1:3, eachindex(tethers)]
        unit_vector(t)[1:3, eachindex(tethers)]
        l_spring(t), c_spring(t), damping(t), m_tether_particle(t)
        len(t)[eachindex(tethers)]
        l0(t)[eachindex(tethers)]
        rel_vel(t)[1:3, eachindex(tethers)]
        spring_vel(t)[eachindex(tethers)]
        spring_force(t)[eachindex(tethers)]
        spring_force_vec(t)[1:3, eachindex(tethers)] # spring force from spring p1 to spring p2
    end
    # basic differential equations
    eqs = [
        vec(D(pos)) ~ vec(vel)
        vec(D(vel)) ~ vec(acc)
        D(pulley_l0) ~ pulley_vel
        D(pulley_vel) ~ pulley_acc
    ]

    for (pulley_idx, pulley) in enumerate(pulleys)
        M = 3.1
        eqs = [
            eqs
            pulley_force[pulley_idx]    ~ spring_force[pulley.tethers[1]] - spring_force[pulley.tethers[2]]
            pulley_acc[pulley_idx]      ~ pulley_force[pulley_idx] / M
        ]
    end

    # loop over all tethers to calculate spring forces
    for (tether_idx, tether) in enumerate(tethers)
        found = false
        for (pulley_idx, pulley) in enumerate(pulleys)
            if tether_idx == pulley.tethers[1] # each tether should only be part of one pulley
                eqs = [
                    eqs
                    l0[tether_idx] ~ pulley_l0[pulley_idx]
                ]
                found = true
                break
            elseif tether_idx == pulley.tethers[2]
                eqs = [
                    eqs
                    l0[tether_idx] ~ pulley.sum_length - pulley_l0[pulley_idx]
                ]
                found = true
                break
            end
        end
        if !found
            eqs = [
                eqs
                l0[tether_idx] ~ tether.l0
            ]
        end
        p1, p2 = tether.points[1], tether.points[2]
        eqs = [
            eqs
            segment[:, tether_idx]       ~ pos[:, p2] - pos[:, p1]
            len[tether_idx]              ~ norm(segment[:, tether_idx])
            unit_vector[:, tether_idx]   ~ segment[:, tether_idx]/len[tether_idx]
            rel_vel[:, tether_idx]       ~ vel[:, p2] - vel[:, p1]
            spring_vel[tether_idx]       ~ rel_vel[:, tether_idx] ⋅ unit_vector[:, tether_idx]
            spring_force[tether_idx]     ~ (tether.stiffness * (len[tether_idx] - l0[tether_idx]) + tether.damping * spring_vel[tether_idx])
            spring_force_vec[:, tether_idx]  ~ spring_force[tether_idx] * unit_vector[:, tether_idx]
        ]
    end

    for (i, point) in enumerate(points)
        if point.fixed
            eqs = [
                eqs
                force[:, i]  ~ zeros(3)
                acc[:, i]    ~ zeros(3)
            ]
        else
            # tether - inverted
            f::Vector{Num} = zeros(Num, 3)
            for (j, tether) in enumerate(tethers)
                if i in tether.points
                    inverted = tether.points[2] == i
                    if inverted
                        f .-= spring_force_vec[:, j]
                    else
                        f .+= spring_force_vec[:, j]
                    end
                end
            end
            eqs = [
                eqs
                force[:, i]  ~ f .- 1 * vel[:, i]
                acc[:, i]    ~ force[:, i] / point.mass + se.g_earth
            ]
        end
    end

    eqs = reduce(vcat, Symbolics.scalarize.(eqs))
    @named sys = ODESystem(eqs, t)
    sys = structural_simplify(sys)
    sys, pos, vel
end

function simulate(se, sys)
    dt = 0.1
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    POS0, VEL0, L0, V0 = calc_initial_state(points, tethers, pulleys)
    prob = ODEProblem(sys, [pos => POS0, vel => VEL0, sys.pulley_l0 => L0, sys.pulley_vel => V0], tspan)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=true); dt, abstol=tol, reltol=tol, saveat=ts)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=true); dt, abstol=tol, reltol=tol, saveat=ts)
    sol, elapsed_time
end

function play(se, sol, pos)
    dt = 0.1
    ylim = (-10.5, 1.5)
    xlim = (-6.0, 6.0)
    mkpath("video")
    z_max = 0.0
    # text position
    xy = (se.l0/4.2, z_max-7)
    start = time_ns()
    i = 1; j = 0
    for time in 0:dt:se.duration
        # while we run the simulation in steps of 20ms, we update the plot only every 150ms
        # therefore we have to skip some steps of the result
        while sol.t[i] < time
            i += 1
        end
        plot2d(sol[pos][i], time; segments=length(tethers), xlim, ylim, xy)
        if se.save
            ControlPlots.plt.savefig("video/"*"img-"*lpad(j, 4, "0"))
        end
        j += 1
        wait_until(start + 0.5 * time * 1e9)
    end
    if se.save
        println("Run the script ./bin/export_gif to create the gif file!")
    end
    nothing
end

function main()
    global sol, pos, vel, len, c_spr
    se = Settings3()
    set_tether_diameter!(se, se.d_tether) # adapt spring and damping constants to tether diameter
    sys, pos, vel = model(se)
    sol, elapsed_time = simulate(se, sys)
    play(se, sol, pos)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    println("Number of evaluations per step: ", round(sol.stats.nf/(se.duration/0.02), digits=1))
    sol, pos, vel, sys
end

if (! @isdefined __BENCH__) || __BENCH__ == false
    sol, pos, vel, sys = main()
end
__BENCH__ = false
nothing
