# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffness
# for l < l_0), n tether segments, tether drag and reel-in and reel-out. 
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D, setp, getu, parameter_index, SymbolicIndexingInterface
using NonlinearSolve
using ControlPlots

struct Point
    type::Symbol
    position::Union{Vector{Float64}, Nothing}
    velocity::Union{Vector{Float64}, Nothing}
    mass::Union{Float64, Nothing}
    force::Vector{Float64}
end

struct Segment
    points::Tuple{Int, Int}
    l0::Union{Float64, Nothing}
    stiffness::Float64
    damping::Float64
end

struct Pulley
    segments::Tuple{Int, Int}
    sum_length::Float64
end

"""
Add or remove points, segments and pulleys to these lists to configure your system
"""
# protune-speed-system.jpg
points = [
    Point(:fixed,  [0, 0, 0],  zeros(3), nothing, zeros(3)),  # Fixed point
    Point(:quasi_static, [-1, 0, 0], zeros(3), 1.0, zeros(3)),
    Point(:quasi_static, [-2, 0, 0], zeros(3), 1.0, zeros(3)),
    Point(:quasi_static, [-3, 0, 0], zeros(3), 1.0, zeros(3)),
    Point(:dynamic, [-4, 0, 0], zeros(3), 100.0, [0., 0., 0.]),
]

stiffness = 614600
damping = 4730
segments = [
    Segment((1, 2), norm(points[1].position - points[2].position), stiffness, damping),
    Segment((2, 3), norm(points[2].position - points[3].position), stiffness, damping),
    Segment((3, 4), norm(points[3].position - points[4].position), stiffness, damping),
    Segment((4, 5), norm(points[4].position - points[5].position), stiffness, damping),
]

pulleys = []

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
    duration = 5.0                                # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
end

function set_tether_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    se.d_tether = d
    se.c_spring = c_spring_4mm * (d/4.0)^2
    se.damping = damping_4mm * (d/4.0)^2
end
                              
function calc_initial_state(points, segments, pulleys)
    POS0 = zeros(3, length(points))
    VEL0 = zeros(3, length(points))
    L0 = zeros(length(pulleys))
    V0 = zeros(length(pulleys))
    for i in eachindex(points)
        POS0[:, i] .= points[i].position
        VEL0[:, i] .= points[i].velocity
    end
    for i in eachindex(pulleys)
        L0[i] = segments[pulleys[i].segments[1]].l0
        V0[i] = 0.0
    end
    POS0, VEL0, L0, V0
end

function calc_spring_forces(pos::AbstractMatrix{T}, vel, pulley_l0) where T
    # loop over all segments to calculate spring forces
    spring_force = zeros(T, length(segments))
    spring_force_vec = zeros(T, 3, length(segments))
    segment = zeros(T, 3)
    unit_vector = zeros(T, 3)
    rel_vel = zeros(T, 3)
    for (tether_idx, tether) in enumerate(segments)
        found = false
        for (pulley_idx, pulley) in enumerate(pulleys)
            if tether_idx == pulley.segments[1] # each tether should only be part of one pulley
                l0 = pulley_l0[pulley_idx]
                found = true
                break
            elseif tether_idx == pulley.segments[2]
                l0 = pulley.sum_length - pulley_l0[pulley_idx]
                found = true
                break
            end
        end
        if !found
            l0 = tether.l0
        end
        p1, p2 = tether.points[1], tether.points[2]

        segment         .= pos[:, p2] .- pos[:, p1]
        len              = norm(segment)
        unit_vector     .= segment ./ len
        rel_vel         .= vel[:, p1] .- vel[:, p2]
        spring_vel       = rel_vel ⋅ unit_vector
        spring_force[tether_idx]          = (tether.stiffness * tether.l0 * (len - l0) - tether.damping * tether.l0 * spring_vel)
        spring_force_vec[:, tether_idx]  .= spring_force[tether_idx] .* unit_vector
    end
    return spring_force_vec, spring_force
end

function calc_acc(se::Settings3, pos::AbstractMatrix{T}, vel, pulley_l0) where T
    spring_force_vec, spring_force = calc_spring_forces(pos, vel, pulley_l0)
    
    pulley_acc = zeros(T, length(pulleys))
    for (pulley_idx, pulley) in enumerate(pulleys)
        M = 3.1
        pulley_force = spring_force[pulley.segments[1]] - spring_force[pulley.segments[2]]
        pulley_acc[pulley_idx] = pulley_force / M
    end

    acc = zeros(T, 3, length(points))
    force = zeros(T, 3)
    for (point_idx, point) in enumerate(points)
        if point.type === :fixed
            acc[:, point_idx] .= 0.0
        else
            force .= 0.0
            for (j, tether) in enumerate(segments)
                if point_idx in tether.points
                    inverted = tether.points[2] == point_idx
                    if inverted
                        force .-= spring_force_vec[:, j]
                    else
                        force .+= spring_force_vec[:, j]
                    end
                end
            end
            force .+= point.force
            acc[:, point_idx] .= force ./ point.mass .+ se.g_earth
        end
    end
    return acc, pulley_acc
end
@register_symbolic calc_acc(se::Settings3, pos, vel, pulley_l0)

function model(se::Settings3)
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

        segment_vec(t)[1:3, eachindex(segments)]
        unit_vector(t)[1:3, eachindex(segments)]
        l_spring(t), c_spring(t), damping(t), m_tether_particle(t)
        len(t)[eachindex(segments)]
        l0(t)[eachindex(segments)]
        rel_vel(t)[1:3, eachindex(segments)]
        spring_vel(t)[eachindex(segments)]
        spring_force(t)[eachindex(segments)]
        spring_force_vec(t)[1:3, eachindex(segments)] # spring force from spring p1 to spring p2
    end

    POS0, VEL0, L0, V0 = calc_initial_state(points, segments, pulleys)

    defaults = Pair{Num, Real}[]
    guesses = Pair{Num, Real}[]
    eqs = [
        D(pulley_l0) ~ pulley_vel
        D(pulley_vel) ~ pulley_acc - 10pulley_vel
    ]

    for (point_idx, point) in enumerate(points)
        if point.type === :fixed
            eqs = [
                eqs
                pos[:, point_idx] ~ point.position
                vel[:, point_idx] ~ zeros(3)
            ]
        elseif point.type === :dynamic
            eqs = [
                eqs
                D(pos[:, point_idx]) ~ vel[:, point_idx]
                D(vel[:, point_idx]) ~ acc[:, point_idx]
            ]
            defaults = [
                defaults
                [pos[j, point_idx] => POS0[j, point_idx] for j in 1:3]
                [vel[j, point_idx] => 0 for j in 1:3]
            ]
        elseif point.type === :quasi_static
            eqs = [
                eqs
                acc[:, point_idx] ~ zeros(3)
                vel[:, point_idx] ~ zeros(3)
            ]
            guesses = [
                guesses
                [pos[j, point_idx] => POS0[j, point_idx] for j in 1:3]
                [vel[j, point_idx] => 0 for j in 1:3]
            ]
        else
            println("wrong type")
        end
    end

    defaults = [
        defaults
        pulley_l0 => L0
        pulley_vel => V0
    ]

    eqs = [
        eqs
        vec(acc) .~ vec(calc_acc(se, pos, vel, pulley_l0)[1])
        vec(pulley_acc) .~ vec(calc_acc(se, pos, vel, pulley_l0)[2])
    ]
    
    eqs = reduce(vcat, Symbolics.scalarize.(eqs))
    @named sys = ODESystem(eqs, t)
    
    @time sys = structural_simplify(sys; simplify=false)
    sys, pos, vel, defaults, guesses
end

function simulate(se, sys, defaults, guesses)
    dt = 0.1
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    @time prob = ODEProblem(sys, defaults, tspan; guesses, simplify=false)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=AutoFiniteDiff()); dt, abstol=tol, reltol=tol, saveat=ts)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=AutoFiniteDiff()); dt, abstol=tol, reltol=tol, saveat=ts)
    sol, elapsed_time
end

function play(se, sol, pos)
    dt = 0.1
    ylim = (-10.5, 2.5)
    xlim = (-6.0, 6.0)
    mkpath("video")
    z_max = 0.0
    # text position
    xy = (se.l0/4.2, z_max-7)
    start = time_ns()
    i = 1; j = 0
    for i in eachindex(sol.t)
        segs = [[segments[k].points[1], segments[k].points[2]] for k in eachindex(segments)]
        pos_ = [sol[pos][i][:, k] for k in eachindex(sol[pos][i][1, :])]
        plot2d(pos_, segs, sol.t[i]; xlim, ylim, xy)
        if se.save
            ControlPlots.plt.savefig("video/"*"img-"*lpad(j, 4, "0"))
        end
        j += 1
        wait_until(start + 0.1 * 0.5 * sol.t[i] * 1e9)
    end
    if se.save
        println("Run the script ./bin/export_gif to create the gif file!")
    end
    nothing
end

function main()
    se = Settings3()
    set_tether_diameter!(se, se.d_tether) # adapt spring and damping constants to tether diameter
    sys, pos, vel, defaults, guesses = model(se)
    sol, elapsed_time = simulate(se, sys, defaults, guesses)
    play(se, sol, pos)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    println("Number of evaluations per step: ", round(sol.stats.nf/(se.duration/0.1), digits=1))
    sol, pos, vel, sys
end

if (! @isdefined __BENCH__) || __BENCH__ == false
    sol, pos, vel, sys = main()
end
__BENCH__ = false
nothing
