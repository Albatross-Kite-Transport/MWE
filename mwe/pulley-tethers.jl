# Tutorial example simulating a 3D mass-spring system with a nonlinear spring (1% stiffness
# for l < l_0), n tether segments, tether drag and reel-in and reel-out. 
using ModelingToolkit, OrdinaryDiffEq, LinearAlgebra, Timers, Parameters, ControlPlots
using ModelingToolkit: t_nounits as t, D_nounits as D
using NonlinearSolve
using ControlPlots


struct Point
    type::Symbol
    position::Union{Vector{Float64}, Nothing}
    velocity::Union{Vector{Float64}, Nothing}
    mass::Union{Float64, Nothing}
    force::Vector{Float64}
end

struct Tether
    points::Tuple{Int, Int}
    l0::Union{Float64, Nothing}
    stiffness::Float64
    damping::Float64
end

struct Pulley
    tethers::Tuple{Int, Int}
    sum_length::Float64
end

# protune-speed-system.jpg
points = [
    Point(:fixed,  [0, 0, 2],  zeros(3), nothing, zeros(3)),  # Fixed point
    Point(:fixed,  [1, 0, 2],  zeros(3), nothing, zeros(3)),  # Fixed point
    Point(:fixed,  [2, 0, 2],   zeros(3), nothing, zeros(3)),  # Fixed point
    Point(:fixed,  [5.5, 0, 2], zeros(3), nothing, zeros(3)),  # Fixed point
   
    Point(:quasi_static, [1, 0, -1], zeros(3), 1.0, zeros(3)),

    Point(:quasi_static, [0.5, 0, -2], zeros(3), 1.0, zeros(3)),
    Point(:quasi_static, [1.5, 0, -2], zeros(3), 1.0, zeros(3)),

    Point(:quasi_static, [1, 0, -3], zeros(3), 1.0, zeros(3)),
    Point(:quasi_static, [2, 0, -3], zeros(3), 1.0, zeros(3)),

    Point(:dynamic,  [1, 0, -10], zeros(3), 1.0, [0., 0., -50]),
    Point(:dynamic, [2, 0, -10], zeros(3), 1.0, [0., 0., -50]),
]

stiffness = 614600
damping = 4730
tethers = [
    Tether((1, 6), norm(points[1].position - points[6].position), stiffness, damping),
    Tether((2, 5), norm(points[2].position - points[5].position), stiffness, damping),
    Tether((3, 7), norm(points[3].position - points[7].position), stiffness, damping),
    Tether((4, 9), norm(points[4].position - points[9].position) - 0.5, stiffness, damping),
    
    Tether((5, 6), norm(points[5].position - points[6].position), stiffness, damping),
    Tether((5, 7), norm(points[5].position - points[7].position), stiffness, damping),
    
    Tether((6, 8), norm(points[6].position - points[8].position), stiffness, damping),
    Tether((7, 8), norm(points[7].position - points[8].position), stiffness, damping),
    Tether((7, 9), norm(points[7].position - points[9].position), stiffness, damping),
    
    Tether((8, 10), norm(points[8].position - points[10].position), stiffness, damping),
    Tether((9, 11), norm(points[9].position - points[11].position), stiffness, damping),
]

pulleys = [
    Pulley((5, 6), (tethers[5].l0 + tethers[6].l0))
    Pulley((8, 9), (tethers[8].l0 + tethers[9].l0))
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
    duration = 5.0                                # duration of the simulation        [s]
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
    POS0, VEL0, L0, V0
end

function calc_spring_forces(pos::AbstractMatrix{T}, pulley_l0) where T
    spring_force = zeros(T, length(tethers))
    spring_force_vec = zeros(T, 3, length(tethers))
    segment = zeros(T, 3)
    unit_vector = zeros(T, 3)
    rel_vel = zeros(T, 3)
    # loop over all tethers to calculate spring forces
    for (tether_idx, tether) in enumerate(tethers)
        found = false
        for (pulley_idx, pulley) in enumerate(pulleys)
            if tether_idx == pulley.tethers[1] # each tether should only be part of one pulley
                l0 = pulley_l0[pulley_idx]
                found = true
                break
            elseif tether_idx == pulley.tethers[2]
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
@register_symbolic calc_spring_forces(pos, pulley_l0)

function calc_acc(se::Settings3, pos::AbstractMatrix{T}, pulley_l0) where T
    spring_force_vec, spring_force = calc_spring_forces(pos, pulley_l0)

    pulley_acc = zeros(T, length(pulleys))
    for (pulley_idx, pulley) in enumerate(pulleys)
        M = 3.1
        pulley_force = spring_force[pulley.tethers[1]] - spring_force[pulley.tethers[2]]
        pulley_acc[pulley_idx] = pulley_force / M
    end

    acc = zeros(T, 3, length(points))
    force = zeros(T, 3)
    for (point_idx, point) in enumerate(points)
        if point.type === :fixed
            acc[:, point_idx] .= 0.0
        else
            force .= 0.0
            for (j, tether) in enumerate(tethers)
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
    return [acc, pulley_acc]
end
@register_symbolic calc_acc(se::Settings3, pos, pulley_l0)

function calc_pos()
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
        D(pulley_vel) ~ pulley_acc - 10pulley_vel
    ]

    eqs = [
        eqs
        vec([acc, pulley_acc]) ~ vec(calc_acc(se, pos, pulley_l0))
    ]
    
    eqs = reduce(vcat, Symbolics.scalarize.(eqs))
    @named sys = ODESystem(eqs, t)
    @time sys = structural_simplify(sys; simplify=false)
    sys, pos, vel
end

function simulate(se, sys)
    dt = 0.1
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    POS0, VEL0, L0, V0 = calc_initial_state(points, tethers, pulleys)
    @time prob = ODEProblem(sys, [pos => POS0, vel => VEL0, sys.pulley_l0 => L0, sys.pulley_vel => V0], tspan; simplify=false)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=true); dt, abstol=tol, reltol=tol, saveat=ts)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=true); dt, abstol=tol, reltol=tol, saveat=ts)
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
        segs = [[tethers[k].points[1], tethers[k].points[2]] for k in eachindex(tethers)]
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
    global sol, pos, vel, len, c_spr
    se = Settings3()
    set_tether_diameter!(se, se.d_tether) # adapt spring and damping constants to tether diameter
    sys, pos, vel = model(se)
    sol, elapsed_time = simulate(se, sys)
    play(se, sol, pos)
    println("Elapsed time: $(elapsed_time) s, speed: $(round(se.duration/elapsed_time)) times real-time")
    println("Number of evaluations per step: ", round(sol.stats.nf/(se.duration/0.1), digits=1))
    println("8 pos: ", sol[pos[:, 8]][end])
    println("9 pos: ", sol[pos[:, 9]][end])
    sol, pos, vel, sys
end

if (! @isdefined __BENCH__) || __BENCH__ == false
    sol, pos, vel, sys = main()
end
__BENCH__ = false
nothing
