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

    Point(:dynamic, [0.5, 0, -2], zeros(3), 1.0, zeros(3)),
    Point(:dynamic, [1.5, 0, -2], zeros(3), 1.0, zeros(3)),

    Point(:dynamic, [1, 0, -3], zeros(3), 1.0, zeros(3)),
    Point(:dynamic, [2, 0, -3], zeros(3), 1.0, zeros(3)),

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
    
    Tether((5, 6), norm(points[5].position - points[6].position) + 0.01, stiffness, damping),
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
    duration = 0.2                                # duration of the simulation        [s]
    save::Bool = false                           # save png files in folder video
end

@with_kw struct Buffer{T}
    spring_force::Vector{T} = zeros(T, length(tethers))
    spring_force_vec::Matrix{T} = zeros(T, 3, length(tethers))
    segment::Vector{T} = zeros(T, 3)
    unit_vector::Vector{T} = zeros(T, 3)
    rel_vel::Vector{T} = zeros(T, 3)
    pulley_acc::Vector{T} = zeros(T, length(pulleys))
    acc::Matrix{T} = zeros(T, 3, length(points))
    force::Vector{T} = zeros(T, 3)
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

function calc_spring_forces!(b::Buffer, pos::AbstractMatrix{T}, vel::AbstractMatrix{T}, pulley_l0) where T
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

        b.segment         .= pos[:, p2] .- pos[:, p1]
        len                = norm(b.segment)
        b.unit_vector     .= b.segment ./ len
        b.rel_vel         .= vel[:, p1] .- vel[:, p2]
        spring_vel         = b.rel_vel ⋅ b.unit_vector
        b.spring_force[tether_idx]          = (tether.stiffness * tether.l0 * (len - l0) - tether.damping * tether.l0 * spring_vel)
        b.spring_force_vec[:, tether_idx]  .= b.spring_force[tether_idx] .* b.unit_vector
    end
    return b.spring_force_vec, b.spring_force
end

function calc_acc!(b::Buffer, se::Settings3, pos::AbstractMatrix{T}, vel::AbstractMatrix{T}, pulley_l0) where T
    calc_spring_forces!(b, pos, vel, pulley_l0)

    for (pulley_idx, pulley) in enumerate(pulleys)
        M = 3.1
        @time pulley_force = b.spring_force[pulley.tethers[1]] - b.spring_force[pulley.tethers[2]]
        b.pulley_acc[pulley_idx] = pulley_force / M
    end

    for (point_idx, point) in enumerate(points)
        if point.type === :fixed
            b.acc[:, point_idx] .= 0.0
        else
            b.force .= 0.0
            for (j, tether) in enumerate(tethers)
                if point_idx in tether.points
                    inverted = tether.points[2] == point_idx
                    if inverted
                        b.force .-= b.spring_force_vec[:, j]
                    else
                        b.force .+= b.spring_force_vec[:, j]
                    end
                end
            end
            b.force .+= point.force
            b.acc[:, point_idx] .= b.force ./ point.mass .+ se.g_earth
        end
    end
    return b.acc, b.pulley_acc
end

function calc_pos(b::Buffer, se::Settings3, idxs, pos_::AbstractMatrix{T}, vel_, pulley_l0) where T
    pos = copy(pos_)
    vel = copy(vel_)
    function f(u, p)
        pos[:, idxs] .= u
        calc_acc!(b, se, pos, vel, pulley_l0)
        return b.acc[:, idxs]
    end
    u0 = pos[:, idxs]
    prob = NonlinearProblem(f, u0, nothing)
    @time sol = solve(prob, NewtonRaphson(autodiff=AutoFiniteDiff(absstep=1e-8, relstep=1e-8)); abstol=1e-5, reltol=1e-5)
    pos[:, idxs] .= sol.u
    return pos
end
@register_array_symbolic calc_pos(
        b::Buffer, 
        se::Settings3, 
        idxs::AbstractVector{Int}, 
        pos::AbstractMatrix, 
        vel::AbstractMatrix, 
        pulley_l0::AbstractVector) begin
    size = size(pos)
    eltype = Float64
end

function calc_pos()
    b = Buffer{Num}()
    se = Settings3()
    s_idxs = [5]
    d_idxs = setdiff(eachindex(points), s_idxs)
    @show d_idxs

    POS0, VEL0, L0, V0 = calc_initial_state(points, tethers, pulleys)
    @parameters begin
        dynamic_pos[1:3, eachindex(d_idxs)] = POS0[:, d_idxs]
        dynamic_acc[1:3, eachindex(d_idxs)] = VEL0[:, d_idxs]
        vel[1:3, eachindex(points)] = VEL0
        pulley_l0[eachindex(pulleys)] = L0
    end
    @variables begin
        static_pos(t)[1:3, eachindex(s_idxs)] = POS0[:, s_idxs]
        static_acc(t)[1:3, eachindex(s_idxs)]
    end
    eqs = []
    for (j, i) in enumerate(s_idxs)
        eqs = [
            eqs
            acc[:, j] ~ zeros(3)
            acc[:, j] ~ calc_acc!(b, se, pos, vel, pulley_l0)[1][:, i] # WRONG POS HERE
        ]
    end
    @mtkbuild ns = NonlinearSystem(eqs)
    prob = NonlinearProblem(ns)
    @time sol = solve(prob, NewtonRaphson())
    @time sol = solve(prob, NewtonRaphson())
end
calc_pos()
@assert false

function model(b::Buffer, se::Settings3)
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

    defaults = Pair{Num, Real}[]
    eqs = [
        # vec(D(pos)) ~ vec(vel)
        # vec(D(vel)) ~ vec(acc)
        D(pulley_l0) ~ pulley_vel
        D(pulley_vel) ~ pulley_acc - 10pulley_vel
    ]

    quasi_idxs = Int[]
    for (point_idx, point) in enumerate(points)
        if point.type === :fixed
            eqs = [
                eqs
                pos[:, point_idx] ~ point.position
                vel[:, point_idx] ~ zeros(3)
                # D(pos[:, point_idx]) ~ vel[:, point_idx]
                # D(vel[:, point_idx]) ~ acc[:, point_idx]
            ]
            # defaults = [
            #     defaults
            #     [pos[j, point_idx] => POS0[j, point_idx] for j in 1:3]
            #     [vel[j, point_idx] => 0 for j in 1:3]
            # ]
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
            push!(quasi_idxs, point_idx)
            defaults = [
                defaults
                [pos[j, point_idx] => POS0[j, point_idx] for j in 1:3]
                [vel[j, point_idx] => 0 for j in 1:3]
            ]
        else
            println("wrong type")
        end
    end
    for i in quasi_idxs
        eqs = [
            eqs
            pos[:, i] ~ calc_pos(b, se, quasi_idxs, pos, vel, pulley_l0)[:, i]
            vel[:, i] ~ zeros(3)
        ]
    end

    defaults = [
        defaults
        pulley_l0 => L0
        pulley_vel => V0
    ]

    b_num = Buffer{Num}()
    eqs = [
        eqs
        vec(acc) ~ vec(calc_acc!(b_num, se, pos, vel, pulley_l0)[1])
        vec(pulley_acc) ~ vec(calc_acc!(b_num, se, pos, vel, pulley_l0)[2])
    ]
    
    eqs = reduce(vcat, Symbolics.scalarize.(eqs))
    @named sys = ODESystem(eqs, t)
    @time sys = structural_simplify(sys; simplify=false)
    sys, pos, vel, defaults
end

function simulate(se, sys, defaults)
    dt = 0.1
    tol = 1e-6
    tspan = (0.0, se.duration)
    ts    = 0:dt:se.duration
    @time prob = ODEProblem(sys, defaults, tspan; simplify=false)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=AutoFiniteDiff(absstep=1e-5, relstep=1e-5)); dt, abstol=tol, reltol=tol, saveat=ts)
    elapsed_time = @elapsed sol = solve(prob, FBDF(autodiff=AutoFiniteDiff(absstep=1e-5, relstep=1e-5)); dt, abstol=tol, reltol=tol, saveat=ts)
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
    b = Buffer{Float64}()
    set_tether_diameter!(se, se.d_tether) # adapt spring and damping constants to tether diameter
    sys, pos, vel, defaults = model(b, se)
    sol, elapsed_time = simulate(se, sys, defaults)
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
