# MWE of tether system creation and simplification
using ModelingToolkit, LinearAlgebra
using ModelingToolkit: t_nounits as t, D_nounits as D

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

function main()
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

    g_earth::Vector{Float64} = [0.0, 0.0, -9.81] # gravitational acceleration     [m/s²]
    l0 = 50                                      # initial tether length             [m]
    c_spring = 614600                            # unit spring constant              [N]
    rel_compression_stiffness = 0.01             # relative compression stiffness    [-]
    damping = 473                                # unit damping constant            [Ns]

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

    function calc_spring_forces(eqs, pos, vel, pulley_l0)
        # loop over all segments to calculate spring forces
        @variables begin
            spring_force(t)[eachindex(segments)]
            spring_force_vec(t)[1:3, eachindex(segments)]
            segment(t)[1:3, eachindex(segments)]
            unit_vector(t)[1:3, eachindex(segments)]
            rel_vel(t)[1:3, eachindex(segments)]
            len(t)[eachindex(segments)]
            spring_vel(t)[eachindex(segments)]
            l0(t)[eachindex(segments)]
        end
        for (segment_idx, tether) in enumerate(segments)
            found = false
            for (pulley_idx, pulley) in enumerate(pulleys)
                if segment_idx == pulley.segments[1] # each tether should only be part of one pulley
                    eqs = [
                        eqs
                        l0[segment_idx] ~ pulley_l0[pulley_idx]
                    ]
                    found = true
                    break
                elseif segment_idx == pulley.segments[2]
                    eqs = [
                        eqs
                        l0[segment_idx] ~ pulley.sum_length - pulley_l0[pulley_idx]
                    ]
                    found = true
                    break
                end
            end
            if !found
                eqs = [
                    eqs
                    l0[segment_idx] ~ tether.l0
                ]
            end
            p1, p2 = tether.points[1], tether.points[2]

            eqs = [
                eqs
                segment[:, segment_idx] ~ pos[:, p2] - pos[:, p1]
                len[segment_idx]        ~ norm(segment[:, segment_idx])
                unit_vector[:, segment_idx] ~ segment[:, segment_idx] / len[segment_idx]
                rel_vel[:, segment_idx]     ~ vel[:, p1] .- vel[:, p2]
                spring_vel[segment_idx]     ~ rel_vel[:, segment_idx] ⋅ unit_vector[:, segment_idx]
                spring_force[segment_idx]   ~ (tether.stiffness * tether.l0 * (len[segment_idx] - l0[segment_idx]) - tether.damping * tether.l0 * spring_vel[segment_idx])
                spring_force_vec[:, segment_idx]  ~ spring_force[segment_idx] * unit_vector[:, segment_idx]
            ]
        end
        return eqs, spring_force_vec, spring_force
    end

    function calc_acc(eqs, pos, vel, pulley_l0)
        eqs, spring_force_vec, spring_force = calc_spring_forces(eqs, pos, vel, pulley_l0)
        
        @variables pulley_acc(t)[eachindex(pulleys)]
        @variables pulley_force(t)[eachindex(pulleys)]
        for (pulley_idx, pulley) in enumerate(pulleys)
            M = 3.1
            eqs = [
                eqs
                pulley_force[pulley_idx] ~ spring_force[pulley.segments[1]] - spring_force[pulley.segments[2]]
                pulley_acc[pulley_idx] ~ pulley_force[pulley_force] / M
            ]
        end

        @variables acc(t)[1:3, eachindex(points)]
        @variables force(t)[1:3, eachindex(points)]
        for (point_idx, point) in enumerate(points)
            if point.type === :fixed
                eqs = [
                    eqs
                    acc[:, point_idx] ~ zeros(3)
                ]
            else
                F = zeros(Num, 3)
                for (j, tether) in enumerate(segments)
                    if point_idx in tether.points
                        inverted = tether.points[2] == point_idx
                        if inverted
                            F .-= spring_force_vec[:, j]
                        else
                            F .+= spring_force_vec[:, j]
                        end
                    end
                end
                eqs = [
                    eqs
                    force[:, point_idx] ~ F
                    acc[:, point_idx] ~ force[:, point_idx] / point.mass + g_earth
                ]
            end
        end
        return eqs, acc, pulley_acc
    end

    function model()
        @parameters rel_compression_stiffness = rel_compression_stiffness
        @variables begin
            pos(t)[1:3, eachindex(points)]
            vel(t)[1:3, eachindex(points)]
            acc(t)[1:3, eachindex(points)]
            force(t)[1:3, eachindex(points)]

            pulley_force(t)[eachindex(pulleys)]
            pulley_l0(t)[eachindex(pulleys)] # first tether length in pulley
            pulley_vel(t)[eachindex(pulleys)]
            pulley_acc(t)[eachindex(pulleys)]
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

        eqs, acc, pulley_acc = calc_acc(eqs, pos, vel, pulley_l0)
        
        eqs = Symbolics.scalarize.(reduce(vcat, Symbolics.scalarize.(eqs)))

        @named sys = ODESystem(eqs, t)     
        @time sys = structural_simplify(sys)
        
        sys, pos, vel, defaults, guesses
    end
    model()
end

sys, pos, vel, defaults, guesses = main()


nothing
