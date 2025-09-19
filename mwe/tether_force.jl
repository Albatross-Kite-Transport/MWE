using SymbolicAWEModels
using OrdinaryDiffEqBDF
using LinearAlgebra
using OrdinaryDiffEqCore
using ControlSystems
using Plots
using UnPack
using ModelingToolkit
import ModelingToolkit: t_nounits as t

@variables begin
    wing_pos(t)[1, 1:3]
    wing_vel(t)[1, 1:3]
    ω_b(t)[1, 1:3],
    Q_b_w(t)[1, 1:3],
    tether_vel(t)[1:3],
    tether_len(t)[1:3],
    tether_acc(t)[1:3],
    twist_ω(t)[1:2],
    free_twist_angle(t)[1:2]

    heading(t)[1]
    elevation(t)[1]
    elevation_vel(t)[1]
    azimuth(t)[1]
    azimuth_vel(t)[1]
    winch_force(t)[1:3]
    angle_of_attack(t)[1]
    wind_scale_gnd(t)
end

outputs = [
    heading[1],
    elevation[1],
    azimuth[1],
    tether_len...,
    tether_vel...,
    tether_acc...,
    elevation_vel[1],
    azimuth_vel[1],

    # unmeasured
    winch_force...,
    angle_of_attack[1],
    wind_scale_gnd
]

vy = [
    "heading [rad]",
    "elevation [rad]",
    "azimuth [rad]",
    "main tether length [m]",
    "left tether length [m]",
    "right tether length [m]",
    "main tether vel [m/s]",
    "left tether vel [m/s]",
    "right tether vel [m/s]",
    "main tether acc [m/s^2]",
    "left tether acc [m/s^2]",
    "right tether acc [m/s^2]",
    "elevation vel [rad/s]",
    "azimuth vel [rad/s]",
    "main winch force [N]",
    "left winch force [N]",
    "right winch force [N]",
    "angle of attack [rad]",
    "wind vel [m/s]",
]

@info "Outputs: $outputs"
ny = length(outputs)
nu = 3

"""
  To use the kalman filter:
  f
      - state = measured y
      - step
  h
      - return complete y
"""

# if !@isdefined simple_sam
    set = Settings("system.yaml")
    dt = 1/set.sample_freq
    sam = SymbolicAWEModel(set, "ram")
    plant_set = Settings("system.yaml")
    plant_sam = SymbolicAWEModel(plant_set, "ram")
    tether_set = Settings("system.yaml")
    tether_sam = SymbolicAWEModel(tether_set, "tether")
    simple_set = Settings("system.yaml")
    simple_sam = SymbolicAWEModel(simple_set, "simple_ram")
# end

init!(sam; outputs)
init!(plant_sam; outputs)
init!(tether_sam)
init!(simple_sam; outputs, create_control_func=true)
[group.damping = 100 for group in simple_sam.sys_struct.groups]

find_steady_state!(sam)
find_steady_state!(plant_sam)
copy_to_simple!(sam, tether_sam, simple_sam; prn=false)
d_tether_len = [swinch.tether_len - winch.tether_len for (winch, swinch) in
                zip(sam.sys_struct.winches, simple_sam.sys_struct.winches)]
find_steady_state!(simple_sam)

struct SAMEstim{GetterType}
    sam::SymbolicAWEModel
    simple_sam::SymbolicAWEModel
    tether_sam::SymbolicAWEModel
    get_x::GetterType
    ŷ::Vector{Float64}
    integ::Vector{Float64}
    function SAMEstim(sam, simple_sam, tether_sam)
        get_x = ModelingToolkit.getu(simple_sam.integrator, simple_sam.control_funcs.dvs)
        ŷ = simple_sam.simple_lin_model.get_y(simple_sam.integrator)
        integ = zeros(2)
        new{typeof(get_x)}(sam, simple_sam, tether_sam, get_x, ŷ, integ)
    end
end

function preparestate!(estim::SAMEstim, y::Vector{<:Real})
    @unpack sam, simple_sam, tether_sam, get_x, ŷ, integ = estim

    # update sam with y
    @unpack winches, transforms, wings, points, groups = simple_sam.sys_struct
    u = [winch.set_value for winch in winches]
    transforms[1].heading = y[1] - simple_sam.sys_struct.wings[1].heading
    # TODO: singularity when elevation = 0.0
    transforms[1].elevation = y[2]
    transforms[1].azimuth = y[3]
    SymbolicAWEModels.reposition!(simple_sam.sys_struct.transforms, simple_sam.sys_struct)
    SymbolicAWEModels.reinit!(simple_sam, simple_sam.prob, FBDF(); lin_vsm=false)
    simple_sam.prob.set_set_values(simple_sam.integrator, u)

    ŷ .= simple_sam.simple_lin_model.get_y(simple_sam.integrator)
    ŷ[4:6] .-= d_tether_len

    # 1. calculate winch force from acc+u
    winch_force = SymbolicAWEModels.calc_winch_force(simple_sam.sys_struct, y[7:9], y[10:12], u)
    # 2. adjust wind_vel based on calculated_winch_force / simple_model_winch_force
    force_fracs = winch_force ./ [norm(winch.force) for winch in winches]
    simple_sam.set.v_wind = 0.5*simple_sam.set.v_wind + 0.5*simple_sam.set.v_wind * force_fracs[1]
    # 3. adjust wind_vel based on estimated tether len vs measured tether len

    # v = y[13:14] ⋅ [cos(y[1]), sin(y[1])]
    # v̂ = ŷ[13:14] ⋅ [cos(ŷ[1]), sin(ŷ[1])]
    # wings[1].drag_frac -= 0.1(norm(y[13:14]) - norm(ŷ[13:14]))
    wings[1].drag_frac = 1.2

    integ .+= 0.01 * (y[5:6] .- ŷ[5:6])
    # TODO: find out why negative elevation is fucked

    return get_x(simple_sam.integrator)
end

function updatestate!(estim::SAMEstim, u::Vector{<:Real}, y::Vector{<:Real})
    @unpack simple_sam, sam, ŷ, integ = estim
    set_values = u # .* [1, 1+integ[1], 1+integ[2]]
    next_step!(simple_sam; set_values, vsm_interval=3)
    # next_step!(sam; set_values=u)
    nothing
end

function updatestate!(plant_sam::SymbolicAWEModel, u::Vector{<:Real})
    next_step!(plant_sam; set_values=u)
    nothing
end

estim = SAMEstim(sam, simple_sam, tether_sam)

function man_sim!(N, ry, u; u_disturb = [1,1.2,1.2])
    U_data, Y_data, Ŷ_data, Ry_data = zeros(nu, N), zeros(ny, N), zeros(ny, N), zeros(ny, N)
    T_data = collect(1:N)*dt
    for i = 1:N
        @show i
        plant_sam.set.v_wind = 15.51 + 4sin(2π * i*dt / 5)
        y = plant_sam.simple_lin_model.get_y(plant_sam.integrator) .* (1 .+ 0.01 .* rand(19))
        x̂ = preparestate!(estim, y)
        ŷ = estim.ŷ
        # TODO: fix u integrator
        u = 0.9*u + 0.1*SymbolicAWEModels.calc_steady_torque(simple_sam)
        # TODO: use trajectorylimiter to directly set tether_acc for a set tether_len
        # and calculate du
        U_data[:,i], Y_data[:,i], Ŷ_data[:,i], Ry_data[:,i] = u, y, ŷ, ry
        updatestate!(estim, u, y) # in the estim: step the complex model
        u_plant = u .* [1, 1+estim.integ[1], 1+estim.integ[2]] .* u_disturb
        @show [1, 1+estim.integ[1], 1+estim.integ[2]] .* u_disturb
        updatestate!(plant_sam, u_plant)  # update plant simulator
    end
    return U_data, Y_data, Ŷ_data, Ry_data, T_data
end


y0 = sam.simple_lin_model.get_y(sam.integrator)
x0 = estim.get_x(estim.simple_sam.integrator)
@show y0
ry = copy(y0)
ry[4:6] .-= 0.5
u0 = SymbolicAWEModels.calc_steady_torque(plant_sam)
U_data, Y_data, Ŷ_data, Ry_data, T_data = man_sim!(400, ry, u0)


function plot_state(plot_idxs)
    # vy = string.(outputs)
    # Prepare the y-axis labels for each pair (actual and estimated) per index
    labels = [["Actual $(vy[idx])" "Estimated $(vy[idx])"] for idx in plot_idxs]
    # Plot all in one call with subplots layout, each subplot having the corresponding data and labels
    p = plot(layout = (length(plot_idxs), 1), size=(900,600))
    for (i, idx) in enumerate(plot_idxs)
        plot!(p[i], T_data,
              [Y_data[idx, :], Ŷ_data[idx, :]],
              label = labels[i])
    end
    xlabel!(p[end], "time [s]")  # add xlabel only once at bottom subplot
    display(p)
    return p
end
plot_idxs = [4,5,6,8,9,19]
plot_state(plot_idxs)
# sl_plant, _ = sim_oscillate!(plant_sam)
# sl_simple, _ = sim_oscillate!(simple_sam)
# plot(sl_plant.syslog.time, [sl_plant.syslog.heading, sl_simple.syslog.heading]; label=["plant" "simple"])
