using SymbolicAWEModels
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
    twist_ω(t)[1:2],
    free_twist_angle(t)[1:2]

    heading(t)[1]
    elevation(t)[1]
    elevation_vel(t)[1]
    elevation_acc(t)[1]
    azimuth(t)[1]
    azimuth_vel(t)[1]
    azimuth_acc(t)[1]
    winch_force(t)[1:3]
    angle_of_attack(t)[1]
end

outputs = [
    heading[1],
    elevation[1],
    azimuth[1],
    tether_len...,
    tether_vel...,
    winch_force...,
    elevation_vel[1],
    azimuth_vel[1],
    elevation_acc[1],
    azimuth_acc[1],

    angle_of_attack[1],
]
@info "Outputs: $outputs"
ny = length(outputs)
nu = 3

# if !@isdefined simple_sam
    set = Settings("system.yaml")
    dt = 1/set.sample_freq
    sam = SymbolicAWEModel(set, "ram")
    plant_set = Settings("system.yaml")
    plant_set.v_wind *= 0.9
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

find_steady_state!(sam)
find_steady_state!(plant_sam)
copy_to_simple!(sam, tether_sam, simple_sam; prn=false)
d_tether_len = [swinch.tether_len - winch.tether_len for (winch, swinch) in
                zip(sam.sys_struct.winches, simple_sam.sys_struct.winches)]
simple_sam.sys_struct.wings[1].drag_frac = 1.0
find_steady_state!(simple_sam)

struct SAMEstim{GetterType}
    sam::SymbolicAWEModel
    simple_sam::SymbolicAWEModel
    tether_sam::SymbolicAWEModel
    get_x::GetterType
    function SAMEstim(sam, simple_sam, tether_sam)
        get_x = ModelingToolkit.getu(simple_sam.integrator, simple_sam.control_funcs.dvs)
        new{typeof(get_x)}(sam, simple_sam, tether_sam, get_x)
    end
end

function preparestate!(estim::SAMEstim, y::Vector{<:Real})
    @unpack sam, simple_sam, tether_sam, get_x = estim

    # update sam with y
    simple_sam.sys_struct.transforms[1].heading = y[1] - simple_sam.sys_struct.wings[1].heading
    simple_sam.sys_struct.transforms[1].elevation = y[2]
    simple_sam.sys_struct.transforms[1].azimuth = y[3]
    SymbolicAWEModels.reposition!(simple_sam.sys_struct.transforms, simple_sam.sys_struct)
    OrdinaryDiffEqCore.reinit!(simple_sam.integrator; reinit_dae=true)

    ŷ = simple_sam.simple_lin_model.get_y(simple_sam.integrator)
    ŷ[4:6] .-= d_tether_len

    simple_sam.set.v_wind += (y[7] - ŷ[7])
    @show simple_sam.set.v_wind
    simple_sam.sys_struct.wings[1].drag_frac += (ŷ[13] - y[13])

    return get_x(simple_sam.integrator), ŷ
end

function updatestate!(estim::SAMEstim, u::Vector{<:Real}, y::Vector{<:Real})
    @unpack simple_sam, sam = estim
    next_step!(simple_sam; set_values=u)
    next_step!(sam; set_values=u)
    nothing
end

function updatestate!(plant_sam::SymbolicAWEModel, u::Vector{<:Real})
    next_step!(plant_sam; set_values=u)
    nothing
end

estim = SAMEstim(sam, simple_sam, tether_sam)

function man_sim!(N, ry, u)
    U_data, Y_data, Ŷ_data, Ry_data = zeros(nu, N), zeros(ny, N), zeros(ny, N), zeros(ny, N)
    T_data = collect(1:N)*dt
    for i = 1:N
        @show i
        y = plant_sam.simple_lin_model.get_y(plant_sam.integrator)
        @time x̂, ŷ = preparestate!(estim, y)
        u = 0.9*u + 0.1*SymbolicAWEModels.calc_steady_torque(simple_sam)
        U_data[:,i], Y_data[:,i], Ŷ_data[:,i], Ry_data[:,i] = u, y, ŷ, ry
        @time updatestate!(estim, u, y) # in the estim: step the complex model
        updatestate!(plant_sam, u)  # update plant simulator
    end
    return U_data, Y_data, Ŷ_data, Ry_data, T_data
end


y0 = sam.simple_lin_model.get_y(sam.integrator)
x0 = estim.get_x(estim.simple_sam.integrator)
@show y0
ry = copy(y0)
u0 = SymbolicAWEModels.calc_steady_torque(plant_sam)
U_data, Y_data, Ŷ_data, Ry_data, T_data = man_sim!(100, ry, u0)


function plot_state(plot_idxs)
    vy = string.(outputs)

    # Prepare the y-axis labels for each pair (actual and estimated) per index
    labels = [["$(vy[idx]) actual" "$(vy[idx]) estimated"] for idx in plot_idxs]

    # Plot all in one call with subplots layout, each subplot having the corresponding data and labels
    p = plot(layout = (length(plot_idxs), 1))

    for (i, idx) in enumerate(plot_idxs)
        plot!(p[i], T_data,
              [Y_data[idx, :], Ŷ_data[idx, :]],
              label = labels[i])
    end

    display(p)
end
plot_idxs = [4, 5, 6, 7, 8, 9]
plot_state(plot_idxs)
