using SymbolicAWEModels
using OrdinaryDiffEqCore
using OrdinaryDiffEqBDF
using LinearAlgebra
using ControlSystems
using Plots
using UnPack
using ModelingToolkit
using Optim # Use Optim for the optimizer
import ModelingToolkit: t_nounits as t

# --- Variable, Output, and Label Definitions (No changes) ---
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
    elevation_acc(t)[1]
    azimuth(t)[1]
    azimuth_vel(t)[1]
    azimuth_acc(t)[1]
    winch_force(t)[1:3]
    angle_of_attack(t)[1]
    wind_scale_gnd(t)
end

outputs = [
    elevation[1],
    elevation_vel[1],
    elevation_acc[1],
    azimuth[1],
    azimuth_vel[1],
    azimuth_acc[1],
    tether_len...,
    tether_vel...,
    tether_acc...,
    # unmeasured
    heading[1],
    winch_force...,
    angle_of_attack[1],
    wind_scale_gnd
]
vy = string.(outputs)
ny = length(outputs)
nu = 3

# --- Model Initialization ---
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

# --- Estimator Struct and Update Functions (No changes) ---
struct SAMEstim{GetterType}
    sam::SymbolicAWEModel
    simple_sam::SymbolicAWEModel
    tether_sam::SymbolicAWEModel
    get_x::GetterType
    ŷ::Vector{Float64}
    function SAMEstim(sam, simple_sam, tether_sam)
        get_x = ModelingToolkit.getu(simple_sam.integrator, simple_sam.control_funcs.dvs)
        ŷ = simple_sam.simple_lin_model.get_y(simple_sam.integrator)
        new{typeof(get_x)}(sam, simple_sam, tether_sam, get_x, ŷ)
    end
end

function updatestate!(estim::SAMEstim, u::Vector{<:Real})
    @unpack simple_sam = estim
    next_step!(simple_sam; set_values=u, vsm_interval=3)
    nothing
end

function updatestate!(plant_sam::SymbolicAWEModel, u::Vector{<:Real})
    next_step!(plant_sam; set_values=u)
    nothing
end

estim = SAMEstim(sam, simple_sam, tether_sam)

# --- `preparestate!` FUNCTION WITH OPTIMIZATION LOGIC ---
function preparestate!(estim::SAMEstim, y::Vector{<:Real}, p_initial_guess::Vector{<:Real})
    
    function cost_function(p::Vector{Float64})
        @unpack simple_sam = estim
        wind_speed, wind_dir, drag_frac = p[1], p[2], p[3]
        # @show p
        
        simple_sam.set.v_wind = wind_speed
        simple_sam.set.upwind_dir = wind_dir
        simple_sam.sys_struct.wings[1].drag_frac = drag_frac
        
        OrdinaryDiffEqCore.reinit!(simple_sam.integrator)
        ŷ_predicted = simple_sam.simple_lin_model.get_y(simple_sam.integrator)
        ŷ_predicted[7:9] .-= d_tether_len
        
        error_indices = [3,6,13]
        error = y[error_indices] .- ŷ_predicted[error_indices]
        # @show error
        return dot(error, error)
    end

    solver = BFGS()
    options = Optim.Options(iterations = 1)
    # options = Optim.Options()
    @time result = optimize(cost_function, p_initial_guess, solver, options)
    
    return Optim.minimizer(result)
end


# --- UPDATED SIMULATION LOOP ---
function man_sim!(N, ry, u)
    U_data, Y_data, Ŷ_data, Ry_data = zeros(nu, N), zeros(ny, N), zeros(ny, N), zeros(ny, N)
    T_data = collect(1:N)*dt
    
    p_estimate = [simple_sam.set.v_wind, 0.0, simple_sam.sys_struct.wings[1].drag_frac]

    for i = 1:N
        @info "Timestep: $i / $N"
        
        plant_sam.set.v_wind = (15.51 + sin(2π * i*dt / 5))
        plant_sam.set.upwind_dir = sin(0.1)
        
        # 1. Get measurements from the "real" system (plant)
        y = plant_sam.simple_lin_model.get_y(plant_sam.integrator)

        # 2. Estimate unmeasured states
        p_estimate = preparestate!(estim, y, p_estimate)

        # 3. Permanently update the simple_sam model with the new best-fit parameters
        estim.simple_sam.set.v_wind = p_estimate[1]
        estim.simple_sam.set.upwind_dir = p_estimate[2]
        estim.simple_sam.sys_struct.wings[1].drag_frac = p_estimate[3]

        # 4. Re-initialize the solver with the new parameters
        # This ensures the changes take effect in the next simulation step.
        OrdinaryDiffEqCore.reinit!(estim.simple_sam.integrator)

        # 5. Recalculate final ŷ for logging and control
        ŷ = estim.simple_sam.simple_lin_model.get_y(estim.simple_sam.integrator)
        ŷ[4:6] .-= d_tether_len

        # 6. Calculate control input
        u = 0.9*u + 0.1*SymbolicAWEModels.calc_steady_torque(simple_sam)
        
        # 7. Store data and step the simulators forward
        U_data[:,i], Y_data[:,i], Ŷ_data[:,i], Ry_data[:,i] = u, y, ŷ, ry
        updatestate!(estim, u)
        updatestate!(plant_sam, u)
    end
    return U_data, Y_data, Ŷ_data, Ry_data, T_data
end

# --- Run Simulation and Plot (no changes) ---
y0 = sam.simple_lin_model.get_y(sam.integrator)
ry = copy(y0)
ry[7:9] .-= 0.5
u0 = SymbolicAWEModels.calc_steady_torque(plant_sam)
U_data, Y_data, Ŷ_data, Ry_data, T_data = man_sim!(200, ry, u0)

function plot_state(plot_idxs)
    labels = [["Actual $(vy[idx])" "Estimated $(vy[idx])"] for idx in plot_idxs]
    p = plot(layout = (length(plot_idxs), 1), size=(900, 600), legend=:outertopright)
    for (i, idx) in enumerate(plot_idxs)
        plot!(p[i], T_data, [Y_data[idx, :], Ŷ_data[idx, :]], label=labels[i])
    end
    xlabel!(p[end], "time [s]")
    display(p)
    return p
end

plot_idxs = [4, 5, 6, 19] # Plot tether lengths and estimated wind speed
plot_state(plot_idxs)
