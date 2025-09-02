# SPDX-FileCopyrightText: 2025 Bart van de Lint <bart@vandelint.net>
#
# SPDX-License-Identifier: MPL-2.0

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end

using SymbolicAWEModels, ControlPlots

set = Settings("system.yaml")
sam = SymbolicAWEModel(set, "ram")
init!(sam)

tether_set = Settings("system.yaml")
tether_sam = SymbolicAWEModel(tether_set, "tether")
init!(tether_sam)

simple_set = Settings("system.yaml")
simple_sam = SymbolicAWEModel(simple_set, "simple_ram")
init!(simple_sam)

bias = 0.3
find_steady_state!(sam)
SymbolicAWEModels.copy_to_simple!(sam, tether_sam, simple_sam)
simple_sam.sys_struct.wings[1].drag_frac = 1.2

sl, _ = sim_oscillate!(sam; total_time=5.0, prn=true, bias) # TODO: add first frac ram model
display(plot(sam.sys_struct, sl; plot_default=false, plot_elevation=true,
             plot_aoa=true, plot_heading=true, plot_aero_force=true))

find_steady_state!(simple_sam)
simple_sl, _ = sim_oscillate!(simple_sam; total_time=5.0, prn=true, bias)
display(plot(simple_sam.sys_struct, simple_sl; plot_default=false, plot_elevation=true,
             plot_aoa=true, plot_heading=true, plot_aero_force=true))
