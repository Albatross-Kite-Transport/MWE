using LinearAlgebra

# (elevation, azimuth)
sphere_pos            = deg2rad.([89.0 89.0; 1.0 -1.0])
sphere_vel            = [3 -2; 1 1]
sphere_acc            = zeros(2, 2)

r = (sphere_pos[:, 2] - sphere_pos[:, 1]) / 2
@show r
@show perp_r = [-r[2], r[1]]
@show (sphere_vel[:, 1] - sphere_vel[:, 2])
@show rot_vel = (sphere_vel[:, 1] - sphere_vel[:, 2]) ⋅ (perp_r / norm(r))
@show rot_vel / norm(r)
# sphere_vel × normalize(sphere_pos[:, 2] - sphere_pos[:, 1])