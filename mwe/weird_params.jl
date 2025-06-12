using ModelingToolkit
nx = 3
ny = 10
jac_ = zeros(1, nx, ny)
@parameters begin
    vsm_jac[1, 1:nx, 1:ny] = jac_
end
vsm_jac[1,:,:]
@parameters vsm_jac[1, 1:nx, 1:ny]
vsm_jac[1,:,:]
