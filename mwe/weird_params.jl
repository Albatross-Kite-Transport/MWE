using ModelingToolkit
nx = 3
ny = 10

@parameters two_dims[1:nx, 1:ny] = zeros(nx,ny)
@show two_dims[1,:]

@parameters three_dims[1, 1:nx, 1:ny]
@show three_dims[1,1,1] # works
@show three_dims[1,1,:] # works
@show three_dims[1,:,:] # works

@parameters three_dims[1, 1:nx, 1:ny] = zeros(1,nx,ny)
@show three_dims[1,1,1] # works
@show three_dims[1,1,:] # doesn't work
@show three_dims[1,:,:] # doesn't work

