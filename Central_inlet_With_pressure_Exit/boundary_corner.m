function [bndry,norm] = boundary_corner(lpb3)
% defines positions of points on the boundary 3
run("Inputs_Grids.m")
if (lpb3 == 1) % top left corner
    bndry = 2;
    norm = 1;
elseif (lpb3 == Nx) % top right corner
    bndry = 2;
    norm = 2;
elseif (lpb3 == 1 + Nx*Ny - Nx) % bottom left corner
    bndry = 2;
    norm = 4;
elseif (lpb3 == Nx*Ny) % bottom right corner
    bndry = 2;
    norm = 3;
elseif (lpb3 > 1) && (lpb3 < Nx) % top boundary
    bndry = 1;
    norm = 2;
elseif (lpb3 > Nx ) && (lpb3 < Nx*Ny) && (rem(lpb3,Nx) == 0)% right boundary
    bndry = 1;
    norm = 3;
elseif (lpb3 > 1) && (lpb3 < Nx*Ny - Nx + 1) && (rem(lpb3,Nx) == 1) % left boundary
    bndry = 1;
    norm = 1;
elseif (lpb3 > Nx*Ny - Nx + 1) && (lpb3 < Nx*Ny)% bottom boundary
    bndry = 1;
    norm = 4;
else
    bndry = 0;
    norm = 0;
end