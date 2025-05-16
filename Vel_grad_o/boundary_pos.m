%contains the scripting for dividing boundarys, copy and use :)
% n is the indexing variable 
n = input("number to find position = ");
% boundary 1 - 
if (n >= 1) && (n <= Nx*Ny)
    fprintf("1")

% for boundary 2 -
elseif(n> Nx*Ny) && (n <= Nx*Nz + Nx*Ny)
    fprintf("2")

% for boundary 3 -
elseif(n> Nx*Nz + Nx*Ny) && (n <= Nx*Nz + 2*Nx*Ny)
    fprintf("3")

% for boundary 4 -
elseif(n> Nx*Nz + 2*Nx*Ny) && (n <= 2*Nx*Nz + 2*Nx*Ny)
    fprintf("4")

% for boundary 5 -
elseif(n> 2*Nx*Nz + 2*Nx*Ny) && (n <= 2*Nx*Nz + 2*Nx*Ny + Ny*Nz)
    fprintf("5")

% for boundary 6 -
elseif(n> 2*Nx*Nz + 2*Nx*Ny + Ny*Nz) && (n <= 2*Nx*Nz + 2*Nx*Ny + 2*Ny*Nz)
    fprintf("6")
end