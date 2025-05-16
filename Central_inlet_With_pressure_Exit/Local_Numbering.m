function D = Local_Numbering(N)
% This function corresponds the global numbering to a local numbering
% system with respect to each surfacec only.
run("Inputs_Grids.m")
% checking for which surface it is using the convention used in local
% numbering
if (N >= 1) && (N <= Nx*Ny)
    surface = 1;
elseif (N > Nx*Ny) && (N <= Nx*Ny + Nx*Nz)
    surface = 2;
elseif (N > Nx*Ny + Nx*Nz) && (N <= 2*Nx*Ny + Nx*Nz)
    surface = 3;
elseif (N > 2*Nx*Ny + Nx*Nz) && (N <= 2*Nx*Ny + 2*Nx*Nz)
    surface = 4;
elseif (N > 2*Nx*Ny + 2*Nx*Nz) && (N <= 2*Nx*Ny + 2*Nx*Nz + Ny*Nz)
    surface = 5;
elseif (N >2*Nx*Ny + 2*Nx*Nz + Ny*Nz) && (N <= 2*Nx*Ny + 2*Nx*Nz + 2*Ny*Nz)
    surface = 6;
else
    surface= 7;
end


% calculating the local number corresponding to the entered number.
if (surface == 1)
    D = N;
elseif (surface == 2)
    N_s1 = Nx*Ny;
    D = N - N_s1;
elseif (surface == 3)
    N_s2_s1 = Nx*Ny + Nx*Nz;
    D = N-N_s2_s1;
elseif (surface == 4)
    N_s3_s2_s1 = 2*Nx*Ny + Nx*Nz;
    D = N-N_s3_s2_s1;
elseif(surface == 5)
    N_s4_s3_s2_s1 = 2*Nx*Ny + 2*Nx*Nz;
    D = N-N_s4_s3_s2_s1;
elseif (surface == 6)
    N_s5_s4_s3_s2_s1 = 2*Nx*Ny + 2*Nx*Nz + Ny*Nz;
    D = N - N_s5_s4_s3_s2_s1;
elseif(surface == 7)
    N_s6_s5_s4_s3_s2_s1 = 2*Nx*Ny + 2*Nx*Nz + 2*Ny*Nz;
    D = N - N_s6_s5_s4_s3_s2_s1;

end

end