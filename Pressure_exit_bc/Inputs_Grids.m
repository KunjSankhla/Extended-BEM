Nx = 8;
Ny = 8;
Nz = 8;
Lx = 1;
Ly = 1;
Lz = 1;
mu = 1;
Nxi = 1;
Nyi = 1;
Nzi = Nz;
P_exit = 1;
h = 0.05*Lz;
% forming the coordinate options for each side;
N_total = 2*Nx*Ny + 2*Ny*Nz + 2*Nx*Nz;

x_diff_start = Lx/(2*Nx);
x_diff = Lx/Nx;
y_diff_start = Ly/(2*Ny);
y_diff = Ly/Ny;
z_diff_start =Lz/(2*Nz);
z_diff = Lz/Nz;
bound_coordinates = zeros(N_total,3);
p = P_exit;
% surface 1 -
m = 0;
for i = 1:Ny
    for j = 1:Nx
        bound_coordinates(j+m, 1) = x_diff_start + x_diff*(j-1);
        bound_coordinates(j+m,2) = y_diff_start + y_diff*(i-1);
        bound_coordinates(j+m,3) = 0;
    end
    m = m + Nx;
end

% surface 2 - 
m = Nx*Ny;
for i = 1:Nz
    for j = 1:Nx
        bound_coordinates(j+m, 1) = x_diff_start + x_diff*(j-1);
        bound_coordinates(j+m,2) = Ly;
        bound_coordinates(j+m,3) = z_diff_start + z_diff*(i-1);
    end
    m = m + Nx;
end

% surface 3 - 
m = Nx*Ny + Nz*Nx;
for i = 1:Ny
    for j = 1:Nx
        bound_coordinates(j+m, 1) = x_diff_start + x_diff*(j-1);
        bound_coordinates(j+m,2) = Ly-(y_diff_start + y_diff*(i-1));
        bound_coordinates(j+m,3) = Lz;
    end
    m = m + Nx;
end

% surface 4 - 
m = 2*Nx*Ny + Nx*Nz;
for i = 1:Nz
    for j = 1:Nx
        bound_coordinates(j+m, 1) = x_diff_start + x_diff*(j-1);
        bound_coordinates(j+m,2) = 0;
        bound_coordinates(j+m,3) = Lz-(z_diff_start + z_diff*(i-1));
    end
    m = m + Nx;
end

% surface 5 - 
m = 2*Nx*Ny + 2*Nx*Nz;
for i = 1:Ny
    for j = 1:Nz
        bound_coordinates(j+m, 1) = 0;
        bound_coordinates(j+m,2) = y_diff_start + y_diff*(i-1);
        bound_coordinates(j+m,3) = Lz-(z_diff_start + z_diff*(j-1));
    end
    m = m + Nz;
end

% surface 6 - 
m = 2*Nx*Ny + 2*Nx*Nz + Ny*Nz;
for i = 1:Ny
    for j = 1:Nz
        bound_coordinates(j+m, 1) = Lx;
        bound_coordinates(j+m,2) = y_diff_start + y_diff*(i-1);
        bound_coordinates(j+m,3) = (z_diff_start + z_diff*(j-1));
    end
    m = m + Nz;
end

Velocity_boundaries = zeros(N_total,3);

for n = 1:N_total
% boundary 1 - 
if (n >= 1) && (n <= Nx*Ny)
    Velocity_boundaries(n,1) = 0; % x-velocity;
    Velocity_boundaries(n,2) = 0; % y-velocity;
    Velocity_boundaries(n,3) = 1; % z-velocity;

% for boundary 2 -
elseif(n> Nx*Ny) && (n <= Nx*Nz + Nx*Ny)
    Velocity_boundaries(n,1) = 0; % x-velocity;
    Velocity_boundaries(n,2) = 0; % y-velocity;
    Velocity_boundaries(n,3) = 0; % z-velocity;

% for boundary 3 -
elseif(n> Nx*Nz + Nx*Ny) && (n <= Nx*Nz + 2*Nx*Ny)
    Velocity_boundaries(n,1) = 0; % x-velocity;
    Velocity_boundaries(n,2) = 0; % y-velocity;
    Velocity_boundaries(n,3) = 0; % z-velocity;

% for boundary 4 -
elseif(n> Nx*Nz + 2*Nx*Ny) && (n <= 2*Nx*Nz + 2*Nx*Ny)
    Velocity_boundaries(n,1) = 0; % x-velocity;
    Velocity_boundaries(n,2) = 0; % y-velocity;
    Velocity_boundaries(n,3) = 0; % z-velocity;
% for boundary 5 -
elseif(n> 2*Nx*Nz + 2*Nx*Ny) && (n <= 2*Nx*Nz + 2*Nx*Ny + Ny*Nz)
    Velocity_boundaries(n,1) = 0; % x-velocity;
    Velocity_boundaries(n,2) = 0; % y-velocity;
    Velocity_boundaries(n,3) = 0; % z-velocity;

% for boundary 6 -
elseif(n> 2*Nx*Nz + 2*Nx*Ny + Ny*Nz) && (n <= 2*Nx*Nz + 2*Nx*Ny + 2*Ny*Nz)
    Velocity_boundaries(n,1) = 0; % x-velocity;
    Velocity_boundaries(n,2) = 0; % y-velocity;
    Velocity_boundaries(n,3) = 0; % z-velocity;
end
end

%for inlet flow from the mid part only
% for n = 1:Nx*Ny
%     t = 0.25*Ny*Nx;
%     if (n > t + 0.25*Nx) && ( n <= t + 0.75*Nx) || (n > t + Nx + 0.25*Nx) && ( n <= t + Nx+  0.75*Nx) ...
%             || (n > t + 2*Nx + 0.25*Nx) && ( n <= t + 2*Nx+  0.75*Nx) || (n > t + 3*Nx + 0.25*Nx) && ( n <= t + 3*Nx+  0.75*Nx) ...
%             || (n > t + 4*Nx + 0.25*Nx) && ( n <= t + 4*Nx+  0.75*Nx) || (n > t + 5*Nx + 0.25*Nx) && ( n <= t + 5*Nx+  0.75*Nx) 
%         Velocity_boundaries(n,3) = 1; % z-velocity;
%     end
% end


%internal/domain numbering

total_domain = Nxi*Nyi*Nzi;
x_diff_starti = Lx/(2*Nxi);
x_diffi = Lx/Nxi;
y_diff_starti = Ly/(2*Nyi);
y_diffi = Ly/Nyi;
z_diff_starti =Lz/(2*Nzi);
z_diffi = Lz/Nzi;
m = 1;
dom_int_coordinates = zeros(total_domain,3);
x_internal  = x_diff_starti:x_diffi:Lx;
y_internal  = y_diff_starti:y_diffi:Ly;
z_internal  = z_diff_starti:z_diffi:Lz;

for a = 1:Nzi
    for b = 1:Nyi
        for c= 1:Nxi
            dom_int_coordinates(m,1) = x_internal(c);
            dom_int_coordinates(m,2) = y_internal(b);
            dom_int_coordinates(m,3) = z_internal(a);
            m = m+1;
        end
    end
end


% create an internal boundary for the 3rd(exit boundary)
N_total_int = N_total + Nx*Ny;
int_coordinates = zeros(N_total_int,3);
for i = 1:N_total_int
    if i <= N_total
        int_coordinates(i,1) = bound_coordinates(i,1);
        int_coordinates(i,2) = bound_coordinates(i,2);
        int_coordinates(i,3) = bound_coordinates(i,3);
    elseif i > N_total
        int_coordinates(i,1) = bound_coordinates(i-(2*Ny*Nz+Nx*Nz+Nx*Ny),1);
        int_coordinates(i,2) = bound_coordinates(i-(2*Ny*Nz+Nx*Nz+Nx*Ny),2);
        int_coordinates(i,3) = bound_coordinates(i-(2*Ny*Nz+Nx*Nz+Nx*Ny),3) - h;
    end
end

% quiver3(bound_coordinates(:,1),bound_coordinates(:,2),bound_coordinates(:,3),Velocity_boundaries(:,1),Velocity_boundaries(:,2),Velocity_boundaries(:,3))
% check the input velocities