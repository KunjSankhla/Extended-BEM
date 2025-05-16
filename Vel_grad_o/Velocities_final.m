Velocities_domain = zeros(total_domain,3);
% % 
% % fprintf("Processing Velocities_final %d of %d",m,total_domain)
% %     fprintf("\n")


parfor m = 1:total_domain
    fprintf("Processing Velocities_final %d of %d",m,total_domain)
    fprintf("\n")
    xo = dom_int_coordinates(m,1);
    yo = dom_int_coordinates(m,2);
    zo = dom_int_coordinates(m,3);
    for n =  1:N_total
        x = int_coordinates(n,1);
        y = int_coordinates(n,2);
        z = int_coordinates(n,3);
        % boundary 1 -
        if (n >= 1) && (n <= Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = y + y_diff_start;
            lu2 = y - y_diff_start;
            for i = 1:3
                for j = 1:3
                    Velocities_domain(m,i) = Velocities_domain(m,i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) + ...
                        Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                end
            end

            % for boundary 2 -
        elseif(n> Nx*Ny) && (n <= Nx*Nz + Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;
            for i = 1:3
                for j = 1:3
                    Velocities_domain(m,i) = Velocities_domain(m,i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2).*R_mat(n,j) - ...
                        Tijk_calc_int(x,xo,y,yo,z,zo,i,j,2,1,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                end
            end

            % for boundary 3 -
        elseif(n> Nx*Nz + Nx*Ny) && (n <= Nx*Nz + 2*Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = y + y_diff_start;
            lu2 = y - y_diff_start;
            for i = 1:3
                for j = 1:3
                    Velocities_domain(m,i) = Velocities_domain(m,i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) - ...
                        Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                end
            end

            % for boundary 4 -
        elseif(n> Nx*Nz + 2*Nx*Ny) && (n <= 2*Nx*Nz + 2*Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;
            for i = 1:3
                for j = 1:3
                    Velocities_domain(m,i) = Velocities_domain(m,i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2).*R_mat(n,j) + ...
                        Tijk_calc_int(x,xo,y,yo,z,zo,i,j,2,1,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                end
            end

            % for boundary 5 -
        elseif(n> 2*Nx*Nz + 2*Nx*Ny) && (n <= 2*Nx*Nz + 2*Nx*Ny + Ny*Nz)
            ll1 = y - x_diff_start;
            lu1 = y + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;
            for i = 1:3
                for j = 1:3
                    Velocities_domain(m,i) = Velocities_domain(m,i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2).*R_mat(n,j) + ...
                        Tijk_calc_int(x,xo,y,yo,z,zo,i,j,1,2,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                end
            end

            % for boundary 6 -
        elseif(n> 2*Nx*Nz + 2*Nx*Ny + Ny*Nz) && (n <= 2*Nx*Nz + 2*Nx*Ny + 2*Ny*Nz)
            ll1 = y - x_diff_start;
            lu1 = y + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;
            for i = 1:3
                for j = 1:3
                    Velocities_domain(m,i) = Velocities_domain(m,i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2).*R_mat(n,j) - ...
                        Tijk_calc_int(x,xo,y,yo,z,zo,i,j,1,2,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                end
            end
        end
    end
end
Ud_x = Velocities_domain(:,1);
Ud_y = Velocities_domain(:,2);
Ud_z = Velocities_domain(:,3);
domx = dom_int_coordinates(:,1);
domy = dom_int_coordinates(:,2);
domz = dom_int_coordinates(:,3);
%%
quiver3(domx,domz,domy,Ud_x,Ud_z,Ud_y)