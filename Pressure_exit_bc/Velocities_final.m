Velocities_domain_1 = zeros(total_domain,1);
Velocities_domain_2 = zeros(total_domain,1);
Velocities_domain_3 = zeros(total_domain,1);

parfor m = 1:total_domain
    fprintf("Processing Velocities_final %d of %d",m,total_domain)
    fprintf("\n")

    xo = dom_int_coordinates(m,1);
    yo = dom_int_coordinates(m,2);
    zo = dom_int_coordinates(m,3);
    Velocities_domain = zeros(3,1);

    for n = 1:N_total
        x = int_coordinates(n,1);
        y = int_coordinates(n,2);
        z = int_coordinates(n,3);
        % boundary 1 -
        if (n >= 1) && (n <= Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = y + y_diff_start;
            lu2 = y - y_diff_start;
            k = 3;
            for i = 1:3
                j = 1;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) + ...
                        Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                j = 2;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) + ...
                        Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                j = 3;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) + ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                
            end
            % for boundary 2 -
        elseif(n> Nx*Ny) && (n <= Nx*Nz + Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;
            for i = 1:3
                j = 1;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2).*R_mat(n,j) - ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,2,1,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                j = 2;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2).*R_mat(n,j) - ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,2,1,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                j = 3;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2).*R_mat(n,j) - ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,2,1,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                 
            end
            % for boundary 3 -
        elseif(n> Nx*Nz + Nx*Ny) && (n <= Nx*Nz + 2*Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = y + y_diff_start;
            lu2 = y - y_diff_start;
            k = 3;
            g = Local_Numbering(n);
            [B,N] = boundary_corner(g);
            N_3ext_3int = Nx*Ny + Nx*Nz + 2*Ny*Nz;
            del_x = 2.*x_diff;
            del_y = 2.*y_diff;
            del_x_d = 1.5*x_diff;
            del_y_d = 1.5*x_diff;
            if B == 2
                if N == 1 % top left corner
                    for i = 1:3
                        j = 1;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 2;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 3;
                        Velocities_domain(i) = Velocities_domain(i) + (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) - p.*(Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) + (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2).*R_mat(n+1,j)./del_x_d)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2).*R_mat(n+Nx,j)./del_y_d)./(8*pi*mu);
                    end
                elseif N == 2 % top right corner
                    for i = 1:3
                        j = 1;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 2;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 3;
                        Velocities_domain(i) = Velocities_domain(i) + (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) - p.*(Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2).*R_mat(n-1,j)./del_x_d)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2).*R_mat(n+Nx,j)./del_y_d)./(8*pi*mu);
                    end
                elseif N == 3 % bottom right corner
                    for i = 1:3
                        j = 1;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 2;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 3;
                        Velocities_domain(i) = Velocities_domain(i) + (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) - p.*(Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2).*R_mat(n-1,j)./del_x_d)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) + (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2).*R_mat(n-Nx,j)./del_y_d)./(8*pi*mu);
                    end

                elseif N == 4 % bottom left corner
                    for i = 1:3
                        j = 1;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 2;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 3;
                        Velocities_domain(i) = Velocities_domain(i) + (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) - p.*(Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) + (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2).*R_mat(n+1,j)./del_x_d)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) + (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2).*R_mat(n-Nx,j)./del_y_d)./(8*pi*mu);
                    end
                end
            elseif B == 1
                if N == 1 % left Boundary
                    for i = 1:3
                        j = 1;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 2;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 3;
                        Velocities_domain(i) = Velocities_domain(i) + (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) - p.*(Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) + (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2).*R_mat(n+1,j)./del_x_d)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) + (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2).*R_mat(n-Nx,j)./del_y)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2).*R_mat(n+Nx,j)./del_y)./(8*pi*mu);
                    end

                elseif N == 2 % top boundary
                    for i = 1:3
                        j = 1;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 2;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 3;
                        Velocities_domain(i) = Velocities_domain(i) + (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) - p.*(Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) + (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2).*R_mat(n+1,j)./del_x)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2).*R_mat(n-1,j)./del_x)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2).*R_mat(n+Nx,j)./del_y_d)./(8*pi*mu);
                    end

                elseif N == 3 % right boundary
                    for i = 1:3
                        j = 1;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 2;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 3;
                        Velocities_domain(i) = Velocities_domain(i) + (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) - p.*(Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2).*R_mat(n-1,j)./del_x_d)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) + (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2).*R_mat(n-Nx,j)./del_y)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2).*R_mat(n+Nx,j)./del_y)./(8*pi*mu);
                    end

                elseif N == 4 % bottom boundary
                    for i = 1:3
                        j = 1;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 2;
                        Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        j = 3;
                        Velocities_domain(i) = Velocities_domain(i) + (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) - p.*(Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)))./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) + (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2).*R_mat(n+1,j)./del_x)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) - (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2).*R_mat(n-1,j)./del_x)./(8*pi*mu);
                        Velocities_domain(i) = Velocities_domain(i) + (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2).*R_mat(n-Nx,j)./del_y_d)./(8*pi*mu);
                    end
                end
            elseif B == 0
                for i = 1:3
                    j = 1;
                    Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                    Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                    j = 2;
                    Velocities_domain(i) = Velocities_domain(i) + (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j))./(8*pi*mu);
                    Velocities_domain(i) = Velocities_domain(i) - (mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                    j = 3;
                    Velocities_domain(i) = Velocities_domain(i) + (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j)./h  - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,1,2,ll1,lu1,ll2,lu2).*R_mat(n,j) - p.*(Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)))./(8*pi*mu);
                    Velocities_domain(i) = Velocities_domain(i) - (2.*mu*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h).*R_mat(n + N_3ext_3int,j)./(8*pi*mu);
                    Velocities_domain(i) = Velocities_domain(i) + (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2).*R_mat(n+1,j)./del_x)./(8*pi*mu);
                    Velocities_domain(i) = Velocities_domain(i) - (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2).*R_mat(n-1,j)./del_x)./(8*pi*mu);
                    Velocities_domain(i) = Velocities_domain(i) + (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2).*R_mat(n-Nx,j)./del_y)./(8*pi*mu);
                    Velocities_domain(i) = Velocities_domain(i) - (mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2).*R_mat(n+Nx,j)./del_y)./(8*pi*mu);
                end
            end

            % for boundary 4 -
        elseif(n> Nx*Nz + 2*Nx*Ny) && (n <= 2*Nx*Nz + 2*Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;
            for i = 1:3
                j = 1;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2).*R_mat(n,j) + ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,2,1,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                j = 2;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2).*R_mat(n,j) + ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,2,1,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                j = 3;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2).*R_mat(n,j) + ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,2,1,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
            end

            % for boundary 5 -
        elseif(n> 2*Nx*Nz + 2*Nx*Ny) && (n <= 2*Nx*Nz + 2*Nx*Ny + Ny*Nz)
            ll1 = y - x_diff_start;
            lu1 = y + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;
            for i = 1:3
                j = 1;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2).*R_mat(n,j) + ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,1,2,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                j = 2;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2).*R_mat(n,j) + ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,1,2,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                j = 3;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2).*R_mat(n,j) + ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,1,2,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
            end

            % for boundary 6 -
        elseif(n> 2*Nx*Nz + 2*Nx*Ny + Ny*Nz) && (n <= 2*Nx*Nz + 2*Nx*Ny + 2*Ny*Nz)
            ll1 = y - x_diff_start;
            lu1 = y + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;
            for i = 1:3
                j = 1;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2).*R_mat(n,j) - ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,1,2,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                j = 2;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2).*R_mat(n,j) - ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,1,2,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
                j = 3;
                Velocities_domain(i) = Velocities_domain(i) + (Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2).*R_mat(n,j) - ...
                    Tijk_calc_int(x,xo,y,yo,z,zo,i,j,1,2,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j))./(8*pi*mu);
            end
        end
    end
    Velocities_domain_1(m) = Velocities_domain(1);
    Velocities_domain_2(m) = Velocities_domain(2);
    Velocities_domain_3(m) = Velocities_domain(3);
    
end
Ud_x = Velocities_domain_1;
Ud_y = Velocities_domain_2;
Ud_z = Velocities_domain_3;
domx = dom_int_coordinates(:,1);
domy = dom_int_coordinates(:,2);
domz = dom_int_coordinates(:,3);
%%
quiver3(domx,domz,domy,Ud_x,Ud_z,Ud_y)
