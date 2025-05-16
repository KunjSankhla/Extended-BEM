A1 = cell(N_total_int);

for m = 1:N_total_int
    for n = 1:N_total_int
        mat = zeros(3,3);
        A1{m,n} = mat;
    end
end

parfor m = 1:N_total_int
    fprintf("Processing MatA %d of %d",m,N_total_int)
    fprintf("\n")
    
    xo = int_coordinates(m,1);
    yo = int_coordinates(m,2);
    zo = int_coordinates(m,3);
    A = cell(N_total_int,1);
        for n = 1:N_total_int
            mat = zeros(3,3);
            A{n} = mat;
        end

    for n = 1:N_total_int
        x = int_coordinates(n,1);
        y = int_coordinates(n,2);
        z = int_coordinates(n,3);
        mat_tmp = zeros(3,3);


        % boundary 1 -
        if (n >= 1) && (n <= Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = y + y_diff_start;
            lu2 = y - y_diff_start;
            for i = 1:3
                j = 1;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2);
                j = 2;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2);
                j = 3;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2);

            end
            A{n} = mat_tmp;

            % for boundary 2 -
        elseif(n> Nx*Ny) && (n <= Nx*Nz + Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;

            for i = 1:3
                j = 1;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2);
                j = 2;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2);
                j = 3;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2);

            end
            A{n} = mat_tmp;



            % for boundary 3 -
        elseif(n> Nx*Nz + Nx*Ny) && (n <= Nx*Nz + 2*Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = y + y_diff_start;
            lu2 = y - y_diff_start;
            g = Local_Numbering(n);
            [B,N] = boundary_corner(g);
            mat_tmp = zeros(3,3);
            mat_tmp_x_u = zeros(3,3);
            mat_tmp_x_d = zeros(3,3);
            mat_tmp_y_u = zeros(3,3);
            mat_tmp_y_d = zeros(3,3);
            mat_tmp_in_3 = zeros(3,3);
            N_3ext_3int = Nx*Ny + Nx*Nz + 2*Ny*Nz;
            del_x = 2.*x_diff;
            del_y = 2.*y_diff;
            del_x_d = 1.5*x_diff;
            del_y_d = 1.5*x_diff;

            if B == 2
                if N == 1 % top left corner
                    for i = 1:3
                        j = 1;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 2;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 3;
                        mat_tmp(i,j) = 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        mat_tmp_x_u(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2)./del_x_d;
                        mat_tmp_y_d(i,j) = -mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2)./del_y_d;
                    end
                    A{n} = A{n} + mat_tmp ;
                    A{n+1} = A{n+1} + mat_tmp_x_u;
                    A{n+Nx} = A{n+Nx} + mat_tmp_y_d;
                    A{n+N_3ext_3int} = A{n+N_3ext_3int} + mat_tmp_in_3;


                elseif N == 2 % top right corner
                    for i = 1:3
                        j = 1;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 2;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 3;
                        mat_tmp(i,j) = 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        mat_tmp_x_d(i,j) = -mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2)./del_x_d;
                        mat_tmp_y_d(i,j) = -mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2)./del_y_d;
                    end
                    A{n} = A{n} + mat_tmp;
                    A{n-1} = A{n-1} + mat_tmp_x_d;
                    A{n+Nx} = A{n+Nx} + mat_tmp_y_d;
                    A{n+N_3ext_3int} = A{n+N_3ext_3int} + mat_tmp_in_3;


                elseif N == 3 % bottom right corner
                    for i = 1:3
                        j = 1;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 2;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 3;
                        mat_tmp(i,j) = 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        mat_tmp_x_d(i,j) = -mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2)./del_x_d;
                        mat_tmp_y_u(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2)./del_y_d;
                    end
                    A{n} = A{n} + mat_tmp;
                    A{n-1} = A{n-1} + mat_tmp_x_d;
                    A{n-Nx} = A{n-Nx} + mat_tmp_y_u;
                    A{n+N_3ext_3int} = A{n+N_3ext_3int} + mat_tmp_in_3;


                elseif N == 4 % bottom left corner
                    for i = 1:3
                        j = 1;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 2;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 3;
                        mat_tmp(i,j) = 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        mat_tmp_x_u(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2)./del_x_d;
                        mat_tmp_y_u(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2)./del_y_d;
                    end
                    A{n} = A{n} + mat_tmp;
                    A{n+1} = A{n+1} + mat_tmp_x_u;
                    A{n-Nx} = A{n-Nx} + mat_tmp_y_u;
                    A{n + N_3ext_3int} = A{n + N_3ext_3int} + mat_tmp_in_3;

                end
            elseif B == 1
                if N == 1 % left boundary
                    for i = 1:3
                        j = 1;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 2;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 3;
                        mat_tmp(i,j) = 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        mat_tmp_x_u(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2)./del_x_d;
                        mat_tmp_y_u(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2)./del_y;
                        mat_tmp_y_d(i,j) = -mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2)./del_y;
                    end
                    A{n} = A{n} + mat_tmp;
                    A{n+1} = A{n+1} + mat_tmp_x_u;
                    A{n+Nx} = A{n+Nx} + mat_tmp_y_d;
                    A{n-Nx} = A{n-Nx} + mat_tmp_y_u;
                    A{n+N_3ext_3int} = A{n+N_3ext_3int} + mat_tmp_in_3;

                elseif N == 2 % top boundary
                    for i = 1:3
                        j = 1;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 2;
                        mat_tmp(i,j) = mat_tmp(i,j) + mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 3;
                        mat_tmp(i,j) = 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        mat_tmp_x_u(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2)./del_x;
                        mat_tmp_x_d(i,j) = -mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2)./del_x;
                        mat_tmp_y_d(i,j) = -mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2)./del_y_d;
                    end
                    A{n} = A{n} + mat_tmp;
                    A{n+1} = A{n+1} + mat_tmp_x_u;
                    A{n-1} = A{n-1} + mat_tmp_x_d;
                    A{n+Nx} = A{n+Nx} + mat_tmp_y_d;
                    A{n + N_3ext_3int} = A{n + N_3ext_3int} + mat_tmp_in_3;


                elseif N == 3 % right boundary
                    for i = 1:3
                        j = 1;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 2;
                        mat_tmp(i,j) = mat_tmp(i,j) + mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 3;
                        mat_tmp(i,j) = 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        mat_tmp_x_d(i,j) = -mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2)./del_x_d;
                        mat_tmp_y_u(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2)./del_y;
                        mat_tmp_y_d(i,j) = -mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2)./del_y;
                    end
                    A{n} = A{n} + mat_tmp;
                    A{n-1} = A{n-1} + mat_tmp_x_d;
                    A{n+Nx} = A{n+Nx} + mat_tmp_y_d;
                    A{n-Nx} = A{n-Nx} + mat_tmp_y_u;
                    A{n+N_3ext_3int} = A{n+N_3ext_3int} + mat_tmp_in_3;


                elseif N == 4 % bottom boundary
                    for i = 1:3
                        j = 1;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 2;
                        mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        j = 3;
                        mat_tmp(i,j) = 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                        mat_tmp_in_3(i,j) = - 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                        mat_tmp_x_u(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2)./del_x;
                        mat_tmp_x_d(i,j) = -mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2)./del_x;
                        mat_tmp_y_u(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2)./del_y_d;

                    end
                    A{n} = A{n} + mat_tmp;
                    A{n+1} = A{n+1} + mat_tmp_x_u;
                    A{n-1} = A{n-1} + mat_tmp_x_d;
                    A{n-Nx} = A{n-Nx} + mat_tmp_y_u;
                    A{n+N_3ext_3int} = A{n+N_3ext_3int} + mat_tmp_in_3;

                end
            elseif B == 0
                for i = 1:3
                    j = 1;
                    mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                    mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                    j = 2;
                    mat_tmp(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                    mat_tmp_in_3(i,j) = - mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                    j = 3;
                    mat_tmp(i,j) = 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                    mat_tmp_in_3(i,j) = - 2.*mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2)./h;
                    mat_tmp_x_u(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2)./del_x;
                    mat_tmp_x_d(i,j) = -mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,1,1,2,ll1,lu1,ll2,lu2)./del_x;
                    mat_tmp_y_u(i,j) = mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2)./del_y;
                    mat_tmp_y_d(i,j) = -mu.*Sij_calc_int(x,xo,y,yo,z,zo,i,2,1,2,ll1,lu1,ll2,lu2)./del_y;
                end
                A{n} = A{n} + mat_tmp;
                A{n+1} = A{n+1} + mat_tmp_x_u;
                A{n-1} = A{n-1} + mat_tmp_x_d;
                A{n+Nx} = A{n+Nx} + mat_tmp_y_d;
                A{n-Nx} = A{n-Nx} + mat_tmp_y_u;
                A{n+N_3ext_3int} = A{n+N_3ext_3int} + mat_tmp_in_3;

            end
            if m == n
                for i = 1:3
                    A{n}(i,i) = A{n}(i,i) - 4*pi*mu;
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
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2);
                j = 2;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2);
                j = 3;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2);
            end
            A{n} = mat_tmp;

            % for boundary 5 -
        elseif(n> 2*Nx*Nz + 2*Nx*Ny) && (n <= 2*Nx*Nz + 2*Nx*Ny + Ny*Nz)
            ll1 = y - y_diff_start;
            lu1 = y + y_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;

            for i = 1:3
                j = 1;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2);
                j = 2;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2);
                j = 3;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2);
            end
            A{n} = mat_tmp;


            % for boundary 6 -
        elseif(n> 2*Nx*Nz + 2*Nx*Ny + Ny*Nz) && (n <= 2*Nx*Nz + 2*Nx*Ny + 2*Ny*Nz)
            ll1 = y - y_diff_start;
            lu1 = y + y_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;

            for i = 1:3
                j = 1;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2);
                j = 2;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2);
                j = 3;
                mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2);
            end
            A{n} = mat_tmp;

        elseif (n>N_total) && (n<=N_total_int) % for boundary internal part
            if m == n
                for i = 1:3

                    A{n}(i,i) = A{n}(i,i) - 8*pi*mu;

                end
            end
        end
    end

    for k = 1:N_total_int
        A1{m,k} = A{k};
    end

end

MatA = cell2mat(A1);
clear('A1');