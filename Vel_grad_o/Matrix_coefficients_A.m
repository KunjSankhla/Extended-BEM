A1= cell(N_total_int);

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
                for j = 1:3
                    mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2);
                end
            end
            A{n} = mat_tmp;

            % for boundary 2 -
        elseif(n> Nx*Ny) && (n <= Nx*Nz + Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;

            for i = 1:3
                for j = 1:3
                    mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2);
                end
            end
            A{n} = mat_tmp;



            % for boundary 3 -
        elseif(n> Nx*Nz + Nx*Ny) && (n <= Nx*Nz + 2*Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = y + y_diff_start;
            lu2 = y - y_diff_start;
            for i = 1:3
                for j = 1:3
                    mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,2,ll1,lu1,ll2,lu2);
                end
            end
            
            A{n} = mat_tmp;
            mat_tmp_in_vel  = zeros(3,3);
            for i = 1:3
                for j = 3
                    mat_tmp_in_vel(i,j) = - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2);
                end
            end
            if m == n
                mat_tmp_in_vel(3,3) = mat_tmp_in_vel(3,3) - 4*pi*mu; 
            end
            A{n+N3_N3in_df} = mat_tmp_in_vel;

            % for boundary 4 -
        elseif(n> Nx*Nz + 2*Nx*Ny) && (n <= 2*Nx*Nz + 2*Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;

            for i = 1:3
                for j = 1:3
                    mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,1,3,ll1,lu1,ll2,lu2);
                end
            end
            A{n} = mat_tmp;

            % for boundary 5 -
        elseif(n> 2*Nx*Nz + 2*Nx*Ny) && (n <= 2*Nx*Nz + 2*Nx*Ny + Ny*Nz)
            ll1 = y - x_diff_start;
            lu1 = y + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;

            for i = 1:3
                for j = 1:3
                    mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2);
                end
            end
            A{n} = mat_tmp;


            % for boundary 6 -
        elseif(n> 2*Nx*Nz + 2*Nx*Ny + Ny*Nz) && (n <= 2*Nx*Nz + 2*Nx*Ny + 2*Ny*Nz)
            ll1 = y - x_diff_start;
            lu1 = y + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;

            for i = 1:3
                for j = 1:3
                    mat_tmp(i,j) = Sij_calc_int(x,xo,y,yo,z,zo,i,j,2,3,ll1,lu1,ll2,lu2);
                end
            end
            A{n} = mat_tmp;
        elseif n > N_total
            if m == n
                A{n}(3,3) = A{n}(3,3) - 8*pi*mu;
            end
        end
    end
    for k = 1:N_total_int
        A1{m,k} = A{k};
    end
end



MatA = cell2mat(A1);