B = cell(N_total_int,1);
for m = 1:N_total_int
    mat = zeros(3,1);
    B{m} = mat;
end


parfor m = 1:N_total_int
    fprintf("Processing MatB %d of %d",m,N_total_int)
    fprintf("\n")
    xo = int_coordinates(m,1);
    yo = int_coordinates(m,2);
    zo = int_coordinates(m,3);

    for n = 1:N_total_int
        x = int_coordinates(n,1);
        y = int_coordinates(n,2);
        z = int_coordinates(n,3);
        mat_tmp = zeros(3,1);


        % boundary 1 -
        if (n >= 1) && ( n <= Nx*Ny)
            ll1 = x - x_diff_start;
            lu1 = x + x_diff_start;
            ll2 = y + y_diff_start;
            lu2 = y - y_diff_start;
            for i = 1:3
                for j = 1:3
                    mat_tmp(i) =  mat_tmp(i) - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j);
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
                    mat_tmp(i) =  mat_tmp(i) + Tijk_calc_int(x,xo,y,yo,z,zo,i,j,2,1,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j);
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
                    mat_tmp(i) =  mat_tmp(i) + Tijk_calc_int(x,xo,y,yo,z,zo,i,j,3,1,2,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j);
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
                    mat_tmp(i) =  mat_tmp(i) - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,2,1,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j);
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
                    mat_tmp(i) =  mat_tmp(i) - Tijk_calc_int(x,xo,y,yo,z,zo,i,j,1,2,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j);
                end
            end


          %  for boundary 6 -
        elseif(n> 2*Nx*Nz + 2*Nx*Ny + Ny*Nz) && (n <= 2*Nx*Nz + 2*Nx*Ny + 2*Ny*Nz)
            ll1 = y - x_diff_start;
            lu1 = y + x_diff_start;
            ll2 = z + z_diff_start;
            lu2 = z - z_diff_start;

            for i = 1:3
                for j = 1:3
                    mat_tmp(i) =  mat_tmp(i) + Tijk_calc_int(x,xo,y,yo,z,zo,i,j,1,2,3,ll1,lu1,ll2,lu2).*Velocity_boundaries(n,j);
                end
            end

        end

        B{m} = B{m} + mat_tmp;
    end
end
for m = 1:N_total
    if m > Nx*Ny + Nx*Nz && m <= 2*Nx*Ny + Nx*Nz
    else
        for i  = 1:3
            B{m}(i) = B{m}(i) + 4*pi*mu*Velocity_boundaries(m,i);
        end
    end
end


MatB = cell2mat(B);