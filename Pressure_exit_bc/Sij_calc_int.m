function d = Sij_calc_int(x,xo,y,yo,z,zo,i,j,int_var1,int_var2,ll1,lu1,ll2,lu2)
% here we integrate the function on a step basis, where the first
% integration will be done analytically and the second integration will be
% done numerically

if i == 1 && j == 1
    f = @(x,xo,y,yo,z,zo) -1./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^0.5 - (x-xo).*(x-xo)./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^1.5;
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,x_0,y,y_0,z,z_0) -log(abs(sqrt((x-x_0).^2+(z_0-z).^2+(y_0-y).^2)+x-x_0));
        d1_int = @(y) d1(lu1,xo,y,yo,z,zo) - d1(ll1,xo,y,yo,z,zo);
        d2 = @(x,x_0,y,y_0,z,z_0) (x-x_0)./sqrt((x-x_0).^2+(z_0-z).^2+(y_0-y).^2)-log(abs(sqrt((x-x_0).^2+(z_0-z).^2+(y_0-y).^2)+x-x_0));
        d2_int = @(y) d2(lu1,xo,y,yo,z,zo) - d2(ll1,xo,y,yo,z,zo);
        d = integral(@(y) d2_int(y),ll2,lu2) + integral(@(y) d1_int(y),ll2,lu2);
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) -log(abs(sqrt((x-x_0).^2+(z_0-z).^2+(y_0-y).^2)+x-x_0));
        d1_int = @(z) d1(lu1,xo,y,yo,z,zo) - d1(ll1,xo,y,yo,z,zo);
        d2 = @(x,x_0,y,y_0,z,z_0) (x-x_0)./sqrt((x-x_0).^2+(z_0-z).^2+(y_0-y).^2)-log(abs(sqrt((x-x_0).^2+(z_0-z).^2+(y_0-y).^2)+x-x_0));
        d2_int = @(z) d2(lu1,xo,y,yo,z,zo) - d2(ll1,xo,y,yo,z,zo);
        d = integral(@(z) d2_int(z),ll2,lu2) + integral(@(z) d1_int(z),ll2,lu2);
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) -log(abs(sqrt((y-y_0).^2./((z_0-z).^2+(x_0-x).^2)+1)+(y-y_0)./sqrt((z_0-z).^2+(x_0-x).^2)))...
            -((x-x_0).^2.*(y-y_0))./(((z_0-z).^2+(x_0-x).^2).^(3./2).*sqrt((y-y_0).^2./((z_0-z).^2+(x_0-x).^2)+1));
        d1_int = @(z) d1(x,xo,lu1,yo,z,zo) - d1(x,xo,ll1,yo,z,zo);
        d = integral(@(z) d1_int(z),ll2,lu2);    
    end

end


if i == 1 && j == 2
    f = @(x,xo,y,yo,z,zo) -(x-xo).*(y-yo)./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^1.5;
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) ((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^0.5;
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) (y-y_0).*log((sqrt((z-z_0).^2+(y_0-y).^2+(x_0-x).^2)+z-z_0));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) (x-x_0).*log((sqrt((z-z_0).^2+(y_0-y).^2+(x_0-x).^2)+z-z_0));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 1 && j == 3
    f = @(x,xo,y,yo,z,zo) -(x-xo).*(z-zo)./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^1.5;
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,x_0,y,y_0,z,z_0) (z-z_0)*log(abs(sqrt((y-y_0)^2+(z_0-z)^2+(x_0-x)^2)+y-y_0));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) sqrt((z-z_0)^2+(y-y_0)^2+(x-x_0)^2);
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) (x-x_0)*log(abs(sqrt((y-y_0)^2+(z_0-z)^2+(x_0-x)^2)+y-y_0));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 2 && j == 1
    f = @(x,xo,y,yo,z,zo) -(x-xo).*(y-yo)./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^1.5;
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) ((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^0.5;
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) (y-y_0).*log((sqrt((z-z_0).^2+(y_0-y).^2+(x_0-x).^2)+z-z_0));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) (x-x_0).*log((sqrt((z-z_0).^2+(y_0-y).^2+(x_0-x).^2)+z-z_0));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 2 && j == 2
    f = @(x,xo,y,yo,z,zo) -1./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^0.5 - (y-yo).*(y-yo)./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^1.5;
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,x_0,y,y_0,z,z_0) (y-y_0)./sqrt((y-y_0).^2+(z_0-z).^2+(x_0-x).^2)-2.*log(abs(sqrt((y-y_0).^2+(z_0-z).^2+(x_0-x).^2)+y-y_0));
        d1_int = @(x) d1(x,xo,lu2,yo,z,zo) - d1(x,xo,ll2,yo,z,zo);
        d = integral(@(x) d1_int(x),ll1,lu1);
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) -log(abs(sqrt((x-x_0).^2./((z_0-z).^2+(y_0-y).^2)+1)+(x-x_0)./sqrt((z_0-z).^2+(y_0-y).^2)))...
    -((y-y_0).^2.*(x-x_0))./(((z_0-z).^2+(y_0-y).^2).^(3./2).*sqrt((x-x_0).^2./((z_0-z).^2+(y_0-y).^2)+1)) ;
        d1_int = @(z) d1(lu1,xo,y,yo,z,zo) - d1(ll1,xo,y,yo,z,zo);
        d = integral(@(z) d1_int(z),ll2,lu2);
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) (y-y_0)./sqrt((y-y_0).^2+(z_0-z).^2+(x_0-x).^2)-2.*log(abs(sqrt((y-y_0).^2+(z_0-z).^2+(x_0-x).^2)+y-y_0));
        d1_int = @(z) d1(x,xo,lu1,yo,z,zo) - d1(x,xo,ll1,yo,z,zo);
        d = integral(@(z) d1_int(z),ll2,lu2);
    end
end

if i == 2 && j == 3
    f = @(x,xo,y,yo,z,zo) -(y-yo).*(z-zo)./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^1.5;
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,x_0,y,y_0,z,z_0) (z-z_0)*log(abs(sqrt((x-x_0)^2+(z_0-z)^2+(y_0-y)^2)+x-x_0));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) (y-y_0)*log(abs(sqrt((x-x_0)^2+(z_0-z)^2+(y_0-y)^2)+x-x_0));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) sqrt((z-z_0)^2+(y-y_0)^2+(x-x_0)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 3 && j == 1
    f = @(x,xo,y,yo,z,zo) -(x-xo).*(z-zo)./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^1.5;
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,x_0,y,y_0,z,z_0) (z-z_0)*log(abs(sqrt((y-y_0)^2+(z_0-z)^2+(x_0-x)^2)+y-y_0));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) sqrt((z-z_0)^2+(y-y_0)^2+(x-x_0)^2);
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) (x-x_0)*log(abs(sqrt((y-y_0)^2+(z_0-z)^2+(x_0-x)^2)+y-y_0));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 3 && j == 2
    f = @(x,xo,y,yo,z,zo) -(y-yo).*(z-zo)./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^1.5;
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,x_0,y,y_0,z,z_0) (z-z_0)*log(abs(sqrt((x-x_0)^2+(z_0-z)^2+(y_0-y)^2)+x-x_0));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) (y-y_0)*log(abs(sqrt((x-x_0)^2+(z_0-z)^2+(y_0-y)^2)+x-x_0));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) sqrt((z-z_0)^2+(y-y_0)^2+(x-x_0)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 3 && j == 3
    f = @(x,xo,y,yo,z,zo) -1./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^0.5 - (z-zo).*(z-zo)./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^1.5;
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,x_0,y,y_0,z,z_0) -log(abs(sqrt((x-x_0).^2./((z_0-z).^2+(y_0-y).^2)+1)+(x-x_0)./sqrt((z_0-z).^2+(y_0-y).^2)))...
    -((z-z_0).^2.*(x-x_0))./(((z_0-z).^2+(y_0-y).^2).^(3./2).*sqrt((x-x_0).^2./((z_0-z).^2+(y_0-y).^2)+1));
        d1_int = @(y) d1(lu1,xo,y,yo,z,zo) - d1(ll1,xo,y,yo,z,zo);
        d = integral(@(y) d1_int(y),ll2,lu2);
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) -log(abs(sqrt((x-x_0).^2./((z_0-z).^2+(y_0-y).^2)+1)+(x-x_0)./sqrt((z_0-z).^2+(y_0-y).^2)))...
    -((z-z_0).^2.*(x-x_0))./(((z_0-z).^2+(y_0-y).^2).^(3./2).*sqrt((x-x_0).^2./((z_0-z).^2+(y_0-y).^2)+1));
        d1_int = @(z) d1(lu1,xo,y,yo,z,zo) - d1(ll1,xo,y,yo,z,zo);
        d = integral(@(z) d1_int(z),ll2,lu2);
    elseif int_var1 == 2 && int_var2 == 3
        %d1 = @(x,x_0,y,y_0,z,z_0) (z-z_0)./sqrt((z-z_0).^2+(y_0-y).^2+(x_0-x).^2)-2.*log(abs(sqrt((z-z_0).^2+(y_0-y).^2+(x_0-x).^2)+z-z_0));
        %d1_int = @(y) d1(x,xo,y,yo,lu2,zo) - d1(x,xo,y,yo,ll2,zo);
        %d = integral(@(y) d1_int(y),ll1,lu1);
        d = quad2d(@(a,b) f(x,xo,a,yo,b,zo),ll1,lu1,ll2,lu2);
    end

end