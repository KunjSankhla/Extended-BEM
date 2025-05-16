function d = Tijk_calc_int(x,xo,y,yo,z,zo,i,j,k,int_var1,int_var2,ll1,lu1,ll2,lu2)
mu = 1; % change mu here if it is to be changed for the code!!
% T = @(x,xo,y,yo,z,zo) 6.*mu.*(x-xo).*(y-yo).*(z-zo)./((x-xo).^2 + (y-yo).^2 + (z-zo).^2).^2.5;
% 6*(x-x0)*(y-y0)*(z-z0)/((x-x0)^2 + (y-y0)^2 + (z-z0)^2)^2.5
% if int_var1 == 1 && int_var2 == 2 % integration wrt x then y
%     d = integral2(@(a,b) T(a,xo,b,yo,z,zo),ll1,lu1,ll2,lu2);
% elseif int_var1 == 1 && int_var2 == 3 % integration wrt x then z
%     d = integral2(@(a,b) T(a,xo,y,yo,b,zo),ll1,lu1,ll2,lu2);
% elseif int_var1 == 2 && int_var2 == 3 % integration wrt y then z
%     d = integral2(@(a,b) T(x,xo,a,yo,b,zo),ll1,lu1,ll2,lu2);
% end
% the integrations for this case will be done analytically


if i == 1 && j == 1 && k == 1
    % 6*(x-x0)*(x-x0)*(x-x0)/((x-x0)^2 + (y-y0)^2 + (z-z0)^2)^2.5
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -4*log(abs(sqrt((y-y_0)^2+(z_0-z)^2+(x_0-x)^2)+y-y_0))-(2*(x_0-x)^2*(y-y_0))/(((z_0-z)^2+(x_0-x)^2)*sqrt((y-y_0)^2+(z_0-z)^2+(x_0-x)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -4.*log(abs(sqrt((z-zo).^2+(yo-y).^2+(xo-x).^2)+z-zo))-(2.*(xo-x).^2.*(z-zo))./(((yo-y).^2+(xo-x).^2).*...
    sqrt((z-zo).^2+(yo-y).^2+(xo-x).^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2.*(xo-x).^3.*(y-yo).*(3.*((y-yo).^2+(zo-z).^2+(xo-x).^2)-y.^2+2.*yo.*y-yo.^2))...
        ./(((zo-z).^2+(xo-x).^2).^2.*((y-yo).^2+(zo-z).^2+(xo-x).^2).^(3./2));
        d1_int = @(z) d1(x,xo,lu1,yo,z,zo) - d1(x,xo,ll1,yo,z,zo);
        d = integral(@(z) d1_int(z),ll2,lu2);
    end
end

if i == 1 && j == 1 && k == 2
    % 6*(x-x0)*(x-x0)*(y-y0)/((x-x0)^2 + (y-y0)^2 + (z-z0)^2)^2.5
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) log(abs(sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2)+xo-x))-log(abs(sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2)-xo+x))+...
            (2*x-2*xo)/sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2);
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) -(2.*(y_0-y).*(x-x_0).^3)./(((z_0-z).^2+(y_0-y).^2).*((x-x_0).^2+(z_0-z).^2+(y_0-y).^2).^(3/2));
        d1_int = @(z) d1(lu1,xo,y,yo,z,zo) - d1(ll1,xo,y,yo,z,zo);
        d = integral(@(z) d1_int(z),ll2,lu2);
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) -(2.*(x-x_0).^2)./((y-y_0).^2+(z-z_0).^2+(x-x_0).^2).^(3/2);
        d1_int = @(z) d1(x,xo,lu1,yo,z,zo) - d1(x,xo,ll1,yo,z,zo);
        d = integral(@(z) d1_int(z),ll2,lu2);
    end
end

if i == 1 && j == 1 && k == 3
     % 6*(x-x0)*(x-x0)*(z-z0)/((x-x0)^2 + (y-y0)^2 + (z-z0)^2)^2.5
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,x_0,y,y_0,z,z_0) -(2.*(z_0-z).*(x-x_0).^3)./(((z_0-z).^2+(y_0-y).^2).*((x-x_0).^2+(z_0-z).^2+(y_0-y).^2).^(3/2));
        d1_int = @(y) d1(lu1,xo,y,yo,z,zo) - d1(ll1,xo,y,yo,z,zo);
        d = integral(@(y) d1_int(y),ll2,lu2);
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+xo-x))-log(abs(sqrt((z-zo)^2+(yo-y)^2+...
            (xo-x)^2)-xo+x))+(2*x-2*xo)/sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2);
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(xo-x)^2*(y-yo))/(((z-zo)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end


if i == 1 && j == 2 && k == 1
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) log(abs(sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2)+xo-x))-log(abs(sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2)-xo+x))+...
            (2*x-2*xo)/sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2);
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) (2*(xo-x)*atan((abs(xo-x)*(z-zo))/((yo-y)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2))))/abs(xo-x)-...
            (2*(xo-x)*(yo-y)*(z-zo))/(((yo-y)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(xo-x)^2*(z-zo))/(((yo-y)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 1 && j == 2 && k == 2
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -2*(log(abs(sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2)+y-yo))+(yo-y)/sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y)^2*(z-zo))/(((yo-y)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) (2*(yo-y)*atan((abs(yo-y)*(z-zo))/((xo-x)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2))))/abs(yo-y)-(2*(xo-x)*(yo-y)*(z-zo))...
            /(((yo-y)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 1 && j == 2 && k == 3

    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -(2*(zo-z))/sqrt((y-yo)^2+(z-zo)^2+(x-xo)^2);
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y))/sqrt((z-zo)^2+(y-yo)^2+(x-xo)^2);
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(xo-x))/sqrt((z-zo)^2+(y-yo)^2+(x-xo)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 1 && j == 3 && k == 1
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) (2*(xo-x)*atan((abs(xo-x)*(y-yo))/((zo-z)*sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2))))/abs(xo-x)-...
            (2*(xo-x)*(zo-z)*(y-yo))/(((zo-z)^2+(xo-x)^2)*sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+xo-x))-log(abs(sqrt((z-zo)^2+(yo-y)^2+...
            (xo-x)^2)-xo+x))+(2*x-2*xo)/sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2);
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(xo-x)^2*(y-yo))/(((z-zo)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 1 && j == 3 && k == 2

    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -(2*(zo-z))/sqrt((y-yo)^2+(z-zo)^2+(x-xo)^2);
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y))/sqrt((z-zo)^2+(y-yo)^2+(x-xo)^2);
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(xo-x))/sqrt((z-zo)^2+(y-yo)^2+(x-xo)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 1 && j == 3 && k == 3
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -(2*(zo-z)^2*(y-yo))/(((zo-z)^2+(xo-x)^2)*sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -2*(log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+z-zo))+(zo-z)/sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) (2*(zo-z)*atan((abs(zo-z)*(y-yo))/((xo-x)*sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2))))/abs(zo-z)-...
            (2*(xo-x)*(zo-z)*(y-yo))/(((zo-z)^2+(xo-x)^2)*sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 2 && j == 1 && k == 1
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) log(abs(sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2)+xo-x))-log(abs(sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2)-xo+x))+...
            (2*x-2*xo)/sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2);
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) (2*(xo-x)*atan((abs(xo-x)*(z-zo))/((yo-y)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2))))/abs(xo-x)-...
            (2*(xo-x)*(yo-y)*(z-zo))/(((yo-y)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(xo-x)^2*(z-zo))/(((yo-y)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 2 && j == 1 && k == 2
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -2*(log(abs(sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2)+y-yo))+(yo-y)/sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y)^2*(z-zo))/(((yo-y)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) (2*(yo-y)*atan((abs(yo-y)*(z-zo))/((xo-x)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2))))/abs(yo-y)-(2*(xo-x)*(yo-y)*(z-zo))...
            /(((yo-y)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 2 && j == 1 && k == 3

    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -(2*(zo-z))/sqrt((y-yo)^2+(z-zo)^2+(x-xo)^2);
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y))/sqrt((z-zo)^2+(y-yo)^2+(x-xo)^2);
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(xo-x))/sqrt((z-zo)^2+(y-yo)^2+(x-xo)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 2 && j == 2 && k == 1
    % 6*(x-x0)*(y-y0)*(y-y0)/((x-x0)^2 + (y-y0)^2 + (z-z0)^2)^2.5
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -2*(log(abs(sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2)+y-yo))+(yo-y)/sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y)^2*(z-zo))/(((yo-y)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) -(2.*(x_0-x).*(y-y_0).^3)./(((z_0-z).^2+(x_0-x).^2).*((y-y_0).^2+(z_0-z).^2+(x_0-x).^2).^(3/2));
        d1_int = @(z) d1(x,xo,lu1,yo,z,zo) - d1(x,xo,ll1,yo,z,zo);
        d = integral(@(z) d1_int(z),ll2,lu2);
    end
end

if i == 2 && j == 2 && k == 2
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,x_0,y,y_0,z,z_0) -4*log(abs(sqrt((x-x_0)^2+(z_0-z)^2+(y_0-y)^2)+x-x_0))-...
            (2*(y_0-y)^2*(x-x_0))/(((z_0-z)^2+(y_0-y)^2)*sqrt((x-x_0)^2+(z_0-z)^2+(y_0-y)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2.*(yo-y).^3.*(x-xo).*(3.*((x-xo).^2+(zo-z).^2+(yo-y).^2)-x.^2+2.*xo.*x-xo.^2))./...
    (((zo-z).^2+(yo-y).^2).^2.*((x-xo).^2+(zo-z).^2+(yo-y).^2).^(3./2));
        d1_int = @(z) d1(lu1,xo,y,yo,z,zo) - d1(ll1,xo,y,yo,z,zo);
        d = integral(@(z) d1_int(z),ll2,lu2);
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -4*log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+z-zo))-(2*(yo-y)^2*(z-zo))/...
            (((yo-y)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end


if i == 2 && j == 2 && k == 3
    % 6*(y-y0)*(y-y0)*(z-z0)/((x-x0)^2 + (y-y0)^2 + (z-z0)^2)^2.5
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,x_0,y,y_0,z,z_0) -(2.*(z_0-z).*(y-y_0).^3)./(((z_0-z).^2+(x_0-x).^2).*((y-y_0).^2+(z_0-z).^2+(x_0-x).^2).^(3/2));
        d1_int = @(x) d1(x,xo,lu2,yo,z,zo) - d1(x,xo,ll2,yo,z,zo);
        d = integral(@(x) d1_int(x),ll1,lu1);
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y)^2*(x-xo))/(((zo-z)^2+(yo-y)^2)*sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+yo-y))-log(abs(sqrt((z-zo)^2+...
            (yo-y)^2+(xo-x)^2)-yo+y))+(2*y-2*yo)/sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 2 && j == 3 && k == 1

    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -(2*(zo-z))/sqrt((y-yo)^2+(z-zo)^2+(x-xo)^2);
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y))/sqrt((z-zo)^2+(y-yo)^2+(x-xo)^2);
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(xo-x))/sqrt((z-zo)^2+(y-yo)^2+(x-xo)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 2 && j == 3 && k == 2
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) (2*(yo-y)*atan((abs(yo-y)*(x-xo))/((zo-z)*sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2))))/...
            abs(yo-y)-(2*(yo-y)*(zo-z)*(x-xo))/(((zo-z)^2+(yo-y)^2)*sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y)^2*(x-xo))/(((zo-z)^2+(yo-y)^2)*sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+yo-y))-log(abs(sqrt((z-zo)^2+...
            (yo-y)^2+(xo-x)^2)-yo+y))+(2*y-2*yo)/sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 2 && j == 3 && k == 3
    % 6*(y-y0)*(z-z0)*(z-z0)/((x-x0)^2 + (y-y0)^2 + (z-z0)^2)^2.5
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,x_0,y,y_0,z,z_0) -(2.*(z-z_0).^2)./((y-y_0).^2+(z-z_0).^2+(x-x_0).^2).^(3/2);
        d1_int = @(x) d1(x,xo,lu2,yo,z,zo) - d1(x,xo,ll2,yo,z,zo);
        d = integral(@(x) d1_int(x),ll1,lu1);
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y)^2*(x-xo))/(((zo-z)^2+(yo-y)^2)*sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+yo-y))-log(abs(sqrt((z-zo)^2+...
            (yo-y)^2+(xo-x)^2)-yo+y))+(2*y-2*yo)/sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 3 && j == 1 && k == 1
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) (2*(xo-x)*atan((abs(xo-x)*(y-yo))/((zo-z)*sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2))))/abs(xo-x)-...
            (2*(xo-x)*(zo-z)*(y-yo))/(((zo-z)^2+(xo-x)^2)*sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+xo-x))-log(abs(sqrt((z-zo)^2+(yo-y)^2+...
            (xo-x)^2)-xo+x))+(2*x-2*xo)/sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2);
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(xo-x)^2*(y-yo))/(((z-zo)^2+(xo-x)^2)*sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 3 && j == 1 && k == 2

    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -(2*(zo-z))/sqrt((y-yo)^2+(z-zo)^2+(x-xo)^2);
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y))/sqrt((z-zo)^2+(y-yo)^2+(x-xo)^2);
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(xo-x))/sqrt((z-zo)^2+(y-yo)^2+(x-xo)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 3 && j == 1 && k == 3
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -(2*(zo-z)^2*(y-yo))/(((zo-z)^2+(xo-x)^2)*sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -2*(log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+z-zo))+(zo-z)/sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) (2*(zo-z)*atan((abs(zo-z)*(y-yo))/((xo-x)*sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2))))/abs(zo-z)-...
            (2*(xo-x)*(zo-z)*(y-yo))/(((zo-z)^2+(xo-x)^2)*sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 3 && j == 2 && k == 1

    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -(2*(zo-z))/sqrt((y-yo)^2+(z-zo)^2+(x-xo)^2);
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y))/sqrt((z-zo)^2+(y-yo)^2+(x-xo)^2);
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(xo-x))/sqrt((z-zo)^2+(y-yo)^2+(x-xo)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 3 && j == 2 && k == 2
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) (2*(yo-y)*atan((abs(yo-y)*(x-xo))/((zo-z)*sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2))))/...
            abs(yo-y)-(2*(yo-y)*(zo-z)*(x-xo))/(((zo-z)^2+(yo-y)^2)*sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -(2*(yo-y)^2*(x-xo))/(((zo-z)^2+(yo-y)^2)*sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+yo-y))-log(abs(sqrt((z-zo)^2+...
            (yo-y)^2+(xo-x)^2)-yo+y))+(2*y-2*yo)/sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 3 && j == 2 && k == 3
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -(2*(zo-z)^2*(x-xo))/(((zo-z)^2+(yo-y)^2)*sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) (2*(zo-z)*atan((abs(zo-z)*(x-xo))/((yo-y)*sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2))))/...
            abs(zo-z)-(2*(yo-y)*(zo-z)*(x-xo))/(((zo-z)^2+(yo-y)^2)*sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -2*(log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+z-zo))+(zo-z)/sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end


if i == 3 && j == 3 && k == 1
    % 6*(x-x0)*(z-z0)*(z-z0)/((x-x0)^2 + (y-y0)^2 + (z-z0)^2)^2.5
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -(2*(zo-z)^2*(y-yo))/(((zo-z)^2+(xo-x)^2)*sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2));
        d = d1(lu1,xo,lu2,yo,z,zo) - d1(lu1,xo,ll2,yo,z,zo) - (d1(ll1,xo,lu2,yo,z,zo) - d1(ll1,xo,ll2,yo,z,zo));
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -2*(log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+z-zo))+(zo-z)/sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) -(2.*(x_0-x).*(z-z_0).^3)./(((y_0-y).^2+(x_0-x).^2).*((z-z_0).^2+(y_0-y).^2+(x_0-x).^2).^(3/2));
        d1_int = @(y) d1(x,xo,y,yo,lu2,zo) - d1(x,xo,y,yo,ll2,zo);
        d = integral(@(y) d1_int(y),ll1,lu1);
    end
end


if i == 3 && j == 3 && k == 2
    % 6*(y-y0)*(z-z0)*(z-z0)/((x-x0)^2 + (y-y0)^2 + (z-z0)^2)^2.5
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,x_0,y,y_0,z,z_0) -(2*(z-z_0)^2)/((y-y_0)^2+(z-z_0)^2+(x-x_0)^2)^(3/2);
        d1_int = @(x) d1(x,xo,lu2,yo,z,zo) - d1(x,xo,ll2,yo,z,zo);
        d = integral(@(x) d1_int(x),ll1,lu1);
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,x_0,y,y_0,z,z_0) -(2.*(y_0-y).*(z-z_0).^3)./(((y_0-y).^2+(x_0-x).^2).*((z-z_0).^2+(y_0-y).^2+(x_0-x).^2).^(3/2));
        d1_int = @(x) d1(x,xo,y,yo,lu2,zo) - d1(x,xo,y,yo,ll2,zo);
        d = integral(@(x) d1_int(x),ll1,lu1);
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) log(abs(sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2)+yo-y))-log(abs(sqrt((z-zo)^2+...
            (yo-y)^2+(xo-x)^2)-yo+y))+(2*y-2*yo)/sqrt((z-zo)^2+(yo-y)^2+(xo-x)^2);
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

if i == 3 && j == 3 && k == 3
    if int_var1 == 1 && int_var2 == 2
        d1 = @(x,xo,y,yo,z,zo) -(2.*(zo-z).^3.*(x-xo).*(3.*((x-xo).^2+(zo-z).^2+(yo-y).^2)-x.^2+2*xo*x-xo.^2))./...
            (((zo-z).^2+(yo-y).^2).^2.*((x-xo).^2+(zo-z).^2+(yo-y).^2).^(3/2));
        d1_int = @(y) d1(lu1,xo,y,yo,z,zo) - d1(ll1,xo,y,yo,z,zo);
        d = integral(@(y) d1_int(y),ll2,lu2);
    elseif int_var1 == 1 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -4*log(abs(sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2)+x-xo))-(2*(zo-z)^2*(x-xo))/...
            (((zo-z)^2+(yo-y)^2)*sqrt((x-xo)^2+(zo-z)^2+(yo-y)^2));
        d = d1(lu1,xo,y,yo,lu2,zo) - d1(lu1,xo,y,yo,ll2,zo) - (d1(ll1,xo,y,yo,lu2,zo) - d1(ll1,xo,y,yo,ll2,zo));
    elseif int_var1 == 2 && int_var2 == 3
        d1 = @(x,xo,y,yo,z,zo) -4*log(abs(sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2)+y-yo))-(2*(zo-z)^2*(y-yo))/...
            (((zo-z)^2+(xo-x)^2)*sqrt((y-yo)^2+(zo-z)^2+(xo-x)^2));
        d = d1(x,xo,lu1,yo,lu2,zo) - d1(x,xo,lu1,yo,ll2,zo) - (d1(x,xo,ll1,yo,lu2,zo) - d1(x,xo,ll1,yo,ll2,zo));
    end
end

d = mu*d;


