% extract the velocity profile of the last internal domain layer, as well
% as the coordinates of those
Vb3= zeros(Nx*Ny,3);
domx_3 = zeros(Nx*Ny,1);
domy_3 = zeros(Nx*Ny,1);
domz_3 = zeros(Nx*Ny,1);

for n = 1:total_domain
    m = (Nz-1)*Nx*Ny;
if n > m
    Vb3(n-m,1) = Ud_x(n);
    Vb3(n-m,2) = Ud_y(n);
    Vb3(n-m,3) = Ud_z(n);
    domx_3(n-m) = domx(n);
    domy_3(n-m) = domy(n);
    domz_3(n-m) = domz(n);

end
end
 
% figure
% quiver3(domx_3,domz_3,domy_3,Vb3(:,1),Vb3(:,3),Vb3(:,2))
% title("EXIT PROFILE")

% create separate velocity matrices for each of the directional velocities
Ud_x_3 = Vb3(:,1);
Ud_y_3 = Vb3(:,2);
Ud_z_3 = Vb3(:,3);

i = 0;
for m = 1:Ny
    t = Nx/2 + i*Nx;
    p1(m) = Ud_z_3(t);
    i = i +1  ;
end


Umx_center = (p1(Ny/2) + p1(Ny/2 +1))/2;

p1 = p1./Umx_center;

p2 = zeros(Ny + 2,1);
for n = 1:Ny
    p2(n+1) = p1(n);
end

bound_y_3 = zeros(18,1);
for p = 1:Ny
   bound_y_3(p+ 1) = y_diff_start + y_diff*(p-1);
end
bound_y_3(Ny+2)= Ly;
bound_y_3 = bound_y_3 - Ly/2;

bound_y_3 = 2*bound_y_3/Ly;

% create the curve for the results
figure
plot(bound_y_3,p2,"o")
title("PROFILE CURVE","Color",'r')
ylabel("U")
hold on


for m = 1:Nx
    t = Nx*Ny/2;
    p1x(m) = Ud_z_3(t+m);
end


Umx_center_2 = (p1x(Nx/2) + p1x(Nx/2 +1))/2;

p1x = p1x./Umx_center_2;

p2x = zeros(Nx + 2,1);
for n = 1:Nx
    p2x(n+1) = p1x(n);
end

bound_x_3 = zeros(Nx,1);
for p = 1:Nx
   bound_x_3(p+ 1) = x_diff_start + x_diff*(p-1);
end
bound_x_3(Nx+2)= Lx;
bound_x_3 = bound_x_3 - Lx/2;

bound_x_3 = 2*bound_x_3/Lx;
plot(bound_x_3,p2x,"o","Color",'r')

%quiver3(domx,domz,domy,Ud_x,Ud_z,Ud_y)

%% 
% the code ahead is to generate the y profile using the analytical equation given in the
% research paper
by = Ly;%giving the y height
hx = Lx;%giving the x breadth
zeta =@(y) 2.*y./by; % create the non dimensional functions of x,y
sai =@(x) 2.*x./hx;
x = 0;
y=  0;

P = 40; %iterations to run of the fourier summation.
V = zeros(101,1);
var1 = -0.5:0.01:0.5;
for o = 1:101
    x = var1(o);
    sum = 0;
    for n = 0:P
        k = 2*n+1;
        sum = sum + (-1)^(0.5*(k-1))*(k^(-3))*cos(0.5*k*pi*zeta(y))*(1-(cosh(0.5*k*pi*sai(x)*by/hx)/(cosh(0.5*k*pi*by/hx))));
    end

    sum2 = 0;
    for n = 0:P
        k = 2*n+1;
        sum2 = sum2 + (-1)^(0.5*(k-1))*(k^(-3))*(1-1/(cosh(0.5*k*pi*by/hx)));
    end
    V(o) = sum/sum2;
end

m = -0.5:0.01:0.5;
m = m./0.5;
plot(m,V)


%%
%%
% the code generates the x profile using the analytical equation given in
% the research paper
V2 = zeros(101,1);
var1 = -0.5:0.01:0.5;
for o = 1:101
    x = var1(o);
    sum = 0;
    for n = 0:P
        k = 2*n+1;
        sum = sum + (-1)^(0.5*(k-1))*(k^(-3))*cos(0.5*k*pi*sai(x))*(1-(cosh(0.5*k*pi*zeta(y)*by/hx)/(cosh(0.5*k*pi*by/hx))));
    end

    sum2 = 0;
    for n = 0:P
        k = 2*n+1;
        sum2 = sum2 + (-1)^(0.5*(k-1))*(k^(-3))*(1-1/(cosh(0.5*k*pi*by/hx)));
    end
    V2(o) = sum/sum2;
end

m2 = -0.5:0.01:0.5;
m2 = m2./0.5;
plot(m2,V2)


