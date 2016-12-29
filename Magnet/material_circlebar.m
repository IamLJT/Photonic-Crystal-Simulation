clear all;
clc;
close all;
% 自由空间参数给出了真实值，介电常数 磁导率 光速 are all not 1 
epsilon0=(1/(36*pi))*1e-9;   %ε0=8.854187817*10^-12F/ m  
% mu0=4*pi*1e-7;               %μ0=4π×10^-7牛顿/安培^2
% c=3e+8;
%epsilon0=5.58;
mu0=4*pi*1e-7;               %μ0=4π×10^-7牛顿/安培^2
c=3e+8;
% epsr_material=3.4;         %介质的介电常数epsr  11.56

delta0=1.98e-1*epsilon0;
% delta0=0;
dx=0.05e-6;    %单位m
dy=0.05e-6;
dz=0.05e-6;
dt=0.5*1/(c*sqrt((1/dx)^2+(1/dy)^2+(1/dz)^2));  %数值色散条件计算dt

Npml=12;  %用网格数来定义PML边界的宽度

nx=199+2*Npml;  %光子晶体区域x方向
ny=199+2*Npml;
nz=31+2*Npml;

wav=1.3;  %入射波长

eps_r=ones(nx,ny,nz);
mu_r=ones(nx,ny,nz);
delta_r=ones(nx,ny,nz);
% epsilon=12.5*epsilon0*eps_r;
mu=mu0*mu_r;
delta=delta_r*delta0;

D=40;L=10;

%---------------------------- 柱形结构 ---------------------
%   首先是一层空气
epsilon=12.5*epsilon0*eps_r;
median_x=(1+nx)/2;      %   中心点坐标
median_y=(1+ny)/2;      
i=median_x-2*D;
j=median_y;
n_layer=0;
while (i-median_x+0.5*D*n_layer)^2+(j-median_y+n_layer*D*sqrt(3)/2)^2<(median_x-1-L)^2+0.5
    n_hole=0;
    while (i-median_x+n_hole*D+0.5*D*n_layer)^2+(j-median_y+n_layer*D*sqrt(3)/2)^2<(median_x-1-L)^2+0.5
        for k=(1+Npml):(nz-Npml)
            epsilon=hole_matrix(epsilon,floor(i+n_hole*D+0.5*D*n_layer),floor(j+n_layer*D*sqrt(3)/2),k,L,epsilon0);
            delta=hole_matrix(delta,floor(i+n_hole*D+0.5*D*n_layer),floor(j+n_layer*D*sqrt(3)/2),k,L,0);
        end
        n_hole=n_hole+1;
    end
    n_hole=0;
    while (i-median_x+n_hole*D+0.5*D*n_layer)^2+(j-median_y-n_layer*D*sqrt(3)/2)^2<(median_x-1-L)^2+0.5
        for k=(1+Npml):(nz-Npml)
            epsilon=hole_matrix(epsilon,floor(i+n_hole*D+0.5*D*n_layer),floor(j-n_layer*D*sqrt(3)/2),k,L,epsilon0);
            delta=hole_matrix(delta,floor(i+n_hole*D+0.5*D*n_layer),floor(j-n_layer*D*sqrt(3)/2),k,L,0);
        end
        n_hole=n_hole+1;
    end
    n_layer=n_layer+1;
end

n_Radiu=0;
while median_x-n_Radiu*D-L>Npml
    n_Radiu=n_Radiu+1;
end

i=median_x-(n_Radiu-1)*D;
j=median_y;
n_layer=0;
while floor(n_layer*D*sqrt(3)+L)<(ny-1)/2-Npml
   n_hole=0;
   while i+n_hole*D+L<nx-Npml
       for k=(1+Npml):(nz-Npml)
           epsilon=hole_matrix(epsilon,floor(i+n_hole*D),floor(j+sqrt(3)*D*n_layer),k,L,epsilon0);
           delta=hole_matrix(delta,floor(i+n_hole*D),floor(j+sqrt(3)*D*n_layer),k,L,0);
           epsilon=hole_matrix(epsilon,floor(i+n_hole*D),floor(j-sqrt(3)*D*n_layer),k,L,epsilon0);
           delta=hole_matrix(delta,floor(i+n_hole*D),floor(j-sqrt(3)*D*n_layer),k,L,0);
       end
       n_hole=n_hole+1;
   end
   n_layer=n_layer+1;
end

n_layer=0;
if i-D/2-L>=Npml
    i=i-D/2;
else
    i=i+D/2;
end
while floor(D*sqrt(3)/2+n_layer*D*sqrt(3))<(ny-1)/2-Npml
    n_hole=0;
    while i+n_hole*D+L<nx-Npml
        for k=(1+Npml):(nz-Npml)
           epsilon=hole_matrix(epsilon,floor(i+n_hole*D),floor(j+D*sqrt(3)/2+sqrt(3)*D*n_layer),k,L,epsilon0);
           delta=hole_matrix(delta,floor(i+n_hole*D),floor(j+D*sqrt(3)/2+sqrt(3)*D*n_layer),k,L,0);
           epsilon=hole_matrix(epsilon,floor(i+n_hole*D),floor(j-D*sqrt(3)/2-sqrt(3)*D*n_layer),k,L,epsilon0);
           delta=hole_matrix(delta,floor(i+n_hole*D),floor(j-D*sqrt(3)/2-sqrt(3)*D*n_layer),k,L,0);
        end
        n_hole=n_hole+1;
    end
    n_layer=n_layer+1;
end