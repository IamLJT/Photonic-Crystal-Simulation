disp('2.初始化常用参量');

nxp1 = nx+1;  nxm1 = nx-1; nxm2 = nx-2;
nyp1 = ny+1;  nym1 = ny-1; nym2 = ny-2; 
nzp1 = nz+1;  nzm1 = nz-1; nzm2 = nz-2; 

Lattice_Num=5;
Lattice_Type=2;      % 1为方形晶格，2为六角晶格
Lattice_Shape=2;     % 介质形状，1为方形，2为圆柱形
Lattice_Style=2;     % 1为空气孔，2为介质柱
Lattice_Const=20;    % 晶格常数S
Lattice_FillingRatio=0.35;   % 晶格填充率
Lattice_IsPointDefect=1;    % 是否有点缺陷
Lattice_DefectPointX=1+Npml+Lattice_Const*floor(Lattice_Num/2)+...
    Lattice_FillingRatio*Lattice_Const;
Lattice_DefectPointY=1+Npml+Lattice_Const*floor(Lattice_Num/2)+...
    Lattice_FillingRatio*Lattice_Const;
Lattice_DefectSize=Lattice_FillingRatio*Lattice_Const;

%PML边界的FDTD更新
%初始化电场磁场 TM波
Ex=zeros(nx,nyp1,nzp1);
Ey=zeros(nxp1,ny,nzp1);
Ez=zeros(nxp1,nyp1,nz);
Hx=zeros(nxp1,ny,nz);
Hy=zeros(nx,nyp1,nz);
Hz=zeros(nx,ny,nzp1);

if Lattice_Style==1
    epsilon=epsr_material*epsilon0*eps_r;   % 基底为磁性材料
    mu=mu_material*mu0*mu_r;
    k=k_material*epsilon0*k_r;
end
if Lattice_Style==2
    epsilon=epsilon0*eps_r;
    mu=mu0*mu_r;
    k=k0*k_r;
end


% 光源参数
source_type=1;      %   1为平面波，2为高斯波
source_shape=2;     %   1为圆，2为方
source_layer=3+Npml;
if Lattice_IsPointDefect==0
    source_mid_x=floor((1+nx)/2)-Lattice_Const/2;
    source_mid_y=floor((1+ny)/2)-floor(sqrt(3)*Lattice_Const/6);
else
    source_mid_x=floor((1+nx)/2);
    source_mid_y=floor((1+ny)/2);
end
source_size=floor(Lattice_Const-2*Lattice_Const*Lattice_FillingRatio)*3;
dst_layer=14+Npml;