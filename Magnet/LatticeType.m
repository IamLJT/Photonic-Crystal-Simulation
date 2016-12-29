disp('3.设置光子晶体类型及结构');

D=Lattice_Const;L=Lattice_Const*Lattice_FillingRatio;                  % 晶格常数为150um

% if Lattice_Style==1
%     epsilon=epsr_material*epsilon0*eps_r;   % 基底为磁性材料
%     mu=mu_material*mu0*mu_r;
%     k=k_material*mu0*k_r;
% end
% if Lattice_Style==2
%     epsilon=epsilon0*eps_r;
%     mu=mu0*mu_r;
%     k=k0*mu_r;
% end

nxp_s=1+Npml;
nxp_e=nx-Npml;
nyp_s=1+Npml;
nyp_e=ny-Npml;

% if Lattice_Style==2
%     epsilon=epsilon0*eps_r;
%     mu=mu0*mu_r;
%     k=k0*mu_r;
%     epsilon0=epsr_material*epsilon0;
%     mu0=mu_material*mu0;
%     k0=k_material*mu0;
% end

if Lattice_Type==1
    if Lattice_Shape==1
        for i=(1+Npml):D:(1+(Lattice_Num-1)*D+Npml)                 
            for j=(1+Npml):D:(1+(Lattice_Num-1)*D+Npml) 
                for z=1+Npml:1:21+Npml
                    if Lattice_IsPointDefect==1
                        if i~=Lattice_DefectPointX-L||j~=Lattice_DefectPointY-L
                            epsilon(i:i+2*L,j:j+2*L,z)=epsilon0;%晶格数和晶格大小
                            mu(i:i+2*L,j:j+2*L,z)=mu0;
                            k(i:i+2*L,j:j+2*L,z)=k0;
                        end
                    elseif Lattice_IsPointDefect==0
                        epsilon(i:i+2*L,j:j+2*L,z)=epsilon0;%晶格数和晶格大小
                        mu(i:i+2*L,j:j+2*L,z)=mu0;
                        k(i:i+2*L,j:j+2*L,z)=k0;
                    end
                end
            end
        end
    elseif Lattice_Shape==2
        for i=(1+Npml+L):D:(1+(Lattice_Num-1)*D+Npml+L)
            for j=(1+Npml+L):D:(1+(Lattice_Num-1)*D+Npml+L)
                for z=1+Npml:1:21+Npml
                    if Lattice_IsPointDefect==1
                        if i~=Lattice_DefectPointX-L||j~=Lattice_DefectPointY-L
                            for a=(i-L):(i+L)   %   实际半径为L+1个像素了
                                for b=(j-L):(j+L)
                                    if (a-i)^2+(b-j)^2<L^2+0.5
                                        epsilon(a,b,z)=epsilon0;
                                        mu(a,b,z)=mu0;
                                        k(a,b,z)=k0;
                                    end
                                end
                            end
                        end
                    elseif Lattice_IsPointDefect==0
                        for a=(i-L):(i+L)   %   实际半径为L+1个像素了
                            for b=(j-L):(j+L)
                                if (a-i)^2+(b-j)^2<L^2+0.5
                                    epsilon(a,b,z)=epsilon0;
                                    mu(a,b,z)=mu0;
                                    k(a,b,z)=k0;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
elseif Lattice_Type==2
    median_x=floor((1+nx)/2);      %   中心点坐标
    median_y=floor((1+ny)/2);  
    num_median_x_1=1;              %   介质个数
    num_median_x_2=0;
    num_median_y_1=1;
    num_median_y_2=0;
    
    %------计算上下左右有多少层----------
    %-------------x方向----------------
    n_hole=1;
    while median_x-n_hole*D+L>=nxp_s
        n_hole=n_hole+1;
    end
    num_median_x_1=num_median_x_1+n_hole-1;
    n_hole=1;
    while median_x+n_hole*D-L<=nxp_e
        n_hole=n_hole+1;
    end
    num_median_x_1=num_median_x_1+n_hole-1;
    
    n_hole=0;
    while median_x-0.5*D-n_hole*D+L>=nxp_s
        n_hole=n_hole+1;
    end
    num_median_x_2=num_median_x_2+n_hole;
    n_hole=0;
    while median_x+0.5*D+n_hole*D-L<=nxp_e
        n_hole=n_hole+1;
    end
    num_median_x_2=num_median_x_2+n_hole;
    %-------------y方向----------------
    n_hole=1;
    while floor(median_y-n_hole*D*sqrt(3)+L)>=nyp_s
        n_hole=n_hole+1;
    end
    num_median_y_1=num_median_y_1+n_hole-1;
    n_hole=1;
    while floor(median_y+n_hole*D*sqrt(3)-L)<=nyp_e
        n_hole=n_hole+1;
    end
    num_median_y_1=num_median_y_1+n_hole-1;
    
    n_hole=0;
    while floor(median_y-D*sqrt(3)/2-n_hole*D*sqrt(3)+L)>=nyp_s
        n_hole=n_hole+1;
    end
    num_median_y_2=num_median_y_2+n_hole;
    n_hole=0;
    while floor(median_y+D*sqrt(3)/2+n_hole*D*sqrt(3)-L)<=nyp_e
        n_hole=n_hole+1;
    end
    num_median_y_2=num_median_y_2+n_hole;
    %----------------------------------------------------------
    %-------------------------挖孔-----------------------------
    i=0;j=0;
    for i=1-floor(num_median_x_1/2)-1:num_median_x_1-floor(num_median_x_1/2)-1
        for j=1-floor(num_median_y_1/2)-1:num_median_y_1-floor(num_median_y_1/2)-1
            if Lattice_IsPointDefect==1&&i==0&&j==0
%                 for z=(4+Npml):(nz-Npml-3)
%                     if Lattice_Style==1
%                         epsilon=hole_matrix(epsilon,floor(median_x+i*D),floor(median_y+j*D*sqrt(3)),z,floor(0.5*L),epsilon0,nxp_s,nxp_e,nyp_s,nyp_e);
%                         mu=hole_matrix(mu,floor(median_x+i*D),floor(median_y+j*D*sqrt(3)),z,floor(0.5*L),mu0,nxp_s,nxp_e,nyp_s,nyp_e);
%                         k=hole_matrix(k,floor(median_x+i*D),floor(median_y+j*D*sqrt(3)),z,floor(0.5*L),k0,nxp_s,nxp_e,nyp_s,nyp_e);
%                     end
%                     if Lattice_Style==2
%                         epsilon=hole_matrix(epsilon,floor(median_x+i*D),floor(median_y+j*D*sqrt(3)),z,floor(0.5*L),epsilon0*epsr_material,nxp_s,nxp_e,nyp_s,nyp_e);
%                         mu=hole_matrix(mu,floor(median_x+i*D),floor(median_y+j*D*sqrt(3)),z,floor(0.5*L),mu0*mu_material,nxp_s,nxp_e,nyp_s,nyp_e);
%                         k=hole_matrix(k,floor(median_x+i*D),floor(median_y+j*D*sqrt(3)),z,floor(0.5*L),mu0*k_material,nxp_s,nxp_e,nyp_s,nyp_e);
%                     end
%                 end
                continue;
            end
            for z=(3+Npml):(nz-Npml-3)
                if Lattice_Style==1
                    epsilon=hole_matrix(epsilon,floor(median_x+i*D),floor(median_y+j*D*sqrt(3)),z,L,epsilon0,nxp_s,nxp_e,nyp_s,nyp_e);
                    mu=hole_matrix(mu,floor(median_x+i*D),floor(median_y+j*D*sqrt(3)),z,L,mu0,nxp_s,nxp_e,nyp_s,nyp_e);
                    k=hole_matrix(k,floor(median_x+i*D),floor(median_y+j*D*sqrt(3)),z,L,k0,nxp_s,nxp_e,nyp_s,nyp_e);
                end
                if Lattice_Style==2
                    epsilon=hole_matrix(epsilon,floor(median_x+i*D),floor(median_y+j*D*sqrt(3)),z,L,epsilon0*epsr_material,nxp_s,nxp_e,nyp_s,nyp_e);
                    mu=hole_matrix(mu,floor(median_x+i*D),floor(median_y+j*D*sqrt(3)),z,L,mu0*mu_material,nxp_s,nxp_e,nyp_s,nyp_e);
                    k=hole_matrix(k,floor(median_x+i*D),floor(median_y+j*D*sqrt(3)),z,L,mu0*k_material,nxp_s,nxp_e,nyp_s,nyp_e);
                end
            end
        end
    end
    for i=1-floor(num_median_x_2/2):num_median_x_2-floor(num_median_x_2/2)
        for j=1-floor(num_median_y_2/2):num_median_y_2-floor(num_median_y_2/2)
            for z=(Npml+3):(nz-Npml-3)
                if Lattice_Style==1
                    epsilon=hole_matrix(epsilon,floor(median_x-D/2+i*D),floor(median_y-D*sqrt(3)/2+j*D*sqrt(3)),z,L,epsilon0,nxp_s,nxp_e,nyp_s,nyp_e);
                    mu=hole_matrix(mu,floor(median_x-D/2+i*D),floor(median_y-D*sqrt(3)/2+j*D*sqrt(3)),z,L,mu0,nxp_s,nxp_e,nyp_s,nyp_e);
                    k=hole_matrix(k,floor(median_x-D/2+i*D),floor(median_y-D*sqrt(3)/2+j*D*sqrt(3)),z,L,k0,nxp_s,nxp_e,nyp_s,nyp_e);
                end
                if Lattice_Style==2
                    epsilon=hole_matrix(epsilon,floor(median_x-D/2+i*D),floor(median_y-D*sqrt(3)/2+j*D*sqrt(3)),z,L,epsilon0*epsr_material,nxp_s,nxp_e,nyp_s,nyp_e);
                    mu=hole_matrix(mu,floor(median_x-D/2+i*D),floor(median_y-D*sqrt(3)/2+j*D*sqrt(3)),z,L,mu0*mu_material,nxp_s,nxp_e,nyp_s,nyp_e);
                    k=hole_matrix(k,floor(median_x-D/2+i*D),floor(median_y-D*sqrt(3)/2+j*D*sqrt(3)),z,L,mu0*k_material,nxp_s,nxp_e,nyp_s,nyp_e);                    
                end
            end
        end
    end
end

% if Lattice_Style==2
%     epsilon0=epsilon0/epsr_material;
%     mu0=mu0/mu_material;
%     k0=0;
% end