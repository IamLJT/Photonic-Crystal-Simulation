disp('3.设置光子晶体类型及结构');

D=Lattice_Const;L=Lattice_Const*Lattice_FillingRatio;                  % 晶格常数为150um
epsilon=epsr_material*epsilon0*eps_r;   % 基底为磁性材料
mu=mu_material*mu_r;
k=k_material*k_r;

if Lattice_Type==1
    if Lattice_Shape==1
        for i=(1+Npml):D:(1+(Lattice_Num-1)*D+Npml)                 
            for j=(1+Npml):D:(1+(Lattice_Num-1)*D+Npml) 
                for z=1+Npml:1:21+Npml
                    if Lattice_IsPointDefect==1
                        if i~=Lattice_DefectPointX-L||j~=Lattice_DefectPointY-L
                            epsilon(i:i+2*L,j:j+2*L,z)=epsilon0;%晶格数和晶格大小
                            mu(i:i+2*L,j:j+2*L,z)=mu0;
                            k(i:i+2*L,j:j+2*L,z)=0;
                        end
                    elseif Lattice_IsPointDefect==0
                        epsilon(i:i+2*L,j:j+2*L,z)=epsilon0;%晶格数和晶格大小
                        mu(i:i+2*L,j:j+2*L,z)=mu0;
                        k(i:i+2*L,j:j+2*L,z)=0;
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
                                        k(a,b,z)=0;
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
                                    k(a,b,z)=0;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
elseif Lattice_Type==2
    median_x=(1+nx)/2;      %   中心点坐标
    median_y=(1+ny)/2;      
    i=median_x-2*D;
    j=median_y;
    n_layer=0;
    while (i-median_x+0.5*D*n_layer)^2+(j-median_y+n_layer*D*sqrt(3)/2)^2<(median_x-1-L)^2+0.5
        n_hole=0;
        while (i-median_x+n_hole*D+0.5*D*n_layer)^2+(j-median_y+n_layer*D*sqrt(3)/2)^2<(median_x-1-L)^2+0.5
            if Lattice_IsPointDefect~=0&&n_layer==0&&n_hole*D+i-median_x==0
                n_hole=n_hole+1;
                continue;
            end
            for z=(1+Npml):(nz-Npml)
                epsilon=hole_matrix(epsilon,floor(i+n_hole*D+0.5*D*n_layer),floor(j+n_layer*D*sqrt(3)/2),z,L,epsilon0);
                mu=hole_matrix(mu,floor(i+n_hole*D+0.5*D*n_layer),floor(j+n_layer*D*sqrt(3)/2),z,L,mu0);
                k=hole_matrix(k,floor(i+n_hole*D+0.5*D*n_layer),floor(j+n_layer*D*sqrt(3)/2),z,L,0);
            end
            n_hole=n_hole+1;
        end
        n_hole=0;
        while (i-median_x+n_hole*D+0.5*D*n_layer)^2+(j-median_y-n_layer*D*sqrt(3)/2)^2<(median_x-1-L)^2+0.5
            if Lattice_IsPointDefect~=0&&n_layer==0&&n_hole*D+i-median_x==0
                n_hole=n_hole+1;
                continue;
            end
            for z=(1+Npml):(nz-Npml)
                epsilon=hole_matrix(epsilon,floor(i+n_hole*D+0.5*D*n_layer),floor(j+n_layer*D*sqrt(3)/2),z,L,epsilon0);
                mu=hole_matrix(mu,floor(i+n_hole*D+0.5*D*n_layer),floor(j+n_layer*D*sqrt(3)/2),z,L,mu0);
                k=hole_matrix(k,floor(i+n_hole*D+0.5*D*n_layer),floor(j+n_layer*D*sqrt(3)/2),z,L,0);
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
    n_layer=1;
    while floor(n_layer*D*sqrt(3)+L)<(ny-1)/2-Npml
       n_hole=0;
       while i+n_hole*D+L<nx-Npml
           for z=(1+Npml):(nz-Npml)
               epsilon=hole_matrix(epsilon,floor(i+n_hole*D),floor(j+sqrt(3)*D*n_layer),z,L,epsilon0);
               mu=hole_matrix(mu,floor(i+n_hole*D),floor(j+sqrt(3)*D*n_layer),z,L,mu0);
               k=hole_matrix(k,floor(i+n_hole*D),floor(j+sqrt(3)*D*n_layer),z,L,0);
               epsilon=hole_matrix(epsilon,floor(i+n_hole*D),floor(j-sqrt(3)*D*n_layer),z,L,epsilon0);
               mu=hole_matrix(mu,floor(i+n_hole*D),floor(j-sqrt(3)*D*n_layer),z,L,mu0);
               k=hole_matrix(k,floor(i+n_hole*D),floor(j-sqrt(3)*D*n_layer),z,L,0);
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
            for z=(1+Npml):(nz-Npml)
               epsilon=hole_matrix(epsilon,floor(i+n_hole*D),floor(j+D*sqrt(3)/2+sqrt(3)*D*n_layer),z,L,epsilon0);
               mu=hole_matrix(mu,floor(i+n_hole*D),floor(j+D*sqrt(3)/2+sqrt(3)*D*n_layer),z,L,mu0);
               k=hole_matrix(k,floor(i+n_hole*D),floor(j+D*sqrt(3)/2+sqrt(3)*D*n_layer),z,L,0);
               epsilon=hole_matrix(epsilon,floor(i+n_hole*D),floor(j-D*sqrt(3)/2-sqrt(3)*D*n_layer),z,L,epsilon0);
               mu=hole_matrix(mu,floor(i+n_hole*D),floor(j-D*sqrt(3)/2-sqrt(3)*D*n_layer),z,L,mu0);
               k=hole_matrix(k,floor(i+n_hole*D),floor(j-D*sqrt(3)/2-sqrt(3)*D*n_layer),z,L,0);
            end
            n_hole=n_hole+1;
        end
        n_layer=n_layer+1;
    end
end