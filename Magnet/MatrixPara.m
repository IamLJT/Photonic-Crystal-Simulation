disp('4.迭代矩阵参数计算');

%初始化电导率
sigmax=zeros(nx,ny,nz);
sigmay=zeros(nx,ny,nz);
sigmaz=zeros(nx,ny,nz);
%系数
factor_n=3;
R0=1e-8;

sigmamax=(-log(R0)*(factor_n+1)*epsilon0*epsr_material*c)/(2*Npml*dx);

%PML传导率分布函数    
for i=1:1:nx
      for x=1:1:Npml
         for j=1:1:nz
             sigmay(i,Npml+1-x,j)=sigmamax*(x/Npml)^factor_n;
             sigmay(i,ny-Npml+x,j)=sigmamax*(x/Npml)^factor_n;
         end
      end
end
for i=1:1:ny
      for y=1:1:Npml
          for j=1:1:nz
              sigmax(Npml+1-y,i,j)=sigmamax*(y/Npml)^factor_n;
              sigmax(nx-Npml+y,i,j)=sigmamax*(y/Npml)^factor_n;
          end
      end
end
for i=1:1:nx
    for z=1:1:Npml
        for j=1:1:ny
            sigmaz(i,j,Npml+1-z)=sigmamax*(z/Npml)^factor_n;
            sigmaz(i,j,nz-Npml+z)=sigmamax*(z/Npml)^factor_n;
        end
    end
end
    
%磁导率矩阵 obtained by Perfectly Matched Layer condition
sigma_starx=(sigmax.*mu)./epsilon;%本构条件，因为已经确定了sigmax,sigamy
sigma_stary=(sigmay.*mu)./epsilon;
sigma_starz=(sigmaz.*mu)./epsilon;

%Ex的迭代系数
A1=-(4*epsilon.^2-2*sigmay.*epsilon*dt+2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2-4*k.^2)./...
    (-4*epsilon.^2-2*sigmay.*epsilon*dt-2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2+4*k.^2);
A2=-1i*k.*sigmay*4*dt./(-4*epsilon.^2-2*sigmay.*epsilon*dt-2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2+4*k.^2);
A3=-(2*epsilon+sigmay*dt)*4*dt./(2*dy*(-4*epsilon.^2-2*sigmay.*epsilon*dt-2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2+4*k.^2));
A4=(2*epsilon+sigmay*dt)*4*dt./(2*dz*(-4*epsilon.^2-2*sigmay.*epsilon*dt-2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2+4*k.^2));
A5=1i*k*4*dt./(dz*(-4*epsilon.^2-2*sigmay.*epsilon*dt-2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2+4*k.^2));
A6=-1i*k*4*dt./(dx*(-4*epsilon.^2-2*sigmay.*epsilon*dt-2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2+4*k.^2));

%Ey的迭代系数
B1=-(4*epsilon.^2-2*sigmax.*epsilon*dt+2*sigmay.*epsilon*dt-sigmax.*sigmay*dt^2-4*k.^2)./...
    (-4*epsilon.^2-2*sigmay.*epsilon*dt-2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2+4*k.^2);
B2=1i*k.*sigmax*4*dt./(-4*epsilon.^2-2*sigmay.*epsilon*dt-2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2+4*k.^2);
B3=-1i*k*4*dt./(dy*(-4*epsilon.^2-2*sigmay.*epsilon*dt-2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2+4*k.^2));
B4=1i*k*4*dt./(dz*(-4*epsilon.^2-2*sigmay.*epsilon*dt-2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2+4*k.^2));
B5=-(2*epsilon+sigmax*dt)*4*dt./(2*dz*(-4*epsilon.^2-2*sigmay.*epsilon*dt-2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2+4*k.^2));
B6=(2*epsilon+sigmax*dt)*4*dt./(2*dx*(-4*epsilon.^2-2*sigmay.*epsilon*dt-2*sigmax.*epsilon*dt-sigmax.*sigmay*dt^2+4*k.^2));

%Ez的迭代系数
C1=(2*epsilon-sigmaz*dt)./(2*epsilon+sigmaz*dt);
C2=(2*dt)./((2*epsilon+sigmax*dt)*dx);
C3=-(2*dt)./((2*epsilon+sigmax*dt)*dy);

%Hx的迭代系数
D1=(2*mu0-sigma_starx*dt)./(2*mu0+sigma_starx*dt);
D2=0*mu;
D3=-(2*dt)./((2*mu0+sigma_starx*dt)*dy);
D4=(2*dt)./((2*mu0+sigma_starx*dt)*dz);
D5=0*mu;
D6=0*mu;

%Hy的迭代系数
E1=(2*mu0-sigma_stary*dt)./(2*mu0+sigma_stary*dt);
E2=0*mu;
E3=0*mu;
E4=0*mu;
E5=-(2*dt)./((2*mu0+sigma_stary*dt)*dz);
E6=(2*dt)./((2*mu0+sigma_stary*dt)*dx);

%Hz的迭代系数
F1=(2*mu0-dt*sigma_starz)./(2*mu0+dt*sigma_starz);
F2=(2*dt)./(dy*(2*mu0+dt*sigma_starz));
F3=-(2*dt)./(dx*(2*mu0+dt*sigma_starz));