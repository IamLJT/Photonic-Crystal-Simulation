clear all;
clc;
close all;
% ���ɿռ������������ʵֵ����糣�� �ŵ��� ���� are all not 1 
epsilon0=(1/(36*pi))*1e-9;   %��0=8.854187817*10^-12F/ m  
% mu0=4*pi*1e-7;               %��0=4�С�10^-7ţ��/����^2
% c=3e+8;
%epsilon0=5.58;
mu0=4*pi*1e-7;               %��0=4�С�10^-7ţ��/����^2
c=3e+8;
% epsr_material=3.4;         %���ʵĽ�糣��epsr  11.56

delta0=1.98e-1*epsilon0;
% delta0=0;
dx=0.05e-6;    %��λm
dy=0.05e-6;
dz=0.05e-6;
dt=0.5*1/(c*sqrt((1/dx)^2+(1/dy)^2+(1/dz)^2));  %��ֵɫɢ��������dt

Npml=12;  %��������������PML�߽�Ŀ��

nx=199+2*Npml;  %���Ӿ�������x����
ny=199+2*Npml;
nz=31+2*Npml;
nt=600;   %����ʱ�䲽��
wav=1.3;  %���䲨��

eps_r=ones(nx,ny,nz);
mu_r=ones(nx,ny,nz);
delta_r=ones(nx,ny,nz);
% epsilon=12.5*epsilon0*eps_r;
mu=mu0*mu_r;
delta=delta_r*delta0;

D=40;L=10;

%---------------------------- ���νṹ ---------------------
%   ������һ�����
epsilon=epsilon0*eps_r;

%��ʼ���絼��
    sigmax=zeros(nx,ny,nz);
    sigmay=zeros(nx,ny,nz);
    sigmaz=zeros(nx,ny,nz);
%ϵ��
factor_n=3;
R0=1e-6;

sigmamax=(-log(R0)*(factor_n+1)*epsilon0*c)/(2*Npml*dx);

%   ��һ��Բ���εĽ��ʣ�������PML�߽�
median_x=(1+nx)/2;      %   ���ĵ�����
median_y=(1+ny)/2;      
for i=(1):(nx)
    for j=(1):(ny)
        for k=(1):(nz)
            if (i-median_x)^2+(j-median_y)^2<(median_x-1)^2+0.5
                epsilon(i,j,k)=12.5*epsilon0;
                delta(i,j,k)=delta0;
%                 ��������PML�����д���Բ���ռ��µ�PML��Ӧ������ѡȡ
                if (i-median_x)^2+(j-median_y)^2>(median_x-1-Npml)^2
                    sigmax(i,j,k)=sigmamax*((abs(i-median_x)-(nx+1-2*Npml)/2)/Npml)^factor_n;
                    sigmay(i,j,k)=sigmamax*((abs(j-median_y)-(ny+1-2*Npml)/2)/Npml)^factor_n;
                end
                if abs(k-(nz+1)/2)>(nz+1-2*Npml)/2
                    sigmaz(i,j,k)=sigmamax*((abs(k-(nz+1)/2)-(nz-2*Npml)/2)/Npml)^factor_n;
                end
            end
        end
    end
end

%   �ڽ��������ڿף������ĵ�Ϊ��׼
i=median_x-2*D;
j=median_y;
a=i;
b=j;
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

%-------------------------------------------------------------
%��ʼ���糡�ų� TM��
    Hz=zeros(nx,ny,nz);
    Hzx=zeros(nx,ny,nz);
    Hzy=zeros(nx,ny,nz);
    Ey=zeros(nx,ny,nz);
    Ex=zeros(nx,ny,nz);
    Ez=zeros(nx,ny,nz);
    Hy=zeros(nx,ny,nz);
    Hx=zeros(nx,ny,nz);
   
%�ŵ��ʾ��� obtained by Perfectly Matched Layer condition
    sigma_starx=(sigmax.*mu)./epsilon;%������������Ϊ�Ѿ�ȷ����sigmax,sigamy
    sigma_stary=(sigmay.*mu)./epsilon;
    sigma_starz=(sigmaz.*mu)./epsilon;
    
M2=2*dt./(2*epsilon+dt*sigmax);
M1=(2*epsilon-dt*sigmax)./(2*epsilon+dt*sigmax);
N1=(2*epsilon-dt*sigmay)./(2*epsilon+dt*sigmay);
N2=2*dt./(2*epsilon+dt*sigmay);

%Ex�ĵ���ϵ��
A1=(M1*(dt^2)-delta.^2.*M2.*N2)./(dt^2-delta.^2.*M2.*N2);
A2=(dt^2*M2/dy)./(dt^2-delta.^2.*M2.*N2);
A3=(-dt^2*M2/dz)./(dt^2-delta.^2.*M2.*N2);
A4=-(1i*delta.*(N1-1).*M2*dt)./(dt^2-delta.^2.*M2.*N2);
A5=-(1i*delta.*M2.*N2*dt/dz)./(dt^2-delta.^2.*M2.*N2);
A6=(1i*delta.*M2.*N2*dt/dx)./(dt^2-delta.^2.*M2.*N2);

%Ey�ĵ���ϵ��
B1=(N1*(dt^2)-delta.^2.*M2.*N2)./(dt^2-delta.^2.*M2.*N2);
B2=(dt^2*N2/dz)./(dt^2-delta.^2.*M2.*N2);
B3=(-dt^2*N2/dx)./(dt^2-delta.^2.*M2.*N2);
B4=(1i*delta.*N2.*(M1-1)*dt)./(dt^2-delta.^2.*M2.*N2);
B5=(1i*delta.*M2.*N2*dt/dy)./(dt^2-delta.^2.*M2.*N2);
B6=-(1i*delta.*M2.*N2*dt/dz)./(dt^2-delta.^2.*M2.*N2);

%Ez�ĵ���ϵ��
C1=(2*epsilon-dt*sigmaz)./(2*epsilon+dt*sigmaz);
C2=2*dt./((2*epsilon+dt*sigmaz)*dx);
C3=-2*dt./((2*epsilon+dt*sigmaz)*dy);

%Hx�ĵ���ϵ��
D1=(2*mu-dt*sigma_starx)./(2*mu+dt*sigma_starx);
D2=2*dt./((2*mu+dt*sigma_starx)*dz);
D3=-2*dt./((2*mu+dt*sigma_starx)*dy);

%Hy�ĵ���ϵ��
E1=(2*mu-dt*sigma_stary)./(2*mu+dt*sigma_stary);
E2=2*dt./((2*mu+dt*sigma_stary)*dx);
E3=-2*dt./((2*mu+dt*sigma_stary)*dz);

%Hz�ĵ���ϵ��
F1=(2*mu-dt*sigma_starz)./(2*mu+dt*sigma_starz);
F2=2*dt./((2*mu+dt*sigma_starz)*dy);
F3=-2*dt./((2*mu+dt*sigma_starz)*dx);

% Ex(nx/2-7:nx/2+7,ny/2-7:ny/2+7,3+Npml)=0.5;
% Sumx1=zeros(nx,ny,nt);    //����̫���ڴ�ռ�ù���
% Sumy1=zeros(nx,ny,nt);
% Sumx2=zeros(nx,ny,nt);
% Sumy2=zeros(nx,ny,nt);

% Sumx=zeros(nx,ny);
% Sumy=zeros(nx,ny);
% Sum1=zeros(nx,ny,600);
% Sum2=zeros(nx,ny,600);
Jx=zeros(nt);
Jy=zeros(nt);
Ix=zeros(nt);
Iy=zeros(nt);

% ѭ�����¿�ʼ
    for n=1:1:nt
        
        Ex(1:nx-1,1:ny-1,1:nz-1)=A1(1:nx-1,1:ny-1,1:nz-1).*Ex(1:nx-1,1:ny-1,1:nz-1)+A2(1:nx-1,1:ny-1,1:nz-1).*(Hz(1:nx-1,2:ny,1:nz-1)-Hz(1:nx-1,1:ny-1,1:nz-1))...
            +A3(1:nx-1,1:ny-1,1:nz-1).*(Hy(1:nx-1,1:ny-1,2:nz)-Hy(1:nx-1,1:ny-1,1:nz-1))+A4(1:nx-1,1:ny-1,1:nz-1).*Ey(1:nx-1,1:ny-1,1:nz-1)...
            +A5(1:nx-1,1:ny-1,1:nz-1).*(Hx(1:nx-1,1:ny-1,2:nz)-Hx(1:nx-1,1:ny-1,1:nz-1))+A6(1:nx-1,1:ny-1,1:nz-1).*(Hz(2:nx,1:ny-1,1:nz-1)-Hz(1:nx-1,1:ny-1,1:nz-1));
        Ey(1:nx-1,1:ny-1,1:nz-1)=B1(1:nx-1,1:ny-1,1:nz-1).*Ey(1:nx-1,1:ny-1,1:nz-1)+B2(1:nx-1,1:ny-1,1:nz-1).*(Hx(1:nx-1,1:ny-1,2:nz)-Hx(1:nx-1,1:ny-1,1:nz-1))...
            +B3(1:nx-1,1:ny-1,1:nz-1).*(Hz(2:nx,1:ny-1,1:nz-1)-Hz(1:nx-1,1:ny-1,1:nz-1))+B4(1:nx-1,1:ny-1,1:nz-1).*Ex(1:nx-1,1:ny-1,1:nz-1)...
            +B5(1:nx-1,1:ny-1,1:nz-1).*(Hz(1:nx-1,2:ny,1:nz-1)-Hz(1:nx-1,1:ny-1,1:nz-1))+B6(1:nx-1,1:ny-1,1:nz-1).*(Hy(1:nx-1,1:ny-1,2:nz)-Hy(1:nx-1,1:ny-1,1:nz-1));
        Ez(1:nx-1,1:ny-1,1:nz-1)=C1(1:nx-1,1:ny-1,1:nz-1).*Ez(1:nx-1,1:ny-1,1:nz-1)+C2(1:nx-1,1:ny-1,1:nz-1).*(Hy(2:nx,1:ny-1,1:nz-1)-Hy(1:nx-1,1:ny-1,1:nz-1))...
            +C3(1:nx-1,1:ny-1,1:nz-1).*(Hx(1:nx-1,2:ny,1:nz-1)-Hx(1:nx-1,1:ny-1,1:nz-1));
        
        Hx(2:nx,2:ny,2:nz)=D1(1:nx-1,1:ny-1,1:nz-1).*Hx(2:nx,2:ny,2:nz)+D2(1:nx-1,1:ny-1,1:nz-1).*(Ey(2:nx,2:ny,2:nz)-Ey(2:nx,2:ny,1:nz-1))...
            +D3(1:nx-1,1:ny-1,1:nz-1).*(Ez(2:nx,2:ny,2:nz)-Ez(2:nx,1:ny-1,2:nz));
        Hy(2:nx,2:ny,2:nz)=E1(1:nx-1,1:ny-1,1:nz-1).*Hy(2:nx,2:ny,2:nz)+E2(1:nx-1,1:ny-1,1:nz-1).*(Ez(2:nx,2:ny,2:nz)-Ez(1:nx-1,2:ny,2:nz))...
            +E3(1:nx-1,1:ny-1,1:nz-1).*(Ex(2:nx,2:ny,2:nz)-Ex(2:nx,2:ny,1:nz-1));
        Hz(2:nx,2:ny,2:nz)=F1(1:nx-1,1:ny-1,1:nz-1).*Hz(2:nx,2:ny,2:nz)+F2(1:nx-1,1:ny-1,1:nz-1).*(Ex(2:nx,2:ny,2:nz)-Ex(2:nx,1:ny-1,2:nz))...
            +F3(1:nx-1,1:ny-1,1:nz-1).*(Ey(2:nx,2:ny,2:nz)-Ey(1:nx-1,2:ny,2:nz));
   

       tstart=1;N_lambda=wav*1e-6/dx;
        
       Ex(floor(median_x+D/2-L/2):floor(median_x+D/2+L/2),floor(median_y+(D/2)/sqrt(3)-L/2):floor(median_y+(D/2)/sqrt(3)+L/2),2+Npml)=0.5*sin(((2*pi*(c/(dx*N_lambda))*(n-tstart)*dt-pi/2)));
%        Ex(2+D+2*L+Npml:2*D+Npml,2+D+2*L+Npml:2*D+Npml,2+Npml)=0.5*sin(((2*pi*(c/(dx*N_lambda))*(n-tstart)*dt-pi/2)));
%        Ey(2+D+2*L+Npml:2*D+Npml,2+D+2*L+Npml:2*D+Npml,3+Npml)=0.5*sin(((2*pi*(c/(dx*N_lambda))*(n-tstart)*dt)));
%         Ex(2+D+2*L+Npml:2*D+Npml,2+D+2*L+Npml:2*D+Npml,3+Npml)=0.5;
        
%        Sumx1(:,:,n)=Ex(:,:,3+Npml);
%        Sumy1(:,:,n)=Ey(:,:,3+Npml);
%        Sumx2(:,:,n)=Ex(:,:,17+Npml);
%        Sumy2(:,:,n)=Ey(:,:,17+Npml);
        
        %��ͼ
        Jx(n)=sqrt(Ex(floor(median_x+D/2),floor(median_y+(D/2)/sqrt(3)),4+Npml)*conj(Ex(floor(median_x+D/2),floor(median_y+(D/2)/sqrt(3)),4+Npml)));
        Jy(n)=sqrt(Ey(floor(median_x+D/2),floor(median_y+(D/2)/sqrt(3)),4+Npml)*conj(Ey(floor(median_x+D/2),floor(median_y+(D/2)/sqrt(3)),4+Npml)));
        
%         Jx(n)=sqrt(Ex(1.5*D+L+1+Npml,1.5*D+L+1+Npml,3+Npml).*conj(Ex(1.5*D+L+Npml+1,1.5*D+L+Npml+1,3+Npml)));
%         Jy(n)=sqrt(Ey(1.5*D+L+1+Npml,1.5*D+L+1+Npml,3+Npml).*conj(Ey(1.5*D+L+Npml+1,1.5*D+L+Npml+1,3+Npml)));
%         Ix(n)=sqrt(Ex(1.5*D+L+1+Npml,1.5*D+L+1+Npml,2+Npml).*conj(Ex(1.5*D+L+Npml+1,1.5*D+L+Npml+1,2+Npml)));
%         Iy(n)=sqrt(Ey(1.5*D+L+1+Npml,1.5*D+L+1+Npml,2+Npml).*conj(Ey(1.5*D+L+Npml+1,1.5*D+L+Npml+1,2+Npml)));
%         Sumx(30+Npml:50+Npml,30+Npml:50+Npml)=sqrt(Ex(30+Npml:50+Npml,30+Npml:50+Npml,8+Npml).*...
%             conj(Ex(30+Npml:50+Npml,30+Npml:50+Npml,8+Npml)));
%         Sumy(30+Npml:50+Npml,30+Npml:50+Npml)=sqrt(Ey(30+Npml:50+Npml,30+Npml:50+Npml,8+Npml).*...
%             conj(Ey(30+Npml:50+Npml,30+Npml:50+Npml,8+Npml)));
%         Sumx(:,:)=sqrt(Ex(:,:,9+Npml).*conj(Ex(:,:,9+Npml)));
%         Sumy(:,:)=sqrt(Ey(:,:,9+Npml).*conj(Ey(:,:,9+Npml)));
        
         if(n>100)
            h=imagesc(1e+6*dy*(1:1:ny),(dx*1e+6*(1:1:nx))',abs(Ex(:,:,4+Npml))',[0,0.5]);colorbar;
            set(h,'AlphaData',10*epsilon(:,:,4+Npml)'/epsilon0);
%             Sum1(:,:,n-3000)=Sumx(:,:);
%             Sum2(:,:,n-3000)=Sumy(:,:);
            shading interp;
            getframe;
         end
    end
%     figure;
%     plot(Jx(:,1),Jy(:,1));
%     figure;
%     plot(Ix(:,1),Iy(:,1));

figure;
plot(Jx(1:nt),Jy(1:nt));
% h=imagesc(1e+6*dy*(1:1:ny),(dx*1e+6*(1:1:nx))',epsilon(:,:,1+Npml)',[0,1]);colorbar;
% set(h,'AlphaData',2*epsilon(:,:,8+Npml)'/epsilon0);