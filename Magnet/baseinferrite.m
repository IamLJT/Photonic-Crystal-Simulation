clear all;clc;close all;

StructureInitial;
ParameterInit;
LatticeType;
MatrixPara;
UpdateCoefficient;

disp('5.开始计算');

Sumx=zeros(nx,ny);
Sumy=zeros(nx,ny);
Sumz=zeros(nx,ny);
Sumhx=zeros(nx,ny);
Sumhy=zeros(nx,ny);
Sumhz=zeros(nx,ny);

% Sumx_t=zeros(nx,ny,nt-1000);
% Sumy_t=zeros(nx,ny,nt-1000);
% Sumz_t=zeros(nx,ny,nt-1000);
% Sumhx_t=zeros(nx,ny,nt-1000);
% Sumhy_t=zeros(nx,ny,nt-1000);
% Sumhz_t=zeros(nx,ny,nt-1000);

% Sum1=zeros(nx,ny,600);
% Sum2=zeros(nx,ny,600);
Jx=zeros(dst_layer-source_layer+1,nt);
Jy=zeros(dst_layer-source_layer+1,nt);
% Jx=zeros(nt,1);
% Jy=zeros(nt,1);
% Ix=zeros(nt,1);
% Iy=zeros(nt,1);

% 循环更新开始
for n=1:1:nt

    source;
    
%     Ex(1:nx-1,1:ny-1,1:nz-1)=A1(1:nx-1,1:ny-1,1:nz-1).*Ex(1:nx-1,1:ny-1,1:nz-1)+A2(1:nx-1,1:ny-1,1:nz-1).*(Hz(1:nx-1,2:ny,1:nz-1)-Hz(1:nx-1,1:ny-1,1:nz-1))...
%         +A3(1:nx-1,1:ny-1,1:nz-1).*(Hy(1:nx-1,1:ny-1,2:nz)-Hy(1:nx-1,1:ny-1,1:nz-1));
%     Ey(1:nx-1,1:ny-1,1:nz-1)=B1(1:nx-1,1:ny-1,1:nz-1).*Ey(1:nx-1,1:ny-1,1:nz-1)+B2(1:nx-1,1:ny-1,1:nz-1).*(Hx(1:nx-1,1:ny-1,2:nz)-Hx(1:nx-1,1:ny-1,1:nz-1))...
%         +B3(1:nx-1,1:ny-1,1:nz-1).*(Hz(2:nx,1:ny-1,1:nz-1)-Hz(1:nx-1,1:ny-1,1:nz-1));
%     Ez(1:nx-1,1:ny-1,1:nz-1)=C1(1:nx-1,1:ny-1,1:nz-1).*Ez(1:nx-1,1:ny-1,1:nz-1)+C2(1:nx-1,1:ny-1,1:nz-1).*(Hy(2:nx,1:ny-1,1:nz-1)-Hy(1:nx-1,1:ny-1,1:nz-1))...
%         +C3(1:nx-1,1:ny-1,1:nz-1).*(Hx(1:nx-1,2:ny,1:nz-1)-Hx(1:nx-1,1:ny-1,1:nz-1));
%     
%     Hx(2:nx,2:ny,2:nz)=D1(1:nx-1,1:ny-1,1:nz-1).*Hx(2:nx,2:ny,2:nz)+D2(1:nx-1,1:ny-1,1:nz-1).*Hy(2:nx,2:ny,2:nz)...
%         +D3(1:nx-1,1:ny-1,1:nz-1).*(Ez(2:nx,2:ny,2:nz)-Ez(2:nx,1:ny-1,2:nz))+D4(1:nx-1,1:ny-1,1:nz-1).*(Ey(2:nx,2:ny,2:nz)-Ey(2:nx,2:ny,1:nz-1))...
%         +D5(1:nx-1,1:ny-1,1:nz-1).*(Ex(2:nx,2:ny,2:nz)-Ex(2:nx,2:ny,1:nz-1))+D6(1:nx-1,1:ny-1,1:nz-1).*(Ez(2:nx,2:ny,2:nz)-Ez(1:nx-1,2:ny,2:nz));
%     Hy(2:nx,2:ny,2:nz)=E1(1:nx-1,1:ny-1,1:nz-1).*Hy(2:nx,2:ny,2:nz)+E2(1:nx-1,1:ny-1,1:nz-1).*Hx(2:nx,2:ny,2:nz)...
%         +E3(1:nx-1,1:ny-1,1:nz-1).*(Ez(2:nx,2:ny,2:nz)-Ez(2:nx,1:ny-1,2:nz))+E4(1:nx-1,1:ny-1,1:nz-1).*(Ey(2:nx,2:ny,2:nz)-Ey(2:nx,2:ny,1:nz-1))...
%         +E5(1:nx-1,1:ny-1,1:nz-1).*(Ex(2:nx,2:ny,2:nz)-Ex(2:nx,2:ny,1:nz-1))+E6(1:nx-1,1:ny-1,1:nz-1).*(Ez(2:nx,2:ny,2:nz)-Ez(1:nx-1,2:ny,2:nz));
%     Hz(2:nx,2:ny,2:nz)=F1(1:nx-1,1:ny-1,1:nz-1).*Hz(2:nx,2:ny,2:nz)+F2(1:nx-1,1:ny-1,1:nz-1).*(Ex(2:nx,2:ny,2:nz)-Ex(2:nx,1:ny-1,2:nz))...
%         +F3(1:nx-1,1:ny-1,1:nz-1).*(Ey(2:nx,2:ny,2:nz)-Ey(1:nx-1,2:ny,2:nz));

    Ex(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)=A1(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*Ex(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)+A2(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*Ey(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)...
        +A3(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*(Hz(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)-Hz(pis+1:pie-1,pjs:pje-2,pks+1:pke-1))+A4(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*(Hy(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)-Hy(pis+1:pie-1,pjs+1:pje-1,pks:pke-2))...
        +A5(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*(Hx(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)-Hx(pis+1:pie-1,pjs+1:pje-1,pks:pke-2))+A6(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*(Hz(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)-Hz(pis:pie-2,pjs+1:pje-1,pks+1:pke-1));
    Ey(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)=B1(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*Ey(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)+B2(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*Ex(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)...
        +B3(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*(Hz(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)-Hz(pis+1:pie-1,pjs:pje-2,pks+1:pke-1))+B4(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*(Hy(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)-Hy(pis+1:pie-1,pjs+1:pje-1,pks:pke-2))...
        +B5(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*(Hx(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)-Hx(pis+1:pie-1,pjs+1:pje-1,pks:pke-2))+B6(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*(Hz(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)-Hz(pis:pie-2,pjs+1:pje-1,pks+1:pke-1));
    Ez(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)=C1(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*Ez(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)+C2(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*(Hy(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)-Hy(pis:pie-2,pjs+1:pje-1,pks+1:pke-1))...
        +C3(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1).*(Hx(pis+1:pie-1,pjs+1:pje-1,pks+1:pke-1)-Hx(pis+1:pie-1,pjs:pje-2,pks+1:pke-1));
    
    PML_Update_E;
    
    Hx(pis:pie-1,pjs:pje-1,pks:pke-1)=D1(pis:pie-1,pjs:pje-1,pks:pke-1).*Hx(pis:pie-1,pjs:pje-1,pks:pke-1)+D2(pis:pie-1,pjs:pje-1,pks:pke-1).*Hy(pis:pie-1,pjs:pje-1,pks:pke-1)...
        +D3(pis:pie-1,pjs:pje-1,pks:pke-1).*(Ez(pis:pie-1,pjs+1:pje,pks:pke-1)-Ez(pis:pie-1,pjs:pje-1,pks:pke-1))+D4(pis:pie-1,pjs:pje-1,pks:pke-1).*(Ey(pis:pie-1,pjs:pje-1,pks+1:pke)-Ey(pis:pie-1,pjs:pje-1,pks:pke-1))...
        +D5(pis:pie-1,pjs:pje-1,pks:pke-1).*(Ex(pis:pie-1,pjs:pje-1,pks+1:pke)-Ex(pis:pie-1,pjs:pje-1,pks:pke-1))+D6(pis:pie-1,pjs:pje-1,pks:pke-1).*(Ez(pis+1:pie,pjs:pje-1,pks:pke-1)-Ez(pis:pie-1,pjs:pje-1,pks:pke-1));
    Hy(pis:pie-1,pjs:pje-1,pks:pke-1)=E1(pis:pie-1,pjs:pje-1,pks:pke-1).*Hy(pis:pie-1,pjs:pje-1,pks:pke-1)+E2(pis:pie-1,pjs:pje-1,pks:pke-1).*Hx(pis:pie-1,pjs:pje-1,pks:pke-1)...
        +E3(pis:pie-1,pjs:pje-1,pks:pke-1).*(Ez(pis:pie-1,pjs+1:pje,pks:pke-1)-Ez(pis:pie-1,pjs:pje-1,pks:pke-1))+E4(pis:pie-1,pjs:pje-1,pks:pke-1).*(Ey(pis:pie-1,pjs:pje-1,pks+1:pke)-Ey(pis:pie-1,pjs:pje-1,pks:pke-1))...
        +E5(pis:pie-1,pjs:pje-1,pks:pke-1).*(Ex(pis:pie-1,pjs:pje-1,pks+1:pke)-Ex(pis:pie-1,pjs:pje-1,pks:pke-1))+E6(pis:pie-1,pjs:pje-1,pks:pke-1).*(Ez(pis+1:pie,pjs:pje-1,pks:pke-1)-Ez(pis:pie-1,pjs:pje-1,pks:pke-1));
    Hz(pis:pie-1,pjs:pje-1,pks:pke-1)=F1(pis:pie-1,pjs:pje-1,pks:pke-1).*Hz(pis:pie-1,pjs:pje-1,pks:pke-1)+F2(pis:pie-1,pjs:pje-1,pks:pke-1).*(Ex(pis:pie-1,pjs+1:pje,pks:pke-1)-Ex(pis:pie-1,pjs:pje-1,pks:pke-1))...
        +F3(pis:pie-1,pjs:pje-1,pks:pke-1).*(Ey(pis+1:pie,pjs:pje-1,pks:pke-1)-Ey(pis:pie-1,pjs:pje-1,pks:pke-1));

    PML_Update;
    
%     AA=Ex(:,:,dst_layer);
%     AA1=Exy_zn(:,:,dst_layer);
%     AA2=Exz_zn(:,:,dst_layer);
%     BB=Ey(:,:,dst_layer);
%     BB1=Eyx_zn(:,:,dst_layer);
%     BB2=Eyz_zn(:,:,dst_layer);
    
%     Hx(1:nx-1,1:ny-1,1:nz-1)=D1(1:nx-1,1:ny-1,1:nz-1).*Hx(1:nx-1,1:ny-1,1:nz-1)+D2(1:nx-1,1:ny-1,1:nz-1).*Hy(1:nx-1,1:ny-1,1:nz-1)...
%         +D3(1:nx-1,1:ny-1,1:nz-1).*(Ez(1:nx-1,2:ny,1:nz-1)-Ez(1:nx-1,1:ny-1,1:nz-1))+D4(1:nx-1,1:ny-1,1:nz-1).*(Ey(1:nx-1,1:ny-1,2:nz)-Ey(1:nx-1,1:ny-1,1:nz-1))...
%         +D5(1:nx-1,1:ny-1,1:nz-1).*(Ex(1:nx-1,1:ny-1,2:nz)-Ex(1:nx-1,1:ny-1,1:nz-1))+D6(1:nx-1,1:ny-1,1:nz-1).*(Ez(2:nx,1:ny-1,1:nz-1)-Ez(1:nx-1,1:ny-1,1:nz-1));
%     Hy(1:nx-1,1:ny-1,1:nz-1)=E1(1:nx-1,1:ny-1,1:nz-1).*Hy(1:nx-1,1:ny-1,1:nz-1)+E2(1:nx-1,1:ny-1,1:nz-1).*Hx(1:nx-1,1:ny-1,1:nz-1)...
%         +E3(1:nx-1,1:ny-1,1:nz-1).*(Ez(1:nx-1,2:ny,1:nz-1)-Ez(1:nx-1,1:ny-1,1:nz-1))+E4(1:nx-1,1:ny-1,1:nz-1).*(Ey(1:nx-1,1:ny-1,2:nz)-Ey(1:nx-1,1:ny-1,1:nz-1))...
%         +E5(1:nx-1,1:ny-1,1:nz-1).*(Ex(1:nx-1,1:ny-1,2:nz)-Ex(1:nx-1,1:ny-1,1:nz-1))+E6(1:nx-1,1:ny-1,1:nz-1).*(Ez(2:nx,1:ny-1,1:nz-1)-Ez(1:nx-1,1:ny-1,1:nz-1));
%     Hz(1:nx-1,1:ny-1,1:nz-1)=F1(1:nx-1,1:ny-1,1:nz-1).*Hz(1:nx-1,1:ny-1,1:nz-1)+F2(1:nx-1,1:ny-1,1:nz-1).*(Ex(1:nx-1,2:ny,1:nz-1)-Ex(1:nx-1,1:ny-1,1:nz-1))...
%         +F3(1:nx-1,1:ny-1,1:nz-1).*(Ey(2:nx,1:ny-1,1:nz-1)-Ey(1:nx-1,1:ny-1,1:nz-1));
%     
%     Ex(2:nx,2:ny,2:nz)=A1(2:nx,2:ny,2:nz).*Ex(2:nx,2:ny,2:nz)+A2(2:nx,2:ny,2:nz).*(Hz(2:nx,2:ny,2:nz)-Hz(2:nx,1:ny-1,2:nz))...
%         +A3(2:nx,2:ny,2:nz).*(Hy(2:nx,2:ny,2:nz)-Hy(2:nx,2:ny,1:nz-1));
%     Ey(2:nx,2:ny,2:nz)=B1(2:nx,2:ny,2:nz).*Ey(2:nx,2:ny,2:nz)+B2(2:nx,2:ny,2:nz).*(Hx(2:nx,2:ny,2:nz)-Hx(2:nx,2:ny,1:nz-1))...
%         +B3(2:nx,2:ny,2:nz).*(Hz(2:nx,2:ny,2:nz)-Hz(1:nx-1,2:ny,2:nz));
%     Ez(2:nx,2:ny,2:nz)=C1(2:nx,2:ny,2:nz).*Ez(2:nx,2:ny,2:nz)+C2(2:nx,2:ny,2:nz).*(Hy(2:nx,2:ny,2:nz)-Hy(1:nx-1,2:ny,2:nz))...
%         +C3(2:nx,2:ny,2:nz).*(Hx(2:nx,2:ny,2:nz)-Hx(2:nx,1:ny-1,2:nz));
    
    
    
%     temp=real(Ez(:,:,dst_layer));  % 计算波形的符号
%     for i=1:nx
%         for j=1:ny
%             if temp(i,j)~=0
%                 temp(i,j)=temp(i,j)/abs(temp(i,j));
%             end
%         end
%     end
    
    for i=source_layer:dst_layer
        temp0=real(Ex(source_mid_x,source_mid_y,i));
%         temp0=real(Ex(source_mid_x,source_mid_y,dst_layer));
%         if temp0~=0
%             temp0=temp0/abs(temp0);
%         end
%         Jx(n)=temp0*sqrt(abs(Ex(source_mid_x,source_mid_y,dst_layer)));
        Jx(i-source_layer+1,n)=temp0;
%         Jx(n)=temp0;
        temp1=real(Ey(source_mid_x,source_mid_y,i));
%         temp1=real(Ey(source_mid_x,source_mid_y,dst_layer));
%         if temp1~=0
%             temp1=temp1/abs(temp1);
%         end
%         Jy(n)=temp1*sqrt(abs(Ey(source_mid_x,source_mid_y,dst_layer)));
        Jy(i-source_layer+1,n)=temp1;
%         Jy(n)=temp1;
    end
    
%     if n>70
%         AA=real(Ex(:,:,dst_layer));
%     end
    
    Sumx(:,:)=real(Ex(:,1:ny,dst_layer));
%     Sumy(:,:)=temp.*sqrt(Ey(:,:,dst_layer).*conj(Ey(:,:,dst_layer)));
%     Sumz(:,:)=temp.*sqrt(Ez(:,:,dst_layer).*conj(Ez(:,:,dst_layer)));
    Sumy(:,:)=real(Ey(1:nx,:,dst_layer));
    Sumz(:,:)=real(Ez(1:nx,1:ny,dst_layer));
    Sumhx(:,:)=real(Hx(1:nx,:,dst_layer));
    Sumhy(:,:)=real(Hy(:,1:ny,dst_layer));
    Sumhz(:,:)=real(Hz(:,:,dst_layer));
    
%     if(n>1000)
%         Sumx_t(:,:,n-1000)=Sumx(:,:);
%         Sumy_t(:,:,n-1000)=Sumy(:,:);
%         Sumz_t(:,:,n-1000)=Sumz(:,:);
%         Sumhx_t(:,:,n-1000)=Sumhx(:,:);
%         Sumhy_t(:,:,n-1000)=Sumhy(:,:);
%         Sumhz_t(:,:,n-1000)=Sumhz(:,:);
%     end

%     if(n>2000)
%         h=imagesc(1e+6*dy*(1:1:ny),(dx*1e+6*(1:1:nx))',Sumx');%,[-0.05,0.05]);
%         colorbar;
%         set(h,'AlphaData',epsilon(:,:,dst_layer)'/epsilon0);
% %         Sum1(:,:,n-3000)=Sumx(:,:);
% %         Sum2(:,:,n-3000)=Sumy(:,:);
%         shading interp;
%         getframe;
%     end

    disp([' n = ' num2str(n)]);
end

for i=source_layer:dst_layer
    figure;
    plot(Jx(i-source_layer+1,nt-1000:nt)',Jy(i-source_layer+1,nt-1000:nt)');
%     plot(Jx(3000:3600),Jy(3000:3600));
    title('Ex--Ey偏振图像');
    xlabel('Ex(N/C)');
    ylabel('Ey(N/C)');
    axis equal;
%     axis([-0.1 0.1 -0.1 0.1]);
    pause(1);
end