%---------更新方程---------------
% Hz
Hzx_xn=Chzxh_xn .* Hzx_xn + Chzxey_xn .* (Ey(2:pis,1:ny,1:nz)-Ey(1:pis-1,1:ny,1:nz));
Hzx_xp=Chzxh_xp .* Hzx_xp + Chzxey_xp .* (Ey(pie+1:nxp1,1:ny,1:nz)-Ey(pie:nx,1:ny,1:nz));
Hzx_yn=Chzxh_yn .* Hzx_yn + Chzxey_yn .* (Ey(pis+1:pie,1:pjs-1,1:nz)-Ey(pis:pie-1,1:pjs-1,1:nz));
Hzx_yp=Chzxh_yp .* Hzx_yp + Chzxey_yp .* (Ey(pis+1:pie,pje:ny,1:nz)-Ey(pis:pie-1,pje:ny,1:nz));
Hzx_zn=Chzxh_zn .* Hzx_zn + Chzxey_zn .* (Ey(pis+1:pie,pjs:pje-1,1:pks-1)-Ey(pis:pie-1,pjs:pje-1,1:pks-1));
Hzx_zp=Chzxh_zp .* Hzx_zp + Chzxey_zp .* (Ey(pis+1:pie,pjs:pje-1,pke:nz)-Ey(pis:pie-1,pjs:pje-1,pke:nz));

Hzy_yn=Chzyh_yn .* Hzy_yn + Chzyey_yn .* (Ex(:,2:pjs,1:nz)-Ex(:,1:pjs-1,1:nz));
Hzy_yp=Chzyh_yp .* Hzy_yp + Chzyey_yp .* (Ex(:,pje+1:nyp1,1:nz)-Ex(:,pje:ny,1:nz));
Hzy_xn=Chzyh_xn .* Hzy_xn + Chzyey_xn .* (Ex(1:pis-1,pjs+1:pje,1:nz)-Ex(1:pis-1,pjs:pje-1,1:nz));
Hzy_xp=Chzyh_xp .* Hzy_xp + Chzyey_xp .* (Ex(pie:nx,pjs+1:pje,1:nz)-Ex(pie:nx,pjs:pje-1,1:nz));
Hzy_zn=Chzyh_zn .* Hzy_zn + Chzyey_zn .* (Ex(pis:pie-1,pjs+1:pje,1:pks-1)-Ex(pis:pie-1,pjs:pje-1,1:pks-1));
Hzy_zp=Chzyh_zp .* Hzy_zp + Chzyey_zp .* (Ex(pis:pie-1,pjs+1:pje,pke:nz)-Ex(pis:pie-1,pjs:pje-1,pke:nz));



% Hy
Hyx_xn=Chyxh_xn .* Hyx_xn + Chyxey_xn .* (Ez(2:pis,1:ny,:)-Ez(1:pis-1,1:ny,:));
Hyx_xp=Chyxh_xp .* Hyx_xp + Chyxey_xp .* (Ez(pie+1:nxp1,1:ny,:)-Ez(pie:nx,1:ny,:));
Hyx_yn=Chyxh_yn .* Hyx_yn + Chyxey_yn .* (Ez(pis+1:pie,1:pjs-1,pks:pke-1)-Ez(pis:pie-1,1:pjs-1,pks:pke-1));
Hyx_yp=Chyxh_yp .* Hyx_yp + Chyxey_yp .* (Ez(pis+1:pie,pje:ny,pks:pke-1)-Ez(pis:pie-1,pje:ny,pks:pke-1));
Hyx_zn=Chyxh_zn .* Hyx_zn + Chyxey_zn .* (Ez(pis+1:pie,1:ny,1:pks-1)-Ez(pis:pie-1,1:ny,1:pks-1));
Hyx_zp=Chyxh_zp .* Hyx_zp + Chyxey_zp .* (Ez(pis+1:pie,1:ny,pke:nz)-Ez(pis:pie-1,1:ny,pke:nz));

Hyz_xn=Chyzh_xn .* Hyz_xn + Chyzey_xn .* (Ex(1:pis-1,1:ny,pks+1:pke)-Ex(1:pis-1,1:ny,pks:pke-1));
Hyz_xp=Chyzh_xp .* Hyz_xp + Chyzey_xp .* (Ex(pie:nx,1:ny,pks+1:pke)-Ex(pie:nx,1:ny,pks:pke-1));
Hyz_yn=Chyzh_yn .* Hyz_yn + Chyzey_yn .* (Ex(pis:pje-1,1:pjs-1,pks+1:pke)-Ex(pis:pie-1,1:pjs-1,pks:pke-1));
Hyz_yp=Chyzh_yp .* Hyz_yp + Chyzey_yp .* (Ex(pis:pje-1,pje:ny,pks+1:pke)-Ex(pis:pie-1,pje:ny,pks:pke-1));
Hyz_zn=Chyzh_zn .* Hyz_zn + Chyzey_zn .* (Ex(:,1:ny,2:pks)-Ex(:,1:ny,1:pks-1));
Hyz_zp=Chyzh_zp .* Hyz_zp + Chyzey_zp .* (Ex(:,1:ny,pke+1:nzp1)-Ex(:,1:ny,pke:nz));



% Hx
Hxy_xn=Chxyh_xn .* Hxy_xn + Chxyey_xn .* (Ez(1:pis-1,pjs+1:pje,pks:pke-1)-Ez(1:pis-1,pjs:pje-1,pks:pke-1));
Hxy_xp=Chxyh_xp .* Hxy_xp + Chxyey_xp .* (Ez(pie:nx,pjs+1:pje,pks:pke-1)-Ez(pie:nx,pjs:pje-1,pks:pke-1));
Hxy_yn=Chxyh_yn .* Hxy_yn + Chxyey_yn .* (Ez(1:nx,2:pjs,:)-Ez(1:nx,1:pjs-1,:));
Hxy_yp=Chxyh_yp .* Hxy_yp + Chxyey_yp .* (Ez(1:nx,pje+1:nyp1,:)-Ez(1:nx,pje:ny,:));
Hxy_zn=Chxyh_zn .* Hxy_zn + Chxyey_zn .* (Ez(1:nx,pjs+1:pje,1:pks-1)-Ez(1:nx,pjs:pje-1,1:pks-1));
Hxy_zp=Chxyh_zp .* Hxy_zp + Chxyey_zp .* (Ez(1:nx,pjs+1:pje,pke:nz)-Ez(1:nx,pjs:pje-1,pke:nz));

Hxz_xn=Chxzh_xn .* Hxz_xn + Chxzey_xn .* (Ey(1:pis-1,pjs:pje-1,pks+1:pke)-Ey(1:pis-1,pjs:pje-1,pks:pke-1));
Hxz_xp=Chxzh_xp .* Hxz_xp + Chxzey_xp .* (Ey(pie:nx,pjs:pje-1,pks+1:pke)-Ey(pie:nx,pjs:pje-1,pks:pke-1));
Hxz_yn=Chxzh_yn .* Hxz_yn + Chxzey_yn .* (Ey(1:nx,1:pjs-1,pks+1:pke)-Ey(1:nx,1:pjs-1,pks:pke-1));
Hxz_yp=Chxzh_yp .* Hxz_yp + Chxzey_yp .* (Ey(1:nx,pje:ny,pks+1:pke)-Ey(1:nx,pje:ny,pks:pke-1));
Hxz_zn=Chxzh_zn .* Hxz_zn + Chxzey_zn .* (Ey(1:nx,:,2:pks)-Ey(1:nx,:,1:pks-1));
Hxz_zp=Chxzh_zp .* Hxz_zp + Chxzey_zp .* (Ey(1:nx,:,pke+1:nzp1)-Ey(1:nx,:,pke:nz));



Hz(1:pis-1,1:pjs-1,1:nz)=Hzx_xn(:,1:pjs-1,:)+Hzy_yn(1:pis-1,:,:);
% Hz(2:pis,2:pjs,:)=Hzx_xn(:,1:pjs-1,:)+Hzy_yn(1:pis-1,:,:);
Hz(1:pis-1,pje:ny,1:nz)=Hzx_xn(:,pje:ny,:)+Hzy_yp(1:pis-1,:,:);
% Hz(2:pis,pje:ny,:)=Hzx_xn(:,pje:ny,:)+Hzy_yp(1:pis-1,:,:);
Hz(pie:nx,1:pjs-1,1:nz)=Hzx_xp(:,1:pjs-1,:)+Hzy_yn(pie:nx,:,:);
% Hz(pie+1:nx,2:pjs,:)=Hzx_xp(:,1:pjs-1,:)+Hzy_yn(pie+1:nxp1-1,:,:);
Hz(pie:nx,pje:ny,1:nz)=Hzx_xp(:,pje:ny,:)+Hzy_yp(pie:nx,:,:);
% Hz(pie+1:nx,pje:ny,:)=Hzx_xp(:,pje:ny,:)+Hzy_yp(pie:nx-1,:,:);
Hz(1:pis-1,pjs:pje-1,1:nz)=Hzx_xn(:,pjs:pje-1,:)+Hzy_xn;
Hz(pie:nx,pjs:pje-1,1:nz)=Hzx_xp(:,pjs:pje-1,:)+Hzy_xp;
Hz(pis:pie-1,1:pjs-1,1:nz)=Hzx_yn+Hzy_yn(pis:pie-1,:,:);
Hz(pis:pie-1,pie:ny,1:nz)=Hzx_yp+Hzy_yp(pis:pje-1,:,:);
Hz(pis:pie-1,pjs:pje-1,1:pks-1)=Hzx_zn+Hzy_zn;
Hz(pis:pie-1,pjs:pje-1,pke:nz)=Hzx_zp+Hzy_zp;


Hy(1:pis-1,1:ny,1:pks-1)=Hyx_xn(:,:,1:pks-1)+Hyz_zn(1:pis-1,:,:);
Hy(1:pis-1,1:ny,pke:nz)=Hyx_xn(:,:,pke:nz)+Hyz_zp(1:pis-1,:,:);
Hy(pie:nx,1:ny,1:pks-1)=Hyx_xp(:,:,1:pks-1)+Hyz_zn(pie:nx,:,:);
Hy(pie:nx,1:ny,pke:nz)=Hyx_xp(:,:,pke:nz)+Hyz_zp(pie:nx,:,:);
Hy(1:pis-1,1:ny,pks:pke-1)=Hyx_xn(:,:,pks:pke-1)+Hyz_xn;
Hy(pie:nx,1:ny,pks:pke-1)=Hyx_xp(:,:,pks:pke-1)+Hyz_xp;
Hy(pis:pie-1,1:ny,1:pks-1)=Hyx_zn+Hyz_zn(pis:pie-1,:,:);
Hy(pis:pie-1,1:ny,pke:nz)=Hyx_zp+Hyz_zp(pis:pie-1,:,:);
Hy(pis:pie-1,1:pjs-1,pks:pke-1)=Hyx_yn+Hyz_yn;
Hy(pis:pie-1,pje:ny,pks:pke-1)=Hyx_yp+Hyz_yp;



Hx(1:nx,1:pis-1,1:pks-1)=Hxy_yn(:,:,1:pks-1)+Hxz_zn(:,1:pjs-1,:);
Hx(1:nx,pje:ny,1:pks-1)=Hxy_yp(:,:,1:pks-1)+Hxz_zn(:,pje:ny,:);
Hx(1:nx,1:pis-1,pke:nz)=Hxy_yn(:,:,pke:nz)+Hxz_zp(:,1:pjs-1,:);
Hx(1:nx,pje:ny,pke:nz)=Hxy_yp(:,:,pke:nz)+Hxz_zp(:,pje:ny,:);
Hx(1:nx,1:pjs-1,pks:pke-1)=Hxy_yn(:,:,pks:pke-1)+Hxz_yn;
Hx(1:nx,pje:ny,pks:pke-1)=Hxy_yp(:,:,pks:pke-1)+Hxz_yp;
Hx(1:nx,pjs:pje-1,1:pks-1)=Hxy_zn+Hxz_zn(:,pjs:pje-1,:);
Hx(1:nx,pjs:pje-1,pke:nz)=Hxy_zp+Hxz_zp(:,pjs:pje-1,:);
Hx(1:pis-1,pjs:pje-1,pks:pke-1)=Hxy_xn+Hxz_xn;
Hx(pie:nx,pjs:pje-1,pks:pke-1)=Hxy_xp+Hxz_xp;