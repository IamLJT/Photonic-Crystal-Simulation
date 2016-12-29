% Ex
Exy_xn=Cexye_xn .* Exy_xn + Cexyhy_xn .* (Hz(2:pis,pjs+1:pje-1,pks+1:pke-1)-Hz(2:pis,pjs:pje-2,pks+1:pke-1));
Exy_xp=Cexye_xp .* Exy_xp + Cexyhy_xp .* (Hz(pie:nx,pjs+1:pje-1,pks+1:pke-1)-Hz(pie:nx,pjs:pje-2,pks+1:pke-1));
Exy_yn=Cexye_yn .* Exy_yn + Cexyhy_yn .* (Hz(2:nx,2:pjs,2:nz)-Hz(2:nx,1:pjs-1,2:nz));
Exy_yp=Cexye_yp .* Exy_yp + Cexyhy_yp .* (Hz(2:nx,pje:ny,2:nz)-Hz(2:nx,pje-1:ny-1,2:nz));
Exy_zn=Cexye_zn .* Exy_zn + Cexyhy_zn .* (Hz(2:nx,pjs+1:pje-1,2:pks)-Hz(2:nx,pjs:pje-2,2:pks));
Exy_zp=Cexye_zp .* Exy_zp + Cexyhy_zp .* (Hz(2:nx,pjs+1:pje-1,pke:nz)-Hz(2:nx,pjs:pje-2,pke:nz));

Exz_xn=Cexze_xn .* Exz_xn + Cexzhy_xn .* (Hy(2:pis,pjs+1:pje-1,pks+1:pke-1)-Hy(2:pis,pjs+1:pje-1,pks:pke-2));
Exz_xp=Cexze_xp .* Exz_xp + Cexzhy_xp .* (Hy(pie:nx,pjs+1:pje-1,pks+1:pke-1)-Hy(pie:nx,pjs+1:pje-1,pks:pke-2));
Exz_yn=Cexze_yn .* Exz_yn + Cexzhy_yn .* (Hy(2:nx,2:pjs,pks+1:pke-1)-Hy(2:nx,2:pjs,pks:pke-2));
Exz_yp=Cexze_yp .* Exz_yp + Cexzhy_yp .* (Hy(2:nx,pje:ny,pks+1:pke-1)-Hy(2:nx,pje:ny,pks:pke-2));
Exz_zn=Cexze_zn .* Exz_zn + Cexzhy_zn .* (Hy(2:nx,2:ny,2:pks)-Hy(2:nx,2:ny,1:pks-1));
Exz_zp=Cexze_zp .* Exz_zp + Cexzhy_zp .* (Hy(2:nx,2:ny,pke:nz)-Hy(2:nx,2:ny,pke-1:nz-1));

% Ey
Eyx_xn=Ceyxe_xn .* Eyx_xn + Ceyxhy_xn .* (Hz(2:pis,2:ny,2:nz)-Hz(1:pis-1,2:ny,2:nz));
Eyx_xp=Ceyxe_xp .* Eyx_xp + Ceyxhy_xp .* (Hz(pie:nx,2:ny,2:nz)-Hz(pie-1:nx-1,2:ny,2:nz));
Eyx_yn=Ceyxe_yn .* Eyx_yn + Ceyxhy_yn .* (Hz(pis+1:pie-1,2:pjs,pks+1:pke-1)-Hz(pis:pie-2,2:pjs,pks+1:pke-1));
Eyx_yp=Ceyxe_yp .* Eyx_yp + Ceyxhy_yp .* (Hz(pis+1:pie-1,pje:ny,pks+1:pke-1)-Hz(pis:pie-2,pje:ny,pks+1:pke-1));
Eyx_zn=Ceyxe_zn .* Eyx_zn + Ceyxhy_zn .* (Hz(pis+1:pie-1,2:ny,2:pks)-Hz(pis:pie-2,2:ny,2:pks));
Eyx_zp=Ceyxe_zp .* Eyx_zp + Ceyxhy_zp .* (Hz(pis+1:pie-1,2:ny,pke:nz)-Hz(pis:pie-2,2:ny,pke:nz));

Eyz_xn=Ceyze_xn .* Eyz_xn + Ceyzhy_xn .* (Hx(2:pis,2:ny,pks+1:pke-1)-Hx(2:pis,2:ny,pks:pke-2));
Eyz_xp=Ceyze_xp .* Eyz_xp + Ceyzhy_xp .* (Hx(pie:nx,2:ny,pks+1:pke-1)-Hx(pie:nx,2:ny,pks:pke-2));
Eyz_yn=Ceyze_yn .* Eyz_yn + Ceyzhy_yn .* (Hx(pis+1:pje-1,2:pjs,pks+1:pke-1)-Hx(pis+1:pie-1,2:pjs,pks:pke-2));
Eyz_yp=Ceyze_yp .* Eyz_yp + Ceyzhy_yp .* (Hx(pis+1:pje-1,pje:ny,pks+1:pke-1)-Hx(pis+1:pie-1,pje:ny,pks:pke-2));
Eyz_zn=Ceyze_zn .* Eyz_zn + Ceyzhy_zn .* (Hx(2:nx,2:ny,2:pks)-Hx(2:nx,2:ny,1:pks-1));
Eyz_zp=Ceyze_zp .* Eyz_zp + Ceyzhy_zp .* (Hx(2:nx,2:ny,pke:nz)-Hx(2:nx,2:ny,pke-1:nz-1));

% Ez
Ezx_xn=Cezxe_xn .* Ezx_xn + Cezxhy_xn .* (Hy(2:pis,2:ny,2:nz)-Hy(1:pis-1,2:ny,2:nz));
Ezx_xp=Cezxe_xp .* Ezx_xp + Cezxhy_xp .* (Hy(pie:nx,2:ny,2:nz)-Hy(pie-1:nx-1,2:ny,2:nz));
Ezx_yn=Cezxe_yn .* Ezx_yn + Cezxhy_yn .* (Hy(pis+1:pie-1,2:pjs,2:nz)-Hy(pis:pie-2,2:pjs,2:nz));
Ezx_yp=Cezxe_yp .* Ezx_yp + Cezxhy_yp .* (Hy(pis+1:pie-1,pje:ny,2:nz)-Hy(pis:pie-2,pje:ny,2:nz));
Ezx_zn=Cezxe_zn .* Ezx_zn + Cezxhy_zn .* (Hy(pis+1:pie-1,pjs+1:pje-1,2:pks)-Hy(pis:pie-2,pjs+1:pje-1,2:pks));
Ezx_zp=Cezxe_zp .* Ezx_zp + Cezxhy_zp .* (Hy(pis+1:pie-1,pjs+1:pje-1,pke:nz)-Hy(pis:pie-2,pjs+1:pje-1,pke:nz));

Ezy_yn=Cezye_yn .* Ezy_yn + Cezyhy_yn .* (Hx(2:nx,2:pjs,2:nz)-Hx(2:nx,1:pjs-1,2:nz));
Ezy_yp=Cezye_yp .* Ezy_yp + Cezyhy_yp .* (Hx(2:nx,pje:ny,2:nz)-Hx(2:nx,pje-1:ny-1,2:nz));
Ezy_xn=Cezye_xn .* Ezy_xn + Cezyhy_xn .* (Hx(2:pis,pjs+1:pje-1,2:nz)-Hx(2:pis,pjs:pje-2,2:nz));
Ezy_xp=Cezye_xp .* Ezy_xp + Cezyhy_xp .* (Hx(pie:nx,pjs+1:pje-1,2:nz)-Hx(pie:nx,pjs:pje-2,2:nz));
Ezy_zn=Cezye_zn .* Ezy_zn + Cezyhy_zn .* (Hx(pis+1:pie-1,pjs+1:pje-1,2:pks)-Hx(pis+1:pie-1,pjs:pje-2,2:pks));
Ezy_zp=Cezye_zp .* Ezy_zp + Cezyhy_zp .* (Hx(pis+1:pie-1,pjs+1:pje-1,pke:nz)-Hx(pis+1:pie-1,pjs:pje-2,pke:nz));

Ex(2:nx,2:pis,2:pks)=Exy_yn(:,:,1:pks-1)+Exz_zn(:,1:pjs-1,:);
Ex(2:nx,pje:ny,2:pks)=Exy_yp(:,:,1:pks-1)+Exz_zn(:,pje-1:ny-1,:);
Ex(2:nx,2:pis,pke:nz)=Exy_yn(:,:,pke-1:nz-1)+Exz_zp(:,1:pis-1,:);
Ex(2:nx,pje:ny,pke:nz)=Exy_yp(:,:,pke-1:nz-1)+Exz_zp(:,pje-1:ny-1,:);
Ex(2:nx,2:pjs,pks+1:pke-1)=Exy_yn(:,:,pks:pke-2)+Exz_yn;
Ex(2:nx,pje:ny,pks+1:pke-1)=Exy_yp(:,:,pks:pke-2)+Exz_yp;
Ex(2:nx,pjs+1:pje-1,2:pks)=Exy_zn+Exz_zn(:,pjs:pje-2,:);
Ex(2:nx,pjs+1:pje-1,pke:nz)=Exy_zp+Exz_zp(:,pjs:pje-2,:);
Ex(2:pis,pjs+1:pje-1,pks+1:pke-1)=Exy_xn+Exz_xn;
Ex(pie:nx,pjs+1:pje-1,pks+1:pke-1)=Exy_xp+Exz_xp;

Ey(2:pis,2:ny,2:pks)=Eyx_xn(:,:,1:pks-1)+Eyz_zn(1:pis-1,:,:);
% Ey(2:pis,:,2:pks)=Eyx_xn(:,:,1:pks-1)+Eyz_zn(1:pis-1,:,:);
Ey(2:pis,2:ny,pke:nz)=Eyx_xn(:,:,pke-1:nz-1)+Eyz_zp(1:pis-1,:,:);
% Ey(2:pis,:,pke+1:nz)=Eyx_xn(:,:,pke+1:nz)+Eyz_yp(1:pis-1,:,:);
Ey(pie:nx,2:ny,2:pks)=Eyx_xp(:,:,1:pks-1)+Eyz_zn(pie-1:nx-1,:,:);
Ey(pie:nx,2:ny,pke:nz)=Eyx_xp(:,:,pke-1:nz-1)+Eyz_zp(pie-1:nx-1,:,:);
% Ey(pie+1:nx,:,2:pks)=Eyx_xp(:,:,1:pks-1)+Eyz_yp(pie:nx-1,:,:);
% Ey(pie+1:nx,:,pke+1:nz)=Eyx_xp(:,:,pke+1:nz)+Eyz_yp(pie:nx-1,:,:);
Ey(2:pis,2:ny,pks+1:pke-1)=Eyx_xn(:,:,pks:pke-2)+Eyz_xn;
Ey(pie:nx,2:ny,pks+1:pke-1)=Eyx_xp(:,:,pks:pke-2)+Eyz_xp;
Ey(pis+1:pie-1,2:ny,2:pks)=Eyx_zn+Eyz_zn(pis:pie-2,:,:);
Ey(pis+1:pie-1,2:ny,pke:nz)=Eyx_zp+Eyz_zp(pis:pie-2,:,:);
Ey(pis+1:pie-1,2:pjs,pks+1:pke-1)=Eyx_yn+Eyz_yn;
Ey(pis+1:pie-1,pje:ny,pks+1:pke-1)=Eyx_yp+Eyz_yp;

Ez(2:pis,2:pjs,2:nz)=Ezx_xn(:,1:pjs-1,:)+Ezy_yn(1:pis-1,:,:);
% Ez(2:pis,2:pjs,:)=Ezx_xn(:,1:pjs-1,:)+Ezy_yn(1:pis-1,:,:);
Ez(2:pis,pje:ny,2:nz)=Ezx_xn(:,pje-1:ny-1,:)+Ezy_yp(1:pis-1,:,:);
% Ez(2:pis,pje+1:ny,:)=Ezx_xn(:,pje+1:ny,:)+Ezy_yp(1:pis-1,:,:);
Ez(pie:nx,2:pjs,2:nz)=Ezx_xp(:,1:pjs-1,:)+Ezy_yn(pie-1:nx-1,:,:);
% Ez(pie+1:nx,2:pjs,:)=Ezx_xp(:,1:pjs-1,:)+Ezy_yn(pie:nx-1,:,:);
Ez(pie:nx,pje:ny,2:nz)=Ezx_xp(:,pje-1:ny-1,:)+Ezy_yp(pie-1:nx-1,:,:);
% Ez(pie+1:nx,pje+1:ny,:)=Ezx_xp(:,pje+1:ny,:)+Ezy_yp(pie:nx-1,:,:);
Ez(2:pis,pjs+1:pje-1,2:nz)=Ezx_xn(:,pjs:pje-2,:)+Ezy_xn;
Ez(pie:nx,pjs+1:pje-1,2:nz)=Ezx_xp(:,pjs:pje-2,:)+Ezy_xp;
Ez(pis+1:pie-1,2:pjs,2:nz)=Ezx_yn+Ezy_yn(pis:pje-2,:,:);
Ez(pis+1:pie-1,pie:ny,2:nz)=Ezx_yp+Ezy_yp(pis:pje-2,:,:);
Ez(pis+1:pie-1,pjs+1:pje-1,2:pks)=Ezx_zn+Ezy_zn;
Ez(pis+1:pie-1,pjs+1:pje-1,pke:nz)=Ezx_zp+Ezy_zp;