if Lattice_Style==1
    eps_0=epsilon0*epsr_material;
elseif Lattice_Style==2
    eps_0=epsilon0;
end

pis = Npml+1;
pie = nx-Npml+1;
pjs = Npml+1;
pje = ny-Npml+1;
pks = Npml+1;
pke = nz-Npml+1;

Ezx_xn=zeros(Npml,nym1,nzm1);
Ezx_xp=zeros(Npml,nym1,nzm1);
Ezx_yn=zeros(nxm1-2*Npml,Npml,nzm1);
Ezx_yp=zeros(nxm1-2*Npml,Npml,nzm1);
Ezx_zn=zeros(nxm1-2*Npml,nym1-2*Npml,Npml);
Ezx_zp=zeros(nxm1-2*Npml,nym1-2*Npml,Npml);

Ezy_xn=zeros(Npml,nym1-2*Npml,nzm1);
Ezy_xp=zeros(Npml,nym1-2*Npml,nzm1);
Ezy_yn=zeros(nxm1,Npml,nzm1);
Ezy_yp=zeros(nxm1,Npml,nzm1);
Ezy_zn=zeros(nxm1-2*Npml,nym1-2*Npml,Npml);
Ezy_zp=zeros(nxm1-2*Npml,nym1-2*Npml,Npml);

Exy_yn=zeros(nxm1,Npml,nzm1);
Exy_yp=zeros(nxm1,Npml,nzm1);
Exy_zn=zeros(nxm1,nym1-2*Npml,Npml);
Exy_zp=zeros(nxm1,nym1-2*Npml,Npml);
Exy_xn=zeros(Npml,nym1-2*Npml,nzm1-2*Npml);
Exy_xp=zeros(Npml,nym1-2*Npml,nzm1-2*Npml);

Exz_yn=zeros(nxm1,Npml,nzm1-2*Npml);
Exz_yp=zeros(nxm1,Npml,nzm1-2*Npml);
Exz_zn=zeros(nxm1,nym1,Npml);
Exz_zp=zeros(nxm1,nym1,Npml);
Exz_xn=zeros(Npml,nym1-2*Npml,nzm1-2*Npml);
Exz_xp=zeros(Npml,nym1-2*Npml,nzm1-2*Npml);

Eyx_xn=zeros(Npml,nym1,nzm1);
Eyx_xp=zeros(Npml,nym1,nzm1);
Eyx_zn=zeros(nxm1-2*Npml,nym1,Npml);
Eyx_zp=zeros(nxm1-2*Npml,nym1,Npml);
Eyx_yn=zeros(nxm1-2*Npml,Npml,nzm1-2*Npml);
Eyx_yp=zeros(nxm1-2*Npml,Npml,nzm1-2*Npml);

Eyz_xn=zeros(Npml,nym1,nzm1-2*Npml);
Eyz_xp=zeros(Npml,nym1,nzm1-2*Npml);
Eyz_zn=zeros(nxm1,nym1,Npml);
Eyz_zp=zeros(nxm1,nym1,Npml);
Eyz_yn=zeros(nxm1-2*Npml,Npml,nzm1-2*Npml);
Eyz_yp=zeros(nxm1-2*Npml,Npml,nzm1-2*Npml);

Hzx_xn=zeros(Npml,ny,nz);
Hzx_xp=zeros(Npml,ny,nz);
Hzx_yn=zeros(nx-2*Npml,Npml,nz);
Hzx_yp=zeros(nx-2*Npml,Npml,nz);
Hzx_zn=zeros(nx-2*Npml,ny-2*Npml,Npml);
Hzx_zp=zeros(nx-2*Npml,ny-2*Npml,Npml);

Hzy_xn=zeros(Npml,ny-2*Npml,nz);
Hzy_xp=zeros(Npml,ny-2*Npml,nz);
Hzy_yn=zeros(nx,Npml,nz);
Hzy_yp=zeros(nx,Npml,nz);
Hzy_zn=zeros(nx-2*Npml,ny-2*Npml,Npml);
Hzy_zp=zeros(nx-2*Npml,ny-2*Npml,Npml);

Hxy_yn=zeros(nx,Npml,nz);
Hxy_yp=zeros(nx,Npml,nz);
Hxy_zn=zeros(nx,ny-2*Npml,Npml);
Hxy_zp=zeros(nx,ny-2*Npml,Npml);
Hxy_xn=zeros(Npml,ny-2*Npml,nz-2*Npml);
Hxy_xp=zeros(Npml,ny-2*Npml,nz-2*Npml);

Hxz_yn=zeros(nx,Npml,nz-2*Npml);
Hxz_yp=zeros(nx,Npml,nz-2*Npml);
Hxz_zn=zeros(nx,ny,Npml);
Hxz_zp=zeros(nx,ny,Npml);
Hxz_xn=zeros(Npml,ny-2*Npml,nz-2*Npml);
Hxz_xp=zeros(Npml,ny-2*Npml,nz-2*Npml);

Hyx_xn=zeros(Npml,ny,nz);
Hyx_xp=zeros(Npml,ny,nz);
Hyx_zn=zeros(nx-2*Npml,ny,Npml);
Hyx_zp=zeros(nx-2*Npml,ny,Npml);
Hyx_yn=zeros(nx-2*Npml,Npml,nz-2*Npml);
Hyx_yp=zeros(nx-2*Npml,Npml,nz-2*Npml);

Hyz_xn=zeros(Npml,ny,nz-2*Npml);
Hyz_xp=zeros(Npml,ny,nz-2*Npml);
Hyz_zn=zeros(nx,ny,Npml);
Hyz_zp=zeros(nx,ny,Npml);
Hyz_yn=zeros(nx-2*Npml,Npml,nz-2*Npml);
Hyz_yp=zeros(nx-2*Npml,Npml,nz-2*Npml);

sigmax_xn=sigmax(1:Npml,1:nym1,1:nzm1);
sigmax_xp=sigmax(pie:nx,1:nym1,1:nzm1);
sigmay_yn=sigmay(1:nxm1,1:Npml,1:nzm1);
sigmay_yp=sigmay(1:nxm1,pje:ny,1:nzm1);
sigmaz_zn=sigmaz(1:nxm1,1:nym1,1:Npml);
sigmaz_zp=sigmaz(1:nxm1,1:nym1,pke:nz);

% Ezx、Ezy
% Ezx_xn
Cezxe_xn=(2*eps_0-dt*sigmax_xn)./(2*eps_0+dt*sigmax_xn);
Cezxhy_xn=(2*dt/dx)./(2*eps_0+dt*sigmax_xn);

% Ezx_xp
Cezxe_xp=(2*eps_0-dt*sigmax_xp)./(2*eps_0+dt*sigmax_xp);
Cezxhy_xp=(2*dt/dx)./(2*eps_0+dt*sigmax_xp);

% Ezx_yn
Cezxe_yn=1;
Cezxhy_yn=dt/(dx*eps_0);

% Ezx_yp
Cezxe_yp=1;
Cezxhy_yp=dt/(dx*eps_0);

% Ezx_zn
Cezxe_zn=1;
Cezxhy_zn=dt/(dx*eps_0);

% Ezx_zp
Cezxe_zp=1;
Cezxhy_zp=dt/(dx*eps_0);

% Ezy_xn
Cezye_xn=1;
Cezyhy_xn=-dt/(dy*eps_0);

% Ezy_xp
Cezye_xp=1;
Cezyhy_xp=-dt/(dy*eps_0);

% Ezy_yn
Cezye_yn=(2*eps_0-dt*sigmay_yn)./(2*eps_0+dt*sigmay_yn);
Cezyhy_yn=-(2*dt/dy)./(2*eps_0+dt*sigmay_yn);

% Ezy_yp
Cezye_yp=(2*eps_0-dt*sigmay_yp)./(2*eps_0+dt*sigmay_yp);
Cezyhy_yp=-(2*dt/dy)./(2*eps_0+dt*sigmay_yp);

% Ezy_zn
Cezye_zn=1;
Cezyhy_zn=-dt/(dy*eps_0);

% Ezy_zp
Cezye_zp=1;
Cezyhy_zp=-dt/(dy*eps_0);

% Exy、Exz
% Exy_xn
Cexye_xn=1;
Cexyhy_xn=dt/(dx*eps_0);

% Exy_xp
Cexye_xp=1;
Cexyhy_xp=dt/(dx*eps_0);

% Exy_yn
Cexye_yn=(2*eps_0-dt*sigmay_yn)./(2*eps_0+dt*sigmay_yn);
Cexyhy_yn=(2*dt/dy)./(2*eps_0+dt*sigmay_yn);

% Exy_yp
Cexye_yp=(2*eps_0-dt*sigmay_yp)./(2*eps_0+dt*sigmay_yp);
Cexyhy_yp=(2*dt/dy)./(2*eps_0+dt*sigmay_yp);

% Exy_zn
Cexye_zn=1;
Cexyhy_zn=dt/(dy*eps_0);

% Exy_zp
Cexye_zp=1;
Cexyhy_zp=dt/(dy*eps_0);

% Exz_xn
Cexze_xn=1;
Cexzhy_xn=-dt/(dz*eps_0);

% Exz_xp
Cexze_xp=1;
Cexzhy_xp=-dt/(dz*eps_0);

% Exz_yn
Cexze_yn=1;
Cexzhy_yn=-dt/(dz*eps_0);

% Exz_yp
Cexze_yp=1;
Cexzhy_yp=-dt/(dz*eps_0);

% Exz_zn
Cexze_zn=(2*eps_0-dt*sigmaz_zn)./(2*eps_0+dt*sigmaz_zn);
Cexzhy_zn=-(2*dt/dz)./(2*eps_0+dt*sigmaz_zn);

% Exz_zp
Cexze_zp=(2*eps_0-dt*sigmaz_zp)./(2*eps_0+dt*sigmaz_zp);
Cexzhy_zp=-(2*dt/dz)./(2*eps_0+dt*sigmaz_zp);

% Eyx、Eyz
% Eyx_xn
Ceyxe_xn=(2*eps_0-dt*sigmax_xn)./(2*eps_0+dt*sigmax_xn);
Ceyxhy_xn=-(2*dt/dx)./(2*eps_0+dt*sigmax_xn);

% Eyx_xp
Ceyxe_xp=(2*eps_0-dt*sigmax_xp)./(2*eps_0+dt*sigmax_xp);
Ceyxhy_xp=-(2*dt/dx)./(2*eps_0+dt*sigmax_xp);

% Eyx_yn
Ceyxe_yn=1;
Ceyxhy_yn=-dt/(dx*eps_0);

% Eyx_yp
Ceyxe_yp=1;
Ceyxhy_yp=-dt/(dx*eps_0);

% Eyx_zn
Ceyxe_zn=1;
Ceyxhy_zn=-dt/(dx*eps_0);

% Eyx_zp
Ceyxe_zp=1;
Ceyxhy_zp=-dt/(dx*eps_0);

% Eyz_xn
Ceyze_xn=1;
Ceyzhy_xn=dt/(dx*eps_0);

% Eyz_xp
Ceyze_xp=1;
Ceyzhy_xp=dt/(dz*eps_0);

% Eyz_yn
Ceyze_yn=1;
Ceyzhy_yn=dt/(dz*eps_0);

% Eyz_yp
Ceyze_yp=1;
Ceyzhy_yp=dt/(dz*eps_0);

% Eyz_zn
Ceyze_zn=(2*eps_0-dt*sigmaz_zn)./(2*eps_0+dt*sigmaz_zn);
Ceyzhy_zn=(2*dt/dz)./(2*eps_0+dt*sigmaz_zn);

% Eyz_zp
Ceyze_zp=(2*eps_0-dt*sigmaz_zp)./(2*eps_0+dt*sigmaz_zp);
Ceyzhy_zp=(2*dt/dz)./(2*eps_0+dt*sigmaz_zp);

sigma_star_xn=sigma_starx(1:Npml,1:ny,1:nz);
sigma_star_xp=sigma_starx(nxp1-Npml:nx,1:ny,1:nz);
sigma_star_yn=sigma_stary(1:nx,1:Npml,1:nz);
sigma_star_yp=sigma_stary(1:nx,nyp1-Npml:ny,1:nz);
sigma_star_zn=sigma_starz(1:nx,1:ny,1:Npml);
sigma_star_zp=sigma_starz(1:nx,1:ny,nzp1-Npml:nz);

% Hzx、Hzy
% Hzx_xn
Chzxh_xn=(2*mu0-dt*sigma_star_xn)./(2*mu0+dt*sigma_star_xn);
Chzxey_xn=-(2*dt/dx)./(2*mu0+dt*sigma_star_xn);

% Hzx_xp
Chzxh_xp=(2*mu0-dt*sigma_star_xp)./(2*mu0+dt*sigma_star_xp);
Chzxey_xp=-(2*dt/dx)./(2*mu0+dt*sigma_star_xp);

% Hzx_yn
Chzxh_yn=1;
Chzxey_yn=-dt/(dx*mu0);

% Hzx_yp
Chzxh_yp=1;
Chzxey_yp=-dt/(dx*mu0);

% Hzx_zn
Chzxh_zn=1;
Chzxey_zn=-dt/(dx*mu0);

% Hzx_zp
Chzxh_zp=1;
Chzxey_zp=-dt/(dx*mu0);

% Hzy_xn
Chzyh_xn=1;
Chzyey_xn=dt/(dy*mu0);

% Hzy_xp
Chzyh_xp=1;
Chzyey_xp=dt/(dy*mu0);

% Hzy_yn
Chzyh_yn=(2*mu0-dt*sigma_star_yn)./(2*mu0+dt*sigma_star_yn);
Chzyey_yn=(2*dt/dy)./(2*mu0+dt*sigma_star_yn);

% Hzy_yp
Chzyh_yp=(2*mu0-dt*sigma_star_yp)./(2*mu0+dt*sigma_star_yp);
Chzyey_yp=(2*dt/dy)./(2*mu0+dt*sigma_star_yp);

% Hzy_zn
Chzyh_zn=1;
Chzyey_zn=dt/(dy*mu0);

% Hzy_zp
Chzyh_zp=1;
Chzyey_zp=dt/(dy*mu0);

% Hxy、Hxz
% Hxy_xn
Chxyh_xn=1;
Chxyey_xn=-dt/(dx*mu0);

% Hxy_xp
Chxyh_xp=1;
Chxyey_xp=-dt/(dx*mu0);

% Hxy_yn
Chxyh_yn=(2*mu0-dt*sigma_star_yn)./(2*mu0+dt*sigma_star_yn);
Chxyey_yn=-(2*dt/dy)./(2*mu0+dt*sigma_star_yn);

% Hxy_yp
Chxyh_yp=(2*mu0-dt*sigma_star_yp)./(2*mu0+dt*sigma_star_yp);
Chxyey_yp=-(2*dt/dy)./(2*mu0+dt*sigma_star_yp);

% Hxy_zn
Chxyh_zn=1;
Chxyey_zn=-dt/(dy*mu0);

% Hxy_zp
Chxyh_zp=1;
Chxyey_zp=-dt/(dy*mu0);

% Hxz_xn
Chxzh_xn=1;
Chxzey_xn=dt/(dz*mu0);

% Hxz_xp
Chxzh_xp=1;
Chxzey_xp=dt/(dz*mu0);

% Hxz_yn
Chxzh_yn=1;
Chxzey_yn=dt/(dz*mu0);

% Hxz_yp
Chxzh_yp=1;
Chxzey_yp=dt/(dz*mu0);

% Hxz_zn
Chxzh_zn=(2*mu0-dt*sigma_star_zn)./(2*mu0+dt*sigma_star_zn);
Chxzey_zn=(2*dt/dz)./(2*mu0+dt*sigma_star_zn);

% Hxz_zp
Chxzh_zp=(2*mu0-dt*sigma_star_zp)./(2*mu0+dt*sigma_star_zp);
Chxzey_zp=(2*dt/dz)./(2*mu0+dt*sigma_star_zp);

% Hyx、Hyz
% Hyx_xn
Chyxh_xn=(2*mu0-dt*sigma_star_xn)./(2*mu0+dt*sigma_star_xn);
Chyxey_xn=(2*dt/dx)./(2*mu0+dt*sigma_star_xn);

% Hyx_xp
Chyxh_xp=(2*mu0-dt*sigma_star_xp)./(2*mu0+dt*sigma_star_xp);
Chyxey_xp=(2*dt/dx)./(2*mu0+dt*sigma_star_xp);

% Hyx_yn
Chyxh_yn=1;
Chyxey_yn=dt/(dx*mu0);

% Hyx_yp
Chyxh_yp=1;
Chyxey_yp=dt/(dx*mu0);

% Hyx_zn
Chyxh_zn=1;
Chyxey_zn=dt/(dx*mu0);

% Hyx_zp
Chyxh_zp=1;
Chyxey_zp=dt/(dx*mu0);

% Hyz_xn
Chyzh_xn=1;
Chyzey_xn=-dt/(dx*mu0);

% Hyz_xp
Chyzh_xp=1;
Chyzey_xp=-dt/(dz*mu0);

% Hyz_yn
Chyzh_yn=1;
Chyzey_yn=-dt/(dz*mu0);

% Hyz_yp
Chyzh_yp=1;
Chyzey_yp=-dt/(dz*mu0);

% Hyz_zn
Chyzh_zn=(2*mu0-dt*sigma_star_zn)./(2*mu0+dt*sigma_star_zn);
Chyzey_zn=-(2*dt/dz)./(2*mu0+dt*sigma_star_zn);

% Hyz_zp
Chyzh_zp=(2*mu0-dt*sigma_star_zp)./(2*mu0+dt*sigma_star_zp);
Chyzey_zp=-(2*dt/dz)./(2*mu0+dt*sigma_star_zp);

