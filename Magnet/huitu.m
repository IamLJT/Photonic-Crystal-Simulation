figure(1);
h=imagesc(1e+6*dy*(1:1:ny),(dx*1e+6*(1:1:nx))',Sumx');%,[-0.05,0.05]);
colorbar;
set(h,'AlphaData',epsilon(:,:,dst_layer)'/epsilon0);
shading interp;
title('Ex����ͼ��');
xlabel('x(um)');
ylabel('y(um)');
figure(2);
h=imagesc(1e+6*dy*(1:1:ny),(dx*1e+6*(1:1:nx))',Sumy');%,[-0.05,0.05]);
colorbar;
set(h,'AlphaData',epsilon(:,:,dst_layer)'/epsilon0);
shading interp;
title('Ey����ͼ��');
xlabel('x(um)');
ylabel('y(um)');
figure(3);
h=imagesc(1e+6*dy*(1:1:ny),(dx*1e+6*(1:1:nx))',Sumz');%,[-0.05,0.05]);
colorbar;
set(h,'AlphaData',epsilon(:,:,dst_layer)'/epsilon0);
shading interp;
title('Ez����ͼ��');
xlabel('x(um)');
ylabel('y(um)');
figure(4);
h=imagesc(1e+6*dy*(1:1:ny),(dx*1e+6*(1:1:nx))',Sumhx');%,[-0.05,0.05]);
colorbar;
set(h,'AlphaData',epsilon(:,:,dst_layer)'/epsilon0);
shading interp;
title('Hx����ͼ��');
xlabel('x(um)');
ylabel('y(um)');
figure(5);
h=imagesc(1e+6*dy*(1:1:ny),(dx*1e+6*(1:1:nx))',Sumhy');%,[-0.05,0.05]);
colorbar;
set(h,'AlphaData',epsilon(:,:,dst_layer)'/epsilon0);
shading interp;
title('Hy����ͼ��');
xlabel('x(um)');
ylabel('y(um)');
figure(6);
h=imagesc(1e+6*dy*(1:1:ny),(dx*1e+6*(1:1:nx))',Sumhz');%,[-0.05,0.05]);
colorbar;
set(h,'AlphaData',epsilon(:,:,dst_layer)'/epsilon0);
shading interp;
title('Hz����ͼ��');
xlabel('x(um)');
ylabel('y(um)');