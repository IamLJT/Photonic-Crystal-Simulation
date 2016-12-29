for i=1:nt-1000
    h=imagesc(1e+6*dy*(1:1:ny),(dx*1e+6*(1:1:nx))',Sumhz_t(:,:,i)');%,[-0.05,0.05]);
        colorbar;
        set(h,'AlphaData',epsilon(:,:,dst_layer)'/epsilon0);
%         Sum1(:,:,n-3000)=Sumx(:,:);
%         Sum2(:,:,n-3000)=Sumy(:,:);
        shading interp;
        getframe;
end