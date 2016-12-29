%   默认正弦波,且只有Ex方向，为线偏光
a=source_mid_x;
b=source_mid_y;
c0=floor(source_size/2);
% if w*n*dt>=2*pi
%     return;
% end
if source_type==1
%     Ex(a-c:a+c,b-c:b+c,source_layer)=Ex(a-c:a+c,b-c:b+c,source_layer)+...
%         exp(1i*((w*n*dt-pi/2)));
    if source_shape==2
        Ex(a-c0:a+c0,b-c0:b+c0,source_layer)=0.5*exp(1i*((w*n*dt-pi/2)));
        Ey(a-c0:a+c0,b-c0:b+c0,source_layer)=0;
    elseif source_shape==1
        for ii=a-c0:a+c0
            for jj=b-c0:b+c0
                if (ii-a)^2+(jj-b)^2<c0^2
                    Ex(ii,jj,source_layer)=0.5*exp(1i*((w*n*dt-pi/2)));
                    Ey(ii,jj,source_layer)=0;
                end
            end
        end
    end
%     Ez(a-c:a+c,b-c:b+c,source_layer)=0;
elseif source_type==2
    if source_shape==2
        [xx,yy]=meshgrid((-c0:c0),(-c0:c0));
        Ex(a-c0:a+c0,b-c0:b+c0,source_layer)=0.5*exp(1i*(w*n*dt-pi/2))*exp(-(xx.^2+yy.^2)/30);
        Ey(a-c0:a+c0,b-c0:b+c0,source_layer)=0;
    elseif source_shape==1
        for ii=a-c0:a+c0
            for jj=b-c0:b+c0
                if (ii-a)^2+(jj-b)^2<c0^2
                    Ex(ii,jj,source_layer)=0.5*exp(1i*((w*n*dt-pi/2)))*exp(-((ii-a)^2+(jj-b)^2)/30);
                    Ey(ii,jj,source_layer)=0;
                end
            end
        end
    end
    
    
%     Ex(a-c:a+c,b-c:b+c,source_layer)=Ex(a-c:a+c,b-c:b+c,source_layer)+...
%         0.5*exp(1i*((w*n*dt-pi/2)))*exp(-(xx.^2+yy.^2)/10);
% if(dt*w*n<=pi)
%     Ex(a-c:a+c,b-c:b+c,source_layer)=0.5*exp(1i*((w*n*dt)))*exp(-(xx.^2+yy.^2)/10);
% end
%     Ex(a-c:a+c,b-c:b+c,source_layer)=0.5*exp(1i*(w*n*N_w/m_T*dt-pi/2))*exp(-(xx.^2+yy.^2)/10);
    
    
    
%     Ex(a-c:a+c,b-c:b+c,source_layer)=0.5*sin(w*n*dt)*exp(-(xx.^2+yy.^2)/10);
%     Ey(a-c:a+c,b-c:b+c,source_layer)=0;
%     Ez(a-c:a+c,b-c:b+c,source_layer)=0;
end
