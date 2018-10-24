%GaussianMove2D.m
%Moves Gaussian curve across time...hopefully
clear;clc;close all;
i=sqrt(-1);
w=rand()*1000;
k=rand()*1000;
theta=rand()*2*pi;
sigma=rand*1000;
absV=w/k;
vx=absV*cos(theta);
vy=absV*sin(theta);
for n=1:201
    x(n)=n-101;
    y(n)=-n+101;
end
%{
Mx=x
My=y'
for m=1:200
    Mx=[Mx;x];
    My=[My,y'];
end
%}
ender=round(100/absV);
if (ender>50)
    ender=50;
end
for t=1:ender
    figure(t)
    for xIndex=1:201
        X=x(xIndex);
        for yIndex=1:201
            Y=y(yIndex);
        Packet(xIndex,yIndex)=exp(-((X-vx.*t).^2)./sigma).*exp(-((Y-vy.*t).^2)./sigma);
        Sinusoid(xIndex,yIndex)=cos((k*(cos(theta)*X+sin(theta)*Y)-w*t));
        end
    end
WavePacket=Packet.*Sinusoid;
contour(x,y,WavePacket)
end