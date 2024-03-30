f=double(imread('01.png'));%double将读的图转化成双精度,灰度图单通道，彩色图三通道。
f=f(1:256,1:256);
f=f/max(f(:));%归一化处理
stream = RandStream('mcg16807', 'Seed',0);%创建单个随机数流
RandStream.setGlobalStream(stream);
stdn=.05;
noise=stdn*randn(size(f));
%Add noise
y=f+noise;%y是256*256*3的
%figure(1), imshow(f,[]); figure(2),imshow(y,[]);
%Denoise color image
lambda=0.039; % regularization parameter
[xST,P,fun_val,ISNR]=proxSTV(y,lambda,'verbose',true,'img',f,'maxiter',100,'kernel',fspecial('gaussian',[3 3],0.5),'L',8/1.25,'snorm','nuclear','project',@(x)BoxProjection(x,[0 1]),'showfig',1);
PSNR=psnr(xST,f);
SSIM=ssim(xST,f);
imshow(xST,[],'border','tight','initialmagnification','fit');
saveas(gcf,'WSTVairplane','epsc');