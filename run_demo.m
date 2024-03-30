f=double(imread('01.png'));%double������ͼת����˫����,�Ҷ�ͼ��ͨ������ɫͼ��ͨ����
f=f(1:256,1:256);
f=f/max(f(:));%��һ������
stream = RandStream('mcg16807', 'Seed',0);%���������������
RandStream.setGlobalStream(stream);
stdn=.05;
noise=stdn*randn(size(f));
%Add noise
y=f+noise;%y��256*256*3��
%figure(1), imshow(f,[]); figure(2),imshow(y,[]);
%Denoise color image
lambda=0.039; % regularization parameter
[xST,P,fun_val,ISNR]=proxSTV(y,lambda,'verbose',true,'img',f,'maxiter',100,'kernel',fspecial('gaussian',[3 3],0.5),'L',8/1.25,'snorm','nuclear','project',@(x)BoxProjection(x,[0 1]),'showfig',1);
PSNR=psnr(xST,f);
SSIM=ssim(xST,f);
imshow(xST,[],'border','tight','initialmagnification','fit');
saveas(gcf,'WSTVairplane','epsc');