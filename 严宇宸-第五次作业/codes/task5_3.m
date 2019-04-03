clear;clc;
%����ͼ�񲢽���Ԥ����
test1 = im2double(imread('test3_corrupt.pgm'));
test2 = im2double(imread('test4 copy.bmp'));
Laplacian(test1)
Unmask(test1,20,2)
Laplacian(test2)
Unmask(test2,20,2)
%% ������˹��ͨ�˲���
function Laplacian(Img_in)
[M,N]=size(Img_in);    M=2*M;N=2*N;
u = -M/2:(M/2-1);   v = -N/2:(N/2-1);
[u,v] = meshgrid(u,v);
D = sqrt(u.^2+v.^2);
%����˲�����
H = -4.*pi.^2*(u.^2+v.^2); 
Img_fft= fft2(Img_in,size(H,1),size(H,2));
Img_fft_shift = fftshift(Img_fft);
Img_Laplacian = Img_fft_shift.*H;
Img_out = ifft2(ifftshift(Img_Laplacian));
Img_out = Img_out(1:size(Img_in,1),1:size(Img_in,2));
%���ͼ��
figure;
subplot(2,2,1);imshow(Img_in);title("ԭʼͼ��");
subplot(2,2,2);imshow(log(1+abs(Img_fft_shift)),[]);title('Ƶ��ͼ��');
subplot(2,2,3);plot3(u,v,H);title('������˹�˲���');
subplot(2,2,4);imshow(Img_out);title('�˲����ͼ��');
end

%% Unmask��ͨ�˲���
function Unmask(Img_in,D0,n)
Img_fft_shift=fftshift(fft2(Img_in));
[M,N]=size(Img_fft_shift);
for i=1:M
    for j=1:N
        d=sqrt((i-fix(M/2))^2+(j-fix(N/2))^2);
        if d==0
            H(i,j)=0;
        else
            H(i,j)=1/(1+0.414*(D0/d)^(2*n));
        end
        Img_out(i,j)=(1+H(i,j))*Img_fft_shift(i,j);
    end
end
Img_out=real(ifft2(ifftshift(Img_out)));
%���ͼ��
figure;
subplot(2,2,1);imshow(Img_in);title("ԭʼͼ��");
subplot(2,2,2);imshow(log(1+abs(Img_fft_shift)),[]);title('Ƶ��ͼ��');
subplot(2,2,3);
u=-M/2:(M/2-1);v=-N/2:(N/2-1);
[u,v]=meshgrid(u,v);plot3(u,v,H);title('Unmask�˲���');
subplot(2,2,4);imshow(Img_out);title(['�˲����ͼ��,D0=',num2str(D0),',n=',num2str(n)]);
end