# Matlab1
(1)粒度测定
I = imread('1.png');
figure;imshow(I)
I=rgb2gray(I);
claheI = adapthisteq(I,'NumTiles',[10,10]);
claheI = imadjust(claheI);
figure;imshow(claheI);
（2）计算粒度大小的总体分布
for counter=0:22
    remain=imopen(claheI,strel('disk',counter));
intensity_area(counter+1)=sum(remain(:));
end
figure
plot(intensity_area,'m - *'),
grid on;
（3）计算不同半径下的粒度分布
intensity_area_prime=diff(intensity_area);
figure;
plot(intensity_area_prime,'m - *'),
grid on;
title('Granulometry(Size Distribution)of Snowflakes');
set(gca,'xtick',[0 2 4 6 8 10 12 14 16 18 20 22]);
xlabel('radius of spore(pixels)');
ylabel('Sum of pixel values in spore as a function of radius');
open5=imopen(claheI,strel('disk',5));
open6=imopen(claheI,strel('disk',6));
rad5=imsubtract(open5,open6);
figure;imshow(rad5,[]);
(4)转换为灰度图像
img= imread('1.png');
figure;imshow(img);
img=rgb2gray(img);     
figure;imshow(img);  
（5）平滑处理
img= imread('1.png');
R=img(:,:,1);
G=img(:,:,2);
B=img(:,:,3);
R1=gaosi(R);
G1=gaosi(G);
B1=gaosi(B);
RGB(:,:,1)=R1(:,:,1);
RGB(:,:,2)=G1(:,:,1);
RGB(:,:,3)=B1(:,:,1);
figure();
subplot(1,2,1);
imshow(img);
title('原始图像');
subplot(1,2,2);
imshow(RGB);
title('低通滤波');
其中gaosi()是一个调用函数，是自定义函数，内容如下：
function [ img] = gaosi(image)
d0=50;%阈值
[M,N]=size(image);
img_f=fft2(double(image));%傅里叶变换得到频谱
img_f=fftshift(img_f);%移到中间
m_mid=floor(M/2);%中心点坐标
    n_mid=floor(N/2);  
    h = zeros(M,N);%高斯低通滤波器构造
    for i = 1:M
        for j = 1:N
            d = ((i-m_mid)^2+(j-n_mid)^2);
            h(i,j) = exp(-d/(2*(d0^2)));      
        end
    end
    img_lpf = h.*img_f;
    img_lpf=ifftshift(img_lpf);    %中心平移回原来状态
    img_lpf=uint8(real(ifft2(img_lpf)));  %反傅里叶变换,
    img = img_lpf;
end
（6）canny算子算法
clc
clear all
img_in=imread('1.png');
img_in=rgb2gray(img_in);
figure,imshow(img_in);
title('原图');
[rows,cols]=size(img_in);
thresh=graythresh(img_in);
img_bw=im2bw(img_in,thresh);
%%step1:高斯滤波
template=fspecial('gaussian',3,0.8);%生成一个3*3的高斯模板，标准差选0.8
img_filt=imfilter(img_bw,template);
 
%%step2:计算梯度（幅度和方向）
%Prewitt梯度模板
dx = [-1 -1 -1;0 0 0;1 1 1];%x方向的梯度模板
dy = [-1 0 1; -1 0 1;-1 0 1];%y方向的梯度模板
img_filt=double(img_filt);
 
grad_x=conv2(img_filt,dx,'same');
grad_y=conv2(img_filt,dy,'same');
grad=sqrt((grad_x.^2)+(grad_y.^2));%梯度幅值图像
figure,imshow(grad);
title('梯度幅值图');
grad_dir=atan2(grad_y,grad_x);%获取梯度方向弧度
grad_dir=grad_dir*180/pi;
%%step3:对梯度幅值进行非极大值抑制
for i = 1:rows
    for j = 1:cols
        if((grad_dir(i,j)>=-22.5 && grad_dir(i,j)<=22.5) || (grad_dir(i,j)>=157.5 && grad_dir(i,j)<=180)...
                                       ||(grad_dir(i,j)<=-157.5 && grad_dir(i,j)>=-180) )
            grad_dir(i,j) = 0;
        elseif((grad_dir(i,j) >= 22.5) && (grad_dir(i,j) < 67.5) || (grad_dir(i,j) <= -112.5) && (grad_dir(i,j) > -157.5))
            grad_dir(i,j) = -45;
        elseif((grad_dir(i,j) >= 67.5) && (grad_dir(i,j) < 112.5) || (grad_dir(i,j) <= -67.5) && (grad_dir(i,j) >- 112.5))
            grad_dir(i,j) = 90;
        elseif((grad_dir(i,j) >= 112.5) && (grad_dir(i,j) < 157.5) || (grad_dir(i,j) <= -22.5) && (grad_dir(i,j) > -67.5))
            grad_dir(i,j) = 45;  
        end
    end
end
 
Nms = zeros(rows,cols);
for i = 2:rows-1
    for j= 2:cols-1
        if (grad_dir(i,j) == 90 && grad(i,j) == max([grad(i,j), grad(i,j+1), grad(i,j-1)]))
            Nms(i,j) = grad(i,j);
        elseif (grad_dir(i,j) == -45 && grad(i,j) == max([grad(i,j), grad(i+1,j-1), grad(i-1,j+1)]))
            Nms(i,j) = grad(i,j);
        elseif (grad_dir(i,j) == 0 && grad(i,j) == max([grad(i,j), grad(i+1,j), grad(i-1,j)]))
            Nms(i,j) = grad(i,j);
        elseif (grad_dir(i,j) == 45 && grad(i,j) == max([grad(i,j), grad(i+1,j+1), grad(i-1,j-1)]))
            Nms(i,j) = grad(i,j);
        end;end;end;
figure,imshow(Nms);
title('非极大值抑制图');
img_out=zeros(rows,cols);%定义一个双阈值图像
YH_L=0.1*max(max(Nms));
YH_H=0.3*max(max(Nms));
for i = 1:rows
    for j = 1:cols
        if(Nms(i,j)<YH_L)
           img_out(i,j)=0;
        elseif(Nms(i,j)>YH_H)
                img_out(i,j)=1;        
        elseif ( Nms(i+1,j) < YH_H || Nms(i-1,j) < YH_H ||Nms(i,j+1)<YH_H||Nms(i,j-1)< YH_H ||...
                Nms(i-1,j-1) < YH_H || Nms(i-1, j+1) < YH_H ||Nms(i+1,j+1) < YH_H || Nms(i+1, j-1) < YH_H)
                   img_out(i,j) = 1;   
        end;end;end
bw=edge(img_bw,'canny');
figure,imshow(img_out);title('结果图');figure,imshow(bw);title('Canny算子效果图');
（7）腐蚀
img=imread('2.png');
B=[0 1 0;1 1 1;0 1 0];
i1=imerode(img,B);
i2=imerode(i1,B);
i3=imerode(i2,B);
imshow(img),figure,imshow(i1),figure,imshow(i2),figure,imshow(i3);
（8）膨胀
X=imread('3.png');
B=[0 1 0;1 1 1;0 1 0];
Z=imdilate(X,B);
imshow(X),figure,imshow(Z);
总结：首先对图像进行一系列预处理，再进行特征提取。
