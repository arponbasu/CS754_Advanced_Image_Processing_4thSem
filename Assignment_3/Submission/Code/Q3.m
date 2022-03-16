% clear all;
addpath('./l1_ls_matlab');

file1='./Images/slice_50.png';
file2='./Images/slice_51.png';
img1=imread(file1);
% img1=double(img1);
img2=imread(file2); 
% img2=double(img2);

% figure();
% imshow(img1);
% 
% figure();
% imshow(img2);

totalAngles=18;
theta= zeros(1,18);
for i=0:totalAngles-1
    theta(i+1) = 10*i;
end

[R,xp] = radon(img1,theta);

I_ramlak = iradon(R,theta,'nearest','Ram-Lak');

figure();
imshow(uint8(I_ramlak));



   






















