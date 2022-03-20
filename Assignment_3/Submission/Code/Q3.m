% clear all;
addpath('./l1_ls_matlab');

file1='./Images/slice_50.png';
file2='./Images/slice_51.png';
img1_temp=imread(file1);
% img1=double(img1);
img2_temp=imread(file2); 
% img2=double(img2);
img1 = zeros(size(img1_temp,2),size(img1_temp,2));
img1((size(img1_temp,2) -size(img1_temp,1))/2:(size(img1_temp,2) +size(img1_temp,1))/2-1,:) = img1_temp;

img2 = zeros(size(img2_temp,2),size(img2_temp,2));
img2((size(img2_temp,2) -size(img2_temp,1))/2:(size(img2_temp,2) +size(img2_temp,1))/2-1,:) = img2_temp;

% figure();
% imshow(uint8(img1));
% 
% figure();
% imshow(uint8(img2));

totalAngles=18;
theta_1= zeros(1,18);
for i=0:totalAngles-1
    theta_1(i+1) = 10*i;
end

% Part a
[R,~] = radon(img1,theta_1);

I_ramlak = iradon(R,theta_1,'nearest','Ram-Lak');

% figure();
% imshow(uint8(I_ramlak));


% Part b
% We need to construct the A matrix now. This matrix will have the number
% of rows as number of rows in R and the number of columns as the number of
% pixels in the original image

m_A = size(R, 1);
n_A = numel(img1);

% A=A_class(m_A,n_A,theta_1);
% At=At_class(n_A,m_A,theta_1);
% 
% y=reshape(R,m_A*numel(theta_1),1);
% lambda = 10;
% [result,~]=l1_ls(A,At,m_A,n_A,y,lambda);
% 
% f=idct2(reshape(result,size(img1,1),size(img1,2)));
% save('result_b.mat','f');
% figure();
% imshow(uint8(f));
% 
% Mean_Squared_Error_2 = mean((f - img1).^2,'all') / mean(img1.^2, 'all')



%Part c
% Here we will be using the similarities of two consecutive frames to
% obtain the tomographic reconstructions of the images.

theta_2= zeros(1,18);
for i=1:totalAngles
    theta_2(i) = 10*i-5;
end

[R_2,~] = radon(img2,theta_2);


num_slices = 2;

A_coupled = A_coupled_class(m_A, n_A,{theta_1,theta_2});
At_coupled = At_coupled_class(n_A, m_A,{theta_1,theta_2});

R = [R R_2];

y=reshape(R,[],1);
m_A_couple=size(y,1);
n_A_couple=2*n_A;
lambda = 4;
[beta_1and2,~]=l1_ls(A_coupled,At_coupled,m_A_couple,n_A_couple,y,lambda);

dctCoeff_1=reshape(beta_1and2(1:n_A,1),sqrt(n_A),sqrt(n_A));
dctCoeff_2=reshape(beta_1and2(n_A+1:2*n_A,1),sqrt(n_A),sqrt(n_A));


f=idct2(reshape(dctCoeff_1,sqrt(n_A),sqrt(n_A)));
figure();
imshow(uint8(f));
save('result_c1.mat','f');
Mean_Squared_Error_3_1 = mean((f - img1).^2,'all') / mean(img1.^2, 'all')


f=idct2(reshape(dctCoeff_1+dctCoeff_2,sqrt(n_A),sqrt(n_A)));
figure();
imshow(uint8(f));
save('result_c2.mat','f');
Mean_Squared_Error_3_2 = mean((f - img2).^2,'all') / mean(img2.^2, 'all')
































   











%% Helper Functions and classes








