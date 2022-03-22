% clear all;
addpath('./l1_ls_matlab');
% Adding Paths
file1='./Images/slice_50.png';
file2='./Images/slice_51.png';

% Reading Images
img1_temp=imread(file1);
img2_temp=imread(file2); 

% Padding Images to get overall square images
img1 = zeros(size(img1_temp,2),size(img1_temp,2));
img1((size(img1_temp,2) -size(img1_temp,1))/2:(size(img1_temp,2) +size(img1_temp,1))/2-1,:) = img1_temp;

img2 = zeros(size(img2_temp,2),size(img2_temp,2));
img2((size(img2_temp,2) -size(img2_temp,1))/2:(size(img2_temp,2) +size(img2_temp,1))/2-1,:) = img2_temp;

% Displaying the original Images
figure('name','Original Image 1');
imshow(uint8(img1));
figure('name','Original Image 2');
imshow(uint8(img2));

% Initialising the theta vector, using angles as 0,10,20...
totalAngles=18;
theta_1= zeros(1,18);
for i=0:totalAngles-1
    theta_1(i+1) = 10*i;
end

% Part (a)
[R,~] = radon(img1,theta_1);
% Calculating the image back from the radon transformation.
I_ramlak = iradon(R,theta_1,'nearest','Ram-Lak');

Mean_Squared_Error_1 = mean((I_ramlak(1:217,1:217) - img1).^2,'all') / mean(img1.^2, 'all');

fig = figure('name', 'RamLak');
imshow(uint8(I_ramlak));
saveas(fig,'RamLak_tomography.png');


% Part b

% We need to construct the A matrix now. This matrix will have the number
% of rows as number of rows in R and the number of columns as the number of
% pixels in the original image

m_A = size(R, 1);
n_A = numel(img1);

% This matrix (and its transpose) have been represented using classes and
% using operator overloading of matrix multiplication.

A=A_class(m_A,n_A,theta_1);
At=At_class(n_A,m_A,theta_1);

% Finally we set value of lambda and vectorise the Radon transform of the
% image to obtain a compressive sensing suitable format.

y=reshape(R,[],1);
lambda = 10;
[result,~]=l1_ls(A,At,m_A,n_A,y,lambda);

% The final result is again a vector and is resized to the original shape
% of the image. idct2 helps convert back to original basis.
I_CS=idct2(reshape(result,size(img1,1),[]));

Mean_Squared_Error_2 = mean((I_CS - img1).^2,'all') / mean(img1.^2, 'all');

fig = figure('name', 'Compressive Sensing');
imshow(uint8(I_CS));
saveas(fig,'CS_tomography.png');



% Part c
% Here we will be using the similarities of two consecutive frames to
% obtain the tomographic reconstructions of the images.

% We obtain a new set of angles for the measurement of second slice
theta_2= zeros(1,18);
for i=1:totalAngles
    theta_2(i) = 10*i-5;
end

% The same procedure as part (a) is repeated for slice 2.
[R_2,~] = radon(img2,theta_2);

num_slices = 2;

% We initialise our new Coupled matrices, for coupled-CS reconstruction.
A_coupled = A_coupled_class(m_A, n_A,{theta_1,theta_2});
At_coupled = At_coupled_class(n_A, m_A,{theta_1,theta_2});

% The radon transform is also combined and vectorised as 1, just like in
% the slides.
R = [R R_2];
y=reshape(R,[],1);

% We define the new sizes, row and column sizes. Initialise lambda.
m_A_couple=size(y,1);
n_A_couple=2*n_A;
lambda = 10;

% Similar to part b, we obtain the result from the l1 lasso minimization
[result,~]=l1_ls(A_coupled,At_coupled,m_A_couple,n_A_couple,y,lambda);

% This is critical. Here we get that both the coefficient for the original
% image as well as the coefficients telling about the difference between
% the two frames.
Coeff_result=reshape(result(1:n_A,1),sqrt(n_A),[]);
DeltaCoeff_result=reshape(result(n_A+1:2*n_A,1),sqrt(n_A),[]);

% What follows is simple basis conversion and displaying of results.
I_coupled_CS_1=idct2(Coeff_result);

Mean_Squared_Error_3_1 = mean((I_coupled_CS_1 - img1).^2,'all') / mean(img1.^2, 'all');
fig = figure('name', 'Coupled Compressive Sensing 1');
imshow(uint8(I_coupled_CS_1));
saveas(fig,'Coupled_CS_tomography1.png');



I_coupled_CS_2=idct2(Coeff_result+DeltaCoeff_result);

Mean_Squared_Error_3_2 = mean((I_coupled_CS_2 - img2).^2,'all') / mean(img2.^2, 'all');
fig = figure('name', 'Coupled Compressive Sensing 2');
imshow(uint8(I_coupled_CS_2));
saveas(fig,'Coupled_CS_tomography2.png');
