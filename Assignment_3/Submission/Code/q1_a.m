rng(1);


orig_image = double(imread("barbara256.png"));
padded_orig_image = padarray(orig_image,[7, 7],0,"both"); %Padded so that averaging overlapping patches becomes easier later on
noise = sqrt(3)*randn(size(padded_orig_image)); %Noise
noisy_image = padded_orig_image + noise;

U = kron(dctmtx(8)',dctmtx(8)'); %Basis Matrix
 
[rows, cols] = size(padded_orig_image);

reconstructed_image = zeros(rows, cols);


alpha = 2; %Since U is orthonormal, the max eigval of U'*U is 1. We take alpha = 1 + 1 to avoid numerical 
%stability issues
lambda = 1;

for i=1:rows-7
    for j=1:cols-7
        noisy_image_patch = noisy_image(i:i+7, j:j+7);
        noisy_image_patch = noisy_image_patch(:);
        theta_est = ISTA(noisy_image_patch,U,lambda,alpha,100);
        reconstructed_image(i:i+7, j:j+7) = reconstructed_image(i:i+7, j:j+7) + reshape(U*theta_est, 8, 8); 
    end
end   

reconstructed_image = reconstructed_image(8:263, 8:263)/64;%Taking averages of overlapping patches


noisy_image = noisy_image(8:263, 8:263); %The padding is removed
orig_error = err(orig_image(:), noisy_image(:)) %Error due to noise
rmse = err(orig_image(:),reconstructed_image(:)) %Error of reconstruction

imshow(uint8(orig_image))
	title('Original Image', 'FontSize', 14);
	drawnow;
	pause(2)
imshow(uint8(noisy_image))
	title('Noisy Image', 'FontSize', 14);
	drawnow;
	pause(2)
    imwrite(uint8(noisy_image), "Noisy_Image_1a.jpg");
imshow(uint8(reconstructed_image))
	title('Reconstructed Image', 'FontSize', 14);
	drawnow;
	pause(2)
    imwrite(uint8(reconstructed_image), "Reconstructed_Image_1a.jpg");

clear;

%Utility functions below

function theta = ISTA(noisy_image,U,lambda,alpha,iter)    
    theta = zeros(size(U,2),1); 
    T = lambda/(2*alpha);
    for i = 1:iter        
        theta = soft(theta + (U'*(noisy_image - U*theta))/alpha, T);
    end
end


function s = soft(x, T)
    s = sign(x).*(max(0, abs(x)-T));
end

function e = err(x, y)
    e = norm(x-y)/norm(x);
end
