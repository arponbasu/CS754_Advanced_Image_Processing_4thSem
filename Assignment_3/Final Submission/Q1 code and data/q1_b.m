rng(4);

orig_image = double(imread("barbara256.png"));
%Padded so that averaging overlapping patches becomes easier later on
padded_orig_image = zeros(size(orig_image,1)+14,size(orig_image,2)+14);
padded_orig_image(8:7+size(orig_image,1),8:7+size(orig_image,2)) = orig_image; 
U = kron(dctmtx(8)',dctmtx(8)');

[rows, cols] = size(padded_orig_image);

phi = randn(32, 64);
reconstructed_image = zeros(size(padded_orig_image));
A = phi*U; %Our LASSO matrix has a measurement component too, this time
max_eigen = max(eig(A'*A));
alpha = max_eigen + 5; %5 is added to maximum eigenvalue to avoid numerical stability issues 
lambda = 1000;

for i=1:rows-7
    for j=1:cols-7
        patch = padded_orig_image(i:i+7, j:j+7);
        compressed_image = phi*patch(:);
        reconstructed_image(i:i+7, j:j+7) = reconstructed_image(i:i+7, j:j+7) + reshape(U*ISTA(compressed_image,A,lambda,alpha,200), 8, 8);
    end
end   

reconstructed_image = reconstructed_image(8:263, 8:263)/64; %because of padding averaging is same for all
rmse = err(orig_image(:),reconstructed_image(:)) %Reconstruction error


imshow(uint8(orig_image))
	title('Original Image', 'FontSize', 14);
	drawnow;
	pause(2)
imshow(uint8(reconstructed_image))
	title('Reconstructed Image', 'FontSize', 14);
	drawnow;
    pause(2)
    imwrite(uint8(reconstructed_image), "Reconstructed_Image_1b.jpg");




clear;
%utility functions below
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
