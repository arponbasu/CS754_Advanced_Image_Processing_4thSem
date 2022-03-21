rng(1000);

im = imread('barbara256.png');

Phi = randn(32,64);
alpha = (svds(Phi,1))^2 + 1;
U = kron(dctmtx(8)',dctmtx(8)');

final = zeros(256,256,"double");
cnt = zeros(256,256,"double");



for i = 1:249
    for j = 1:249
        
        r = reconstruct_patch(imcrop(im,[i,j,7,7]),Phi,U,alpha,1,100);
        for k1 = i:i+7
            for k2 = j:j+7
                cnt(k1,k2) = cnt(k1,k2) + 1;
                final(k1,k2) = final(k1,k2) + r(k1-i+1,k2-j+1); 
            end
        end
    
    
    end
end

final = final./cnt;
final = transpose(final);
imshow(uint8(final))
err(double(im(:)),double(final(:)))
imwrite(uint8(final),'1b.jpg')
clear;

function r = reconstruct_patch (patch, Phi, U, alpha, lambda, iter)
    p = patch(:); 
    p = double(p);
    
    y = Phi*p;
    theta = ISTA(y,Phi*U,lambda,alpha,iter);
    Ut = U*theta;
    
    r = reshape(Ut,size(patch,1),size(patch,2));
end


function theta = ISTA(y,A,lambda,alpha,iter)
    theta = 0*A'*y; 
    T = lambda/(2*alpha);
    for k = 1:iter
        theta = soft(theta + (A'*(y - A*theta))/alpha, T);
    end
end

function s = soft(x, T)
    s = sign(x).*(max(0, abs(x)-T));
end

function e = err(x, y)
    e = norm(x-y)/norm(x);
end