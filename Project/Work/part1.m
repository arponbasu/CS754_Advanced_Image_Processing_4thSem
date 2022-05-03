n = 256;
m = 100;

k_arr = 10:2:60;
l = size(k_arr,2);
prob_arr = zeros(l,1,'double');
for i = 1:l
    disp(i)
    prob_arr(i) = rp(k_arr(i),0.01,m,n,2,1000/l);
end

figure();
xi = [k_arr(1):2:k_arr(end)];
vid=interp1(k_arr,prob_arr,xi,'spline');
plot(k_arr,prob_arr,'b*')
hold on
plot(xi,vid,'r')




function recv_prob = rp(k,eps,m,n,iter,trials)
    succ = 0;
    for i = 1:trials
        x = init_vector(zeros(n,1,'double'),k);
        x_r = expt_run(x, eps, iter, m, n);
        d = max(abs(x-x_r));
        if d <= 0.001
            succ = succ + 1;
        end 
    end
    recv_prob = succ/trials;
end




function x_recons = expt_run (x, eps, iter, m, n)
    Phi = randn(m, n);
    x_recons = zeros(n,1,'double');
    y = Phi*x;
    W = eye(n);
    for i = 1:iter
        l = min_x(y,Phi,W);
        %norm(l - x_recons)
        x_recons = l;
        W = diag(invert(x_recons,eps,n));
    end
end

function wt_arr = invert(x, eps, n)
    wt_arr = zeros(n,1,'double');
    for i = 1:n
        wt_arr(i) = 1/(abs(x(i)) + eps);
    end
end

function x = min_x(y, Phi, W)
    [~,n] = size(Phi);
    cvx_begin quiet  
        variables x(n) 
        minimize( norm (W * x, 1) );
        subject to
            y == Phi*x;
    cvx_end
end


function x_0 = init_vector(x_0, k)
    y = randsample(size(x_0, 1), k);
    for i = 1:size(y,1)
        x_0(y(i)) = randn;
    end
end

