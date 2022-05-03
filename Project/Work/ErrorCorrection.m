n = 128;
m = 512;

k_arr = 128:4:228;

t = cputime;
gen_array(k_arr,0.001,m,n,4,105,'log');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.01,m,n,4,105,'log');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,4,105,'log');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,1,m,n,4,105,'log');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,10,m,n,4,105,'log');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,1,105,'log'); %unweighted l1, epsilon value irrelevant here
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.001,m,n,4,105,'atan');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.01,m,n,4,105,'atan');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,4,105,'atan');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,1,m,n,4,105,'atan');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,10,m,n,4,105,'atan');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,1,105,'atan'); %unweighted l1, epsilon value irrelevant here
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.001,m,n,4,105,'tanh');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.01,m,n,4,105,'tanh');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,4,105,'tanh');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,1,m,n,4,105,'tanh');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,10,m,n,4,105,'tanh');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,1,105,'tanh'); %unweighted l1, epsilon value irrelevant here
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.001,m,n,4,105,'sigmoid');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.01,m,n,4,105,'sigmoid');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,4,105,'sigmoid');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,1,m,n,4,105,'sigmoid');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,10,m,n,4,105,'sigmoid');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,1,105,'sigmoid'); %unweighted l1, epsilon value irrelevant here
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.1,m,n,2,105,'log'); 
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,8,105,'log');
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.1,m,n,2,105,'atan'); 
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,8,105,'atan');
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.1,m,n,2,105,'tanh'); 
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,8,105,'tanh');
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.1,m,n,2,105,'sigmoid'); 
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,8,105,'sigmoid');
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function prob_arr = gen_array(k_arr,beta,m,n,iter,trials,mode)
    l = size(k_arr,2);
    prob_arr = zeros(l,1,'double');
    for i = 1:l
        disp(i)
        prob_arr(i) = rp(k_arr(i),beta,m,n,iter,trials/l,mode);
    end
    a1 = num2str(beta);
    a2 = num2str(iter);
    a3 = num2str(trials);
    a4 = num2str(mode);
    a = append(a1, '-', a2, '-', a3, '-', a4, '-part5', '.mat');
    save(a,'prob_arr','-mat');
end


function recv_prob = rp(k,beta,m,n,iter,trials,mode)
    succ = 0;
    for i = 1:trials
        x = randn(n,1);
        x_r = expt_run(x, beta, iter, m, n, mode, k);
        d = max(abs(x-x_r));
        if d <= 0.001
            succ = succ + 1;
        end 
    end
    recv_prob = succ/trials;
end




function x_recons = expt_run (x, beta, iter, m, n, mode, k)
    Phi = randn(m, n);
    x_recons = zeros(n,1,'double');
    y = toggle(Phi*x, k);
    W = eye(m);
    eps = beta*std(y);
    for i = 1:iter
        l = min_x(y,Phi,W);
        x_recons = l;
        W = diag(invert(y, Phi, x_recons,eps,n,mode));
    end
end

function wt_arr = invert(y, Phi, x, eps, n, mode)
    wt_arr = zeros(n,1,'double');
    r = y - Phi*x;
    if mode == "log"
        wt_arr = 1./(abs(r) + eps);
    elseif mode == "atan"
        wt_arr = 1./(r.^2 + eps^2);    
    elseif mode == "tanh"
        wt_arr = (sech(r/eps)).^2;
    elseif mode == "sigmoid"
        z = sgm(r/eps);
        wt_arr = z.*(1-z);        
    end
end

function y = sgm(x)
    y = 1./(1 + exp(-x));
end

function x = min_x(y, Phi, W)
    [~,n] = size(Phi);
    cvx_begin quiet  
        variables x(n) 
        minimize( norm (W * (y - Phi*x), 1) );
    cvx_end
end

function x_0 = toggle(x_0, k)
    y = randsample(size(x_0, 1), k);
    for i = 1:size(y,1)
        x_0(y(i)) = -x_0(y(i));
    end
end

