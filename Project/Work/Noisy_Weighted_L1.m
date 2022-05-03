n = 256;
m = 100;

k_arr = 10:2:60;
t = cputime;
gen_array(k_arr,0.001,m,n,4,520,'log');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.01,m,n,4,520,'log');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,4,520,'log');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,1,m,n,4,520,'log');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,10,m,n,4,520,'log');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,1,520,'log'); %unweighted l1, epsilon value irrelevant here
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.001,m,n,4,520,'atan');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.01,m,n,4,520,'atan');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,4,520,'atan');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,1,m,n,4,520,'atan');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,10,m,n,4,520,'atan');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,1,520,'atan'); %unweighted l1, epsilon value irrelevant here
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.001,m,n,4,520,'tanh');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.01,m,n,4,520,'tanh');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,4,520,'tanh');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,1,m,n,4,520,'tanh');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,10,m,n,4,520,'tanh');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,1,520,'tanh'); %unweighted l1, epsilon value irrelevant here
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.001,m,n,4,520,'sigmoid');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.01,m,n,4,520,'sigmoid');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,4,520,'sigmoid');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,1,m,n,4,520,'sigmoid');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,10,m,n,4,520,'sigmoid');
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,1,520,'sigmoid'); %unweighted l1, epsilon value irrelevant here
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.1,m,n,2,520,'log'); 
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,8,520,'log');
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.1,m,n,2,520,'atan'); 
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,8,520,'atan');
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.1,m,n,2,520,'tanh'); 
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,8,520,'tanh');
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;
gen_array(k_arr,0.1,m,n,2,520,'sigmoid'); 
time_taken = cputime - t;
disp(time_taken);
t = cputime;
gen_array(k_arr,0.1,m,n,8,520,'sigmoid');
time_taken = cputime - t;
disp(time_taken);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function prob_arr = gen_array(k_arr,eps,m,n,iter,trials,mode)
    l = size(k_arr,2);
    prob_arr = zeros(l,1,'double');
    for i = 1:l
        disp(i)
        prob_arr(i) = rp(k_arr(i),eps,m,n,iter,trials/l,mode);
    end
    a1 = num2str(eps);
    a2 = num2str(iter);
    a3 = num2str(trials);
    a4 = num2str(mode);
    a = append(a1, '-', a2, '-', a3, '-', a4, '.mat');
    save(a,'prob_arr','-mat');
end


function recv_prob = rp(k,eps,m,n,iter,trials,mode)
    succ = 0;
    for i = 1:trials
        x = init_vector(zeros(n,1,'double'),k);
        x_r = expt_run(x, eps, iter, m, n, mode);
        d = max(abs(x-x_r));
        if d <= 0.001
            succ = succ + 1;
        end 
    end
    recv_prob = succ/trials;
end




function x_recons = expt_run (x, eps, iter, m, n, mode)
    Phi = randn(m, n);
    x_recons = zeros(n,1,'double');
    y = Phi*x;
    noise_vector=randn(m,1);
    b = 0.2;
    y = y + noise_vector*b*norm(y)/norm(noise_vector);
    W = eye(n);
    for i = 1:iter
        l = min_x(y,Phi,W);
        x_recons = l;
        W = diag(invert(x_recons,eps,n,mode));
    end
end

function wt_arr = invert(x, eps, n, mode)
    wt_arr = zeros(n,1,'double');
    if mode == "log"
        wt_arr = 1./(abs(x) + eps);
    elseif mode == "atan"
        wt_arr = 1./(x.^2 + eps^2);    
    elseif mode == "tanh"
        wt_arr = (sech(x/eps)).^2;
    elseif mode == "sigmoid"
        z = sgm(x/eps);
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
