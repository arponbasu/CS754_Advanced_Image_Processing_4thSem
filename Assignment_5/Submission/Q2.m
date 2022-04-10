rng('default')

U = randn(128,128,"double");
U = orth(U.').';


m = [40, 50, 64, 80, 100, 120];
alpha = [0, 3];
lm = length(m);
la = length(alpha);
rmse = zeros(lm, la);
for i = 1:lm
    for j = 1:la
        rmse(i,j) = log(gen_err(m(i),U,alpha(j))); % log is taken for better visual representation
    end
    
end




figure(1)
plot(m, rmse,'-o')                                           
grid
ylabel('log(RMSE)');
xlabel('# of measurements(m)')
legend('alpha=0', 'alpha=3');
savefig('plot.fig')

function xMAP = MAP(Phi, sigma, Sigmax, y)
    xMAP = (inv(Phi'*Phi+sigma^2*inv(Sigmax)))*Phi'*y;
end

function err = gen_err (m, U, alpha)
    Phi = genPhi(m, 128);
    err = 0;
    Sigmax = gen_Sigma(U, alpha, 128);
    x_vec = gen_x(Sigmax);
    for i = 1:10
        x = x_vec(i,:).';
        sigma = calc_sigma(Phi, x);
        y = comp_meas(Phi, x, m, sigma);
        xMAP = MAP(Phi, sigma, Sigmax, y);
        err = err + norm(x-xMAP)/norm(x);
    end
    err = err/10;
end

function sigma = calc_sigma(Phi, x)
    sigma = 0.01*mean2(abs(Phi*x));
end

function x_vec = gen_x (Sigmax)
    sz = size(Sigmax, 1);
    mu = zeros(sz, 1);
    x_vec = mvnrnd(mu,Sigmax,10);

end

function y = comp_meas(Phi, x, m, sigma)
    eta = sigma*randn(m, 1);
    y = Phi*x + eta;
end

function Phi = genPhi (m, n)
    Phi = sqrt(1/m)*randn(m, n);
end

function Sigmax = gen_Sigma(U, alpha, n)
    Lambda = diag(gen_vec(alpha,n));
    Sigmax = U*Lambda*(U.');
end

function v = gen_vec(alpha, n)
    v = zeros(n,1);
    for i = 1:n
        v(i) = i^(-alpha);
    end
end
