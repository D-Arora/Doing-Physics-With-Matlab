% mpDecay.m

N0 = 1e5;  % initial number of atoms
lambda1 = 0.0001; %decays/s
lambda2 = 0.00005;
dt = 100; % time step (in seconds)
M = 1500; % total number of time steps
t = 0:dt:(M-1)*dt; % time series
p1 = lambda1*dt;  % probability of having a decay in the time dt
p2 = lambda2*dt;
% Define the array that will be populated with the number of atoms at each
% time step
N1 = zeros(M,1); N2 = zeros(M,1); N3 = zeros(M, 1); 
N1(1) = N0;
N_1 = zeros(M,1); N_2 = zeros(M,1); N_3 = zeros(M,1);
N_1(2)=N0;
N_1 = N0*exp(-lambda1*(t));
N_2 = N0*(lambda1/(lambda1-lambda2))*(exp(-lambda2*t)-exp(-lambda1*t));
N_3 = ((N0)/(lambda1-lambda2))*(lambda2*exp(-lambda1*(t))-lambda1*exp(-lambda2*(t)+(lambda1-lambda2)));
% for k=2:M
% N_1(k) = N_1(k-1);
% N_2(k) = N_2(k-1) + n_decaysN1toN2 - n_decaysN2toN3;
% N_3(k) = N_3(k-1) + n_decaysN2toN3;
% N_1 = N0*exp(-lambda1*(t));
% N_2 = N0*(lambda1/(lambda1-lambda2))*(exp(-lambda2*t)-exp(-lambda1*t));
% N_3 = ((N0)/(lambda1-lambda2))*(lambda2*exp(-lambda1*(t))-lambda1*exp(-lambda2*(t)+(lambda1-lambda2)));
% end
% figure(2)
% plot(t, [N_1, N_2, N_3]), box off
% legend({'# atoms of type 1', '# atoms of type 2', '#atoms of type 3'})
% title('Analytical Calculation of decay')
% xlabel('time')
% ylabel('number of particles')
% legend boxoff
% tic
for j=2:M
    n_decaysN1toN2 = nnz(rand(N1(j-1),1) <= p1);
    n_decaysN2toN3 = nnz(rand(N2(j-1),1) <= p2);
      N1(j) = N1(j-1) - n_decaysN1toN2; % number of type 1 atoms that are left
      N2(j) = N2(j-1) + n_decaysN1toN2 - n_decaysN2toN3;
      N3(j) = N3(j-1) + n_decaysN2toN3;
end

figure(1)
plot(t, [N1, N2, N3]), box off
legend({'# atoms of type 1', '# atoms of type 2', '# atoms of type 3'})
title('Monte Carlo Calculation of Decay')
xlabel('time (s)')
ylabel('number of particles')
legend boxoff