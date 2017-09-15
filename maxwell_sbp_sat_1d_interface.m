clear all, close all;

ordning_v =  [2 4 6 10];
grd_pts_v = [31 62 124 248 496]; % number of grid points

e1_4 = [1 0 0 0];
e2_4 = [0 1 0 0];
e3_4 = [0 0 1 0];
e4_4 = [0 0 0 1];
% e1 = [1 0];
% em = [0 1];

rr = 0.1; % width of Gaussian
xl = -1;  % left boundary
xr = 1;   % right boundary
x_0 = -1; % initial position of Gaussian
L = xr-xl;  % lenngth of interval
A = kron(eye(2), [0 1; 1 0]); % system matrix ~
t_start = 0; % start time
t_end = 1.8; % end time




% penalty parameters
taul = [-1; 1; 0; 0];
taur = [0; 0; -1; -1];
sigmar = [-1; 1; 1; -1];
sigmal = [-1; -1;-1; 0];


ordning = ordning_v(2);
grd_pts = grd_pts_v(1);
t = t_start;
h = L / (grd_pts - 1);
dt = 0.1*h;
m = grd_pts;
n_steps = floor(t_end/dt);
x = linspace(xl, xr, grd_pts);
Val_operator_ANM;



% SBP - SAT approximation for interfaced system
PP_l = kron(A, D1) + kron(taul, HI*e_1)*kron(e1_4, e_1')...
    + kron(sigmar, HI*e_m)*(kron(e1_4, e_m')- kron(e3_4, e_1'))...
    + kron(sigmar, HI*e_m)*(kron(e2_4, e_m')- kron(e4_4, e_1'));
%+ kron(taur, HI*e_m)*kron(e1_4, e_m')...

PP_r = ...
    + kron(taur, HI*e_m)*kron(e3_4,e_m')...
    + kron(sigmal, HI*e_1)*(kron(e3_4, e_1')- kron(e1_4, e_m'))...
    + kron(sigmar, HI*e_1)*(kron(e4_4, e_1')- kron(e2_4, e_m'));

%kron(A, D1)
%+ kron(taul, HI*e_1)*kron(e1_4, e_1')...

PP = PP_l + PP_r;
P = sparse(PP);
%V = zeros(4*m, 1);
V = [init_cond(x,t_start,rr,x_0); init_cond(x,t_start,rr,x_0);zeros(m, 1);...
    zeros(m, 1)];
temp = zeros(4*m, 1);
w1 = zeros(4*m, 1);
w2 = zeros(4*m, 1);
w3 = zeros(4*m, 1);
w4 = zeros(4*m, 1);

figure(1)
plot(real(eig(PP)), imag(eig(PP)), '*')

for k = 1:1:9
    
    w1 = P*V;
    temp = V + w1/2;
    
    w2 = P*temp;
    temp = V + w2/2;
    
    w3 = P*temp;
    temp = V + w3;
    
    w4 = P*temp;
    
    V = V + (w1 + 2*w2 + 2*w3 + w4)/6;
    t = t + dt;
    %     figure(2)
    %     plot(x, V(1:m), 'r', x, V(m+1:2*m), 'g')%, x, V(2*m+1:3*m), x, V(3*m+1:4*m))
    %     xlim([-1 1]);
    %     ylim([-2 2]);
    %     pause(0.01)
    %     hold
end




function ic = init_cond(x, t, rr, x_0)
ic = exp(-((x-t-x_0)/rr).^2)';
end
