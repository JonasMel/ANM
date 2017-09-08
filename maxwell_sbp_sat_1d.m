clear all, close all;

ordning = 4;
grd_pts = 500; % number of grid points
mm = grd_pts*2;
m = grd_pts;
e1 = [1 0];
em = [0 1];
rr = 0.1; % width of Gaussian
xl = -1;  % left boundary
xr = 1;   % right boundary
L = xr-xl;  % lenngth of interval
x = linspace(xl, xr, grd_pts); % discrete x-values


% penalty parameters
taul = [-0.5; 1];
taur = [-0.5; -1];

t_end = 1.8; % end time
t = 0;
h = L / grd_pts; % step size
dt = 0.1*h; % step size in time
n_steps = round(t_end/dt)+1;
Val_operator_ANM; % initializing difference operator and shit..

A = [0 1; 1 0]; % system matrix ~

% SBP-SAT operator with Dirichlet BCs.
PP = kron(A, D1) + kron(taul, HI)*e_1*kron(e1,e_1') + ...
        kron(taur, HI)*e_m*kron(em,e_m');
P = dt*sparse(PP);

% initializing vectors for RK4
tmp = zeros(mm,1);
w1 = zeros(mm,1);
w2 = zeros(mm,1);
w3 = zeros(mm,1);
w4 = zeros(mm,1);


V = zeros(mm,1);
V_ex = zeros(mm,1);
V(1:grd_pts) = -exp(-((x+t)/rr).^2) - exp(-((x-t)/rr).^2);
V(grd_pts+1:mm) = -exp(-((x+t)/rr).^2) + exp(-((x-t)/rr).^2);

V_ex(1:grd_pts) = exp(-((x-(L-t))/rr).^2) + exp(-((x+(L-t))/rr).^2);
V_ex(grd_pts+1:mm) = exp(-((x-(L-t))/rr).^2) - exp(-((x+(L-t))/rr).^2);

err_E = zeros(n_steps,1);
err_H = zeros(n_steps,1);
for k = 1:1:n_steps
    
    w1 = P*V;
    temp = V + w1/2;
    
    w2 = P*temp;
    temp = V + w2/2;
    
    w3 = P*temp;
    temp = V + w3;
    
    w4 = P*temp;

    V = V + (w1 + 2*w2 + 2*w3 + w4)/6;
    t = t + dt;

    V_ex(1:grd_pts) = exp(-((x-(L-t))/rr).^2) + exp(-((x+(L-t))/rr).^2);
    V_ex(grd_pts+1:mm) = exp(-((x-(L-t))/rr).^2) - exp(-((x+(L-t))/rr).^2);
    plot(x, V(1:m), 'r', x, V_ex(1:m), 'g', x, V(m+1:mm), x, V_ex(m+1:mm))
    xlim([-1 1]);
    ylim([-2 2]);
    pause(0.000001)
    hold

    err_E(k) = sqrt(h)*norm(V(1:m) - V_ex(1:m));
    err_H(k) = sqrt(h)*norm(V(m+1:mm) - V_ex(m+1:mm));
end
figure(2)
plot(real(eig(PP)),imag(eig(PP)), '*');