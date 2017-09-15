clear all; close all;

ordning_v = [2 4 6 10];
grd_pts_v = [31 62 124 248 496];




e1 = [1 0];
e2 = [0 1];

t_start = 0;
t_end = 3.6;

x_l = -1;
x_r = 1;
x_0 = -0.5;
L = x_r - x_l;
rr = 0.1;

eta_l = 2;
eta_r = 1;
T = 2*eta_l/(eta_l+eta_r);
R = (eta_l - eta_r)/(eta_l + eta_r);
A = [0 1; 1 0];
[S, lambda] = eig(A);
A_p = S*((lambda + abs(lambda))*0.5)/S;
A_m = S*((lambda - abs(lambda))*0.5)/S;
C_l = [eta_l 0; 0 1];
C_r = [eta_r 0; 0 1];
I2 = eye(2);

sigmal = 1;
sigmar = -1;
taul = [-1; 1];
taur = [-1; -1];

t = t_start;
ordning = ordning_v(2);
m = grd_pts_v(3);
h = L/(2*m - 1);
dt = h*0.1;
n_steps = floor(t_end/dt);

%x = linspace(x_l, x_r, m);
x_L = linspace(x_l, 0, m);
x_R = linspace(0, x_r, m);

Val_operator_ANM;

Ap_Cl = C_l\A_p;
Am_Cr = C_r\A_m;

PP_l = kron(C_l\A, D1) + sigmar*kron(Ap_Cl, HI*e_m*e_m')...
    + kron(C_l\taul*e1, HI*e_1*e_1');

PP_r = kron(C_r\A, D1) + sigmal*kron(Am_Cr, HI*e_1*e_1')...
    + kron(C_r\taur*e1, HI*e_m*e_m');



P_l = dt*sparse(PP_l);
P_r = dt*sparse(PP_r);

IT_l = -dt*sigmar*sparse(kron(Ap_Cl, HI)*kron(I2, e_m*e_1'));
IT_r = -dt*sigmal*sparse(kron(Am_Cr, HI)*kron(I2, e_1*e_m'));

% IT_l = dt*sigmar*sparse(-kron(Ap_Cl, HI))*kron(eye(2), e_m*e_1');    % Penalty data
% IT_r = dt*sigmal*sparse(-kron(Am_Cr, HI))*kron(eye(2), e_1*e_m');

V_l = [theta_2(x_L, x_0, 0, rr) - theta_1(x_L, x_0, 0, rr);...
        theta_2(x_L, x_0, 0, rr) + theta_1(x_L, x_0, 0, rr)];
V_r = zeros(2*m,1);

% Pre-allocate RK-vectors and errors
templ = zeros(2*m, 1);    % Temporary vector in RK4
w1l = zeros(2*m, 1);      % Step 1 vector in RK4
w2l = zeros(2*m, 1);      % Step 2 vector in RK4
w3l = zeros(2*m, 1);      % Step 3 vector in RK4
w4l = zeros(2*m, 1);      % Step 4 vector in RK4
tempr = zeros(2*m, 1);    % Temporary vector in RK4
w1r = zeros(2*m, 1);      % Step 1 vector in RK4
w2r = zeros(2*m, 1);      % Step 2 vector in RK4
w3r = zeros(2*m, 1);      % Step 3 vector in RK4
w4r = zeros(2*m, 1);      % Step 4 vector in RK4

% Setup video
vidObj = VideoWriter('Test_1D.avi');
open(vidObj);
update_movie = 10; % How often to update movie

for k = 1:n_steps
    
    w1l = P_l*V_l + IT_l*V_r;
    w1r = P_r*V_r + IT_r*V_l;
    templ = V_l + w1l*0.5;
    tempr = V_r + w1r*0.5;
    
    w2l = P_l*templ + IT_l*tempr;
    w2r = P_r*tempr + IT_r*templ;
    templ = V_l + w2l*0.5;
    tempr = V_r + w2r*0.5;
    
    w3l = P_l*templ + IT_l*tempr;
    w3r = P_r*tempr + IT_r*templ;
    templ = V_l + w3l;
    tempr = V_r + w3r;
    
    w4l = P_l*templ + IT_l*tempr;
    w4r = P_r*tempr + IT_r*templ;
    
    V_l = V_l + (w1l + 2*w2l + 2*w3l + w4l)/6;
    V_r = V_r + (w1r + 2*w2r + 2*w3r + w4r)/6;
    %V_l(1)
    
    t = t+dt;
    
    %PLOTTA SKITEN
    if mod(k, update_movie) == 0
        plot(x_L, V_l(1:m), 'b', x_L, V_l(m+1:end), 'r', ...
            x_R, V_r(1:m), 'b', x_R, V_r(m+1:end), 'r')
        ylim([-2 2]);
        xlim([-1 1])
        currFrame = getframe;
        writeVideo(vidObj, currFrame);
    end
    
end
close(vidObj)



% function ic = init_cond1(x, t, rr, x_0)
% ic = exp(-((x-t-x_0)/rr).^2)';
% end
%
% function ic = init_cond2(x, t, rr, x_0)
% ic = -exp(-((x+t-x_0)/rr).^2)';
% end


% Help functions
function theta = theta_1(x, x0, t, rr)
theta = exp(-((x - x0 - t)/rr).^2)';
end

function theta = theta_2(x, x0, t, rr)
theta = -exp(-((x - x0 + t)/rr).^2)';
end










