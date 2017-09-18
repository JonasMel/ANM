clear all, close all;

ordning_v =  2;%[2 4 6 10];
grd_pts_v = 400;%[31 62 124 248 496]; % number of grid points

% e1_4 = [1 0 0 0];
% e2_4 = [0 1 0 0];
% e3_4 = [0 0 1 0];
% e4_4 = [0 0 0 1];
e1 = [1 0];
em = [0 1];
rr = 0.1; % width of Gaussian
xl = -1;  % left boundary
xr = 1;   % right boundary
L = xr-xl;  % lenngth of interval
A = [0 1; 1 0]; % system matrix ~
[S, Lambda] = eig(A);
Lambda_pos = (Lambda + abs(Lambda))*0.5;
Lambda_neg = (Lambda - abs(Lambda))*0.5;
A_p = S*Lambda_pos/S;
A_n = S*Lambda_neg/S;
II = eye(2);
I_u = [1 0; 0 0];
I_l = [0 0; 0 1];
% penalty parameters
taul = [0; 1];
taur = [0; -1];

t_end = 1*1.8; % end time
err_E = zeros(length(grd_pts_v),length(ordning_v));
err_H = zeros(length(grd_pts_v),length(ordning_v));
q = zeros(length(ordning_v),1);
h_v = zeros(1, length(grd_pts_v));

for ii = 1:length(ordning_v)
    ordning = ordning_v(ii);
    for jj = 1:22
        grd_pts = grd_pts_v(1);
        t = 0;
        h = L / (grd_pts-1); % step size
        h_v(jj) = h;
        dt = jj*0.1*h; % step size in time
        x = linspace(xl, xr, grd_pts); % discrete x-values
        mm = grd_pts*2;
        m = grd_pts;
        n_steps = floor(t_end/dt);
        Val_operator_ANM; % initializing difference operator and shit..
        
        
        
        % SBP-SAT operator with Dirichlet BCs.
        PP = kron(A, D1) + kron(taul, HI*e_1)*kron(e1,e_1') + ...
            kron(taur, HI*e_m)*kron(e1,e_m');
        P = dt*sparse(PP);
        
        % SBP = -SAT approximation for characteristic
%                 PP = kron(A, D1) + kron(A_n, HI*e_1*e_1') - kron(A_p, HI*e_m*e_m');
%                 P = dt*sparse(PP);
        
%         % SBP - SAT approximation for interfaced system
%         PP_l = kron(A, D1) + kron(taul, HI*e_1)*kron(e1, e_1')...
%                 + kron(taur, HI*e_m)*kron(e1,e_m')...
%                 + kron(sigmar, HI*e_m)*kron(II,e_m');
%         
%         PP_r = kron(A, D1) + kron(taul, HI*e_1)*kron(e1, e_1')...
%                 + kron(taur, HI*e_m)*kron(e1,e_m')...
%                 + kron(sigmal, HI*e_1)*kron(II,e_1');
%         
%         PP = kron(I_u, PP_l) + kron(I_l, PP_r);

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
        
        
        for k = 1:1:7000
            
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
            
            if (mod(k,40) == 0)
                        figure(1)
                        plot(x, V(1:m), 'r', x, V_ex(1:m), 'g', x, V(m+1:mm), x, V_ex(m+1:mm))
                        xlim([-1 1]);
                        ylim([-2 2]);
                        pause(0.000001)
                        hold
            end
            if k == floor(1*n_steps)
%                                 plot(x, V(1:m), 'r', x, V_ex(1:m), 'g', x, V(m+1:mm), x, V_ex(m+1:mm))
%                                 xlim([-1 1]);
%                                 ylim([-2 2]);
%                                 pause(0.000001)
%                                 hold
                 err_E(jj,ii) = sqrt(h)*norm(V(1:m) - V_ex(1:m))
%                 err_H(jj,ii) = sqrt(h)*norm(V(m+1:mm) - V_ex(m+1:mm));
            end
        end
        
    end
%     q(ii,1) = log10(err_E(1,ii)/err_E(end,ii))...
%         /log10(grd_pts_v(1)/grd_pts_v(end));
    
end

% figure(2)
% plot(real(eig(PP)),imag(eig(PP)), '*');
% figure(3)
% plot(ordning_v, q, '^-.')
% xlabel('Stated order of method');
% ylabel('q');
