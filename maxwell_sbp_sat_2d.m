clear all; close all;

ordning_v= 10;% [2 4 6 10];
grd_pts_v = 70; %[20 30 40 50 60 70];


e1 = [1 0 0];
e2 = [0 1 0];
e3 = [0 0 1];

rr = 0.1;
xl = -1;
xr = 1;
x0 = 0;
yl = -1;
yr = 1;
y0 = 0;
L = xr - xl;
area = L^2;

t_start = 0;
t_end = 2; %1.8;


A = [0 0 0; 0 0 -1; 0 -1 0];
B = [0 1 0; 1 0 0; 0 0 0];
tauw = [0; 0; -2];
taue = [0; 0; 2];
taus = [2; 0; 0];
taun = [-2; 0; 0];


q = zeros(length(ordning_v),1);
theta_n  = zeros(length(grd_pts_v),1);
for j = 1 : length(ordning_v)
    ordning = ordning_v(j);
    for i = 1:length(grd_pts_v)
        m = grd_pts_v(i);
        x = linspace(-1, 1, m);
        y = linspace(-1, 1, m);
        h = L / (m-1);
        dt = 0.1*h;
        n_steps = floor(t_end/dt);

        Val_operator_ANM;
        I_m = eye(m);
        %     Dx = kron(D1, I_m);
        %     Dy = kron(I_m, D1);
        Hx = kron(H, I_m);
        Hy = kron(I_m, H);
        %     HIx = inv(Hx);
        %     HIy = inv(Hy);
        H = Hx*Hy;
        %     ew = kron(e_1, I_m);
        %     ee = kron(e_m, I_m);
        %     es = kron(I_m, e_1);
        %     en = kron(I_m, e_m);
        
        
        %     satw = kron(tauw, HIx)*ew*kron(e2, ew');
        %     sate = kron(taue, HIx)*ee*kron(e2, ee');
        %     sats = kron(taus, HIy)*es*kron(e2, es');
        %     satn = kron(taun, HIy)*en*kron(e2, en');
        PP = kron(A, kron(D1, I_m)) + kron(B, kron(I_m, D1))...
            + kron(tauw, inv(Hx))*kron(e_1, I_m)*kron(e2, kron(e_1, I_m)')...
            + kron(taue, inv(Hx))*kron(e_m, I_m)*kron(e2, kron(e_m, I_m)')...
            + kron(taus, inv(Hy))*kron(I_m, e_1)*kron(e2, kron(I_m, e_1)')...
            + kron(taun, inv(Hy))*kron(I_m, e_m)*kron(e2, kron(I_m, e_m)');
        
        %     PP = kron(A, Dx) + kron(B, Dy) + satw + sate + sats + satn;
        P = dt*sparse(PP);
        
        V = [zeros(m*m,1); init_cond1(x, x0, y, y0, rr); zeros(m*m,1)];
        
        temp = zeros(3*m*m, 1);    % Temporary vector in RK4
        w1 = zeros(3*m*m, 1);      % Step 1 vector in RK4
        w2 = zeros(3*m*m, 1);      % Step 2 vector in RK4
        w3 = zeros(3*m*m, 1);      % Step 3 vector in RK4
        w4 = zeros(3*m*m, 1);      % Step 4 vector in RK4
        
        t = t_start;
        for k = 1:n_steps
            w1 = P*V;
            temp = V + w1*0.5;
            
            w2 = P*temp;
            temp = V + w2*0.5;
            
            w3 = P*temp;
            temp = V + w3;
            
            w4 = P*temp;
            
            V = V + (w1 + 2*w2 + 2*w3 + w4)/6;
            Vmat_1 = vec2mat(V(1:m*m), m);
            Vmat_2 = vec2mat(V(m*m+1:2*m*m), m);
            Vmat_3 = vec2mat(V(2*m*m+1:3*m*m), m);
            
            
            
            t = t + dt;
            % Plot
            if mod(k,15) == 0
                s = surf(x, y, Vmat_2);
                s.EdgeColor = 'none';
                xlim([-1 1])
                ylim([-1 1])
                zlim([-1 1])
                pause(0.001)
            end
        end
        theta_n(i) = sqrt((kron(D1, I_m)*V(1:m*m)...
            +kron(I_m, D1)*V(2*m*m+1:3*m*m))'*H*(kron(D1, I_m)*V(1:m*m)...
            +kron(I_m, D1)*V(2*m*m+1:3*m*m)));
    end
    q(j) = log10(theta_n(1)/theta_n(end))...
        /log10(grd_pts_v(end)/grd_pts_v(1));
end


plot(ordning_v, q);

function ic = init_cond1(x, x_0, y, y_0, rr)
mm = length(x);
nn = length(y);
ic = zeros( mm*nn , 1);
for i = 0:(mm-1)
    for j = 1:nn
        ic(i*mm + j, 1) = exp(-((x(i+1)-x_0)/rr).^2 - ((y(j)-y_0)/rr).^2);
    end
end
end

