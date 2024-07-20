clearvars;clc;
% Similarity solutions attained from function sim_sol_Crane
% specify viscosity sensitivity parameter & Prandtl number

%%% parameters
% sensitivity parameter
m = -0.4;
% Prandtl number
Pr = 0.72;
% Eckert number
Ec = 0.01;
%%%%%%%%%

% velocity of the sheet
vel =@(x) x;

% define mesh
L = 20; %  stream wise length
H = 40; % spanwise length
nx = 100*L; % Density
ny = 15*H;

dx = L/nx;
dy = H/ny;

% tranformation to concentrate spanwise coordinate near the sheet
y = [0:dy:H]; yy = y;
y = y.*exp(1/y(round(end/2))*(y - y(end)));
dy = [0,diff(y)];

% meshgrid
[X,Y] = meshgrid(0:dx:L,y);


% For plotting
x = X(1,:);

% initialise variables
a = zeros(size(X));
b = a; c = a; d = a; e =a;

% parameters for similarity solution
tspan = [Y(1,1), Y(end,1)];

% initial guess for f''(0) & T'(0)
y0 = [0.2,-.4];

% computes similarity solution at a given x value
[~,bcs] = sim_sol_crane(m,Ec, y,y0,x(1),Pr);

% use as Left BC @ x = dx
a(:,1) = -bcs(:,1); % v = - f
b(:,1) = bcs(:,2)*x(1); % x * f'
c(:,1) = bcs(:,3)*x(1); % x * f''
d(:,1) = bcs(:,4); % T
e(:,1) = bcs(:,5); % T_y

% function handles to simplify code
mid = @(a,j,i) 0.25*(a(j,i) + a(j-1,i) + a(j,i-1) + a(j-1,i-1) );
del_x = @(a,j,i) 1/2/dx*( a(j,i) - a(j,i-1)+ a(j-1,i) - a(j-1,i-1) );
del_y = @(a,j,i) 1/2/dy(j)*( a(j,i) + a(j,i-1) - a(j-1,i) - a(j-1,i-1) );
j12 = @(a,j,i) 0.5*(a(j,i) + a(j-1,i));

mid_mu =@(d,j,i) 1/(1 + m*mid(d,j,i));

% Errors associated with current step
r1 = @ (a,b,j,i)  -del_y(a,j,i) - del_x(b,j,i) ;
r2 =  0; % prevents oscillations from developing in corrections
r3 = @(a,b,c,d,e,j,i) -mid_mu(d,j,i)*(del_y(c,j,i) - m*mid_mu(d,j,i)*mid(c,j,i)*mid(e,j,i) ) ...
    + mid(b,j,i)*del_x(b,j,i)  + mid(a,j,i)*mid(c,j,i) ;
r4 = 0; % prevents oscillations from developing in corrections
r5 = @(a,b,c,d,e,j,i) -1/Pr*del_y(e,j,i) + mid(b,j,i)*del_x(d,j,i) ...
    + mid(a,j,i)*mid(e,j,i) - Ec*mid_mu(d,j,i)*(mid(c,j,i)).^2;

% ce coefficients
% a
alpha1 = @(j) 0.5/dy(j); % delta a_{j+1}
alpha2 = @(j) -0.5/dy(j); % delta a_j
% b
alpha3 = 0.5/dx;
alpha4 = alpha3;

% momentum eq
% a
beta1 = @(c,j,i) -0.25*mid(c,j,i); % delta a_{j+1}
beta2 = beta1;                     % delta a_j and so on
% b
beta3 = @(b,j,i) - 0.25*del_x(b,j,i) - 0.5/dx*mid(b,j,i);
beta4 = beta3;
% c
beta5 =@(a,d,e,j,i) mid_mu(d,j,i)*(1/2/dy(j) - 0.25*m*mid_mu(d,j,i)*mid(e,j,i)) -0.25*mid(a,j,i);
beta6 = @(a,d,e,j,i) mid_mu(d,j,i)*(-1/2/dy(j) - 0.25*m*mid_mu(d,j,i)*mid(e,j,i)) -0.25*mid(a,j,i);
% d
beta7 = @(c,d,e,j,i) -0.25*m*( mid_mu(d,j,i) )^2*(del_y(c,j,i) - 2*m*mid_mu(d,j,i)*mid(e,j,i)*mid(c,j,i));
beta8 =  beta7;
% e
beta9 = @(c,d,j,i) -0.25*(mid_mu(d,j,i))^2*m*mid(c,j,i);
beta10 = beta9;

% energy eq
% a
gamma1 = @(e,j,i) - 0.25*mid(e,j,i);
gamma2  = gamma1;
% b
gamma3 = @(d,j,i) - 0.25*del_x(d,j,i);
gamma4 = gamma3;
% c
gamma5 =@ (c,d,j,i) Ec/2*mid_mu(d,j,i)*mid(c,j,i);
gamma6 = gamma5;
% d
gamma7 = @(b,c,d,j,i) -1/2/dx*mid(b,j,i) - 0.25*Ec*m*( mid_mu(d,j,i) )^2*(mid(c,j,i)).^2;
gamma8= gamma7;
% e
gamma9 = @(a,d,j,i)   1/Pr/2/dy(j) -0.25*mid(a,j,i);
gamma10 = @(a,d,j,i) -1/Pr/2/dy(j) -0.25*mid(a,j,i);

% initialise matrices
A = zeros(5,5,length(y));
B = A; C = A;
R = zeros(5,1, length(y));

% keller box loop

for i = 2:length(x)
    % Initialise next step
    % set previous solution as initial guess
    a(:,i) = a(:,i-1); % v
    b(:,i) = b(:,i-1); % u
    c(:,i) = c(:,i-1);% u_Y
    d(:,i) = d(:,i-1); % T
    e(:,i) = e(:,i-1); % T_y

    count = 0; H = 1;


    % H = max | correction |
    % exits loop if corrections grow too large
    while H > 1e-10 && H <= 10

        for j = 1:length(y)
            if j == 1
                B(:,:,j) = [1, 0, 0, 0, 0;
                    0, 1, 0, 0, 0;
                    0, 0, 0, 1, 0;
                    0,0, 0, -1/dy(j+1), -1/2;
                    0, -1/dy(j+1), -1/2, 0, 0];
                C(:,:,j) = [0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0;
                    0,0, 0, 1/dy(j+1), -1/2;
                    0, 1/dy(j+1), -1/2, 0, 0];
                % BC's @ y = 0
                % v = 0         u = x          T = 1   @ y = 0
                R(:,:,j) = [-a(j,i); vel(x(i)) - b(j,i) ;1 - d(j,i); 0; 0];
            elseif j >= 2 && j < length(y)
                A(:,:,j) = [beta2(c,j,i), beta4(b,j,i), beta6(a,d,e,j,i), beta8(c,d,e,j,i), beta10(c,d,j,i);
                    alpha2(j), alpha4, 0, 0, 0;
                    gamma2(e,j,i), gamma4(d,j,i), gamma6(c,d,j,i), gamma8(b,c,d,j,i), gamma10(a,d,j,i);
                    0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0] ;
                B(:,:,j) = [beta1(c,j,i), beta3(b,j,i), beta5(a,d,e,j,i), beta7(c,d,e,j,i), beta9(c,d,j,i)
                    alpha1(j), alpha3, 0, 0, 0;
                    gamma1(e,j,i), gamma3(d,j,i), gamma5(c,d,j,i), gamma7(b,c,d,j,i), gamma9(a,d,j,i);
                    0,0, 0, -1/dy(j), -1/2;
                    0, -1/dy(j), -1/2, 0, 0];
                C(:,:,j) = [0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0;
                    0,0, 0, 1/dy(j), -1/2;
                    0, 1/dy(j), -1/2, 0, 0];
                R(:,:,j) = [ r3(a,b,c,d,e,j,i); r1(a,b,j,i); r5(a,b,c,d,e,j,i); 0; 0];
            elseif j == length(y)
                A(:,:,j) = [beta2(c,j,i), beta4(b,j,i), beta6(a,d,e,j,i), beta8(c,d,e,j,i), beta10(c,d,j,i);
                    alpha2(j), alpha4, 0, 0, 0;
                    gamma2(e,j,i), gamma4(d,j,i), gamma6(c,d,j,i), gamma8(b,c,d,j,i), gamma10(a,d,j,i);
                    0, 0, 0, 0, 0;
                    0, 0, 0, 0, 0];
                B(:,:,j) = [beta1(c,j,i), beta3(b,j,i), beta5(a,d,e,j,i), beta7(c,d,e,j,i), beta9(c,d,j,i)
                    alpha1(j), alpha3, 0, 0, 0;
                    gamma1(e,j,i), gamma3(d,j,i), gamma5(c,d,j,i), gamma7(b,c,d,j,i), gamma9(a,d,j,i);
                    0, 1, 0, 0, 0;
                    0, 0, 0, 1 0];
                % BC's @ y = Y_max
                % u = 0    T = 0 as y -> \infty
                R(:,:,j) = [r3(a,b,c,d,e,j,i); r1(a,b,j,i); r5(a,b,c,d,e,j,i);-b(j,i); -d(j,i)];
            end

        end
        % bdtma function implements block version of TDMA algoritm
        % written following Physical and Computational Aspects of
        % Convective Heat Transfer, Cebici & Bradshaw
        correction = reshape(btdma(A,B,C,R), 5,size(a,1));

        % update solution
        a(:,i) = a(:,i) + correction(1,:)'; % v
        b(:,i) = b(:,i) + correction(2,:)'; % u
        c(:,i) = c(:,i) + correction(3,:)'; % u_Y
        d(:,i) = d(:,i) + correction(4,:)'; % T
        e(:,i) = e(:,i) + correction(5,:)'; % T_y

        H = max(max(abs(correction)));

        count  = count +1;

    end
    % Exit loop if solution didn't converge
    if H > 1
        break
    end
end

% Plotting
% compared similarity solution without dissipation to
% the dissipitive regime at different x locations

[~,f_1] = sim_sol_crane(m,Ec, y,[0,0], x(100),Pr);


figure(1), clf
h = plot(y,d(:,100), y,d(:,round(end/2)),y, d(:,end),y, f_1(:,4));
set(h,{'LineStyle'},{'-','--',':','-.'}')
xlabel('$y$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')
legend("$x = 1$","$x = 10$","$x = 20$","$T$", 'Interpreter', 'latex')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment', 'right')
xlim([0,10])

figure(2),clf
h = plot(y,b(:,100)/x(100), y, b(:,round(end/2))/x(round(end/2)) , y, b(:,end)/x(end), y, f_1(:,2) );
set(h,{'LineStyle'},{'-','--',':','-.'}')
xlabel('$y$','Interpreter','latex')
ylabel('$\frac{u_B}{x}$','Interpreter','latex')
legend("$x = 1$","$x = 10$","$x = 20$","$f'$", 'Interpreter', 'latex')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment', 'right')
xlim([0,10])

figure(3), clf
h = plot(y,a(:,100),y,a(:,round(end/2)) ,y, a(:,end),y,- f_1(:,1));
set(h,{'LineStyle'},{'-','--',':','-.'}')
xlabel('$y$','Interpreter','latex')
ylabel('$v$','Interpreter','latex')
legend("$x = 1$","$x = 10$","$x = 20$","$-f$", 'Interpreter', 'latex')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment', 'right')
xlim([0,10])

