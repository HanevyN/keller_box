function [t,y] = sim_sol_crane(m,Ec, y,y0,x,Pr)
% m - viscocity parameter
% Ec - Eckert Number
% y0 - guessed initial conditions
% x - x location (for Ec != 0 )
% Pr = prandtl number 
% y - mesh you want your final solution on 

a = x;

% function handle for Cranes ODE

NRF = @(x,y)[y(2); 
             y(3); 
           ( y(2).^2 - y(1)*y(3)  + m.*(1 + m.*y(4) ).^(-2)*y(3).*y(5) ).*(1 + m.*y(4) );
             y(5) ; 
       -Pr*( y(1)*y(5) + Ec/( 1 + m*y(4) )*(a*y(3)).^2 ) ];

tspan = 2:y(end);
ICS = y0;

% change  in guessed I.C
dp = 0.000001;

K = zeros(2);

% Newton Loop
% Allow size of the  domain to increase improve probability of convergence
    for i = 1:length(tspan)
        H = [1,1];
        count = 0;
        while  max(abs(H)) > 1e-10  
            count = count+ 1;              %   f f' f''         T  T'
            [~,y1] = ode45(NRF, [0,tspan(i)], [0;1;ICS(1);      1; ICS(2)     ] );
            [~,y2] = ode45(NRF, [0,tspan(i)], [0;1;ICS(1) - dp; 1; ICS(2)     ] );
            [~,y3] = ode45(NRF, [0,tspan(i)], [0;1;ICS(1);      1; ICS(2) - dp] );
           
            
            K(1,1) = (y2(end,2) - y1(end,2))/dp; 
            K(2,1) = (y2(end,4) - y1(end,4))/dp;
            K(1,2) = (y3(end,2) - y1(end,2))/dp;
            K(2,2) = (y3(end,4) - y1(end,4))/dp;
            
            H = K\[y1(end,2); y1(end,4)];

            ICS = ICS + H';
            
        end
     
    end

% make sure output is on correct y mesh    
[t,y] = ode45(NRF, y, [0;1;ICS(1); 1; ICS(2) ]); 

end