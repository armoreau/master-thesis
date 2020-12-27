

%sol = ode45(@eq3, [1, 2] , [0]); %point de solve_ivp

A = 1;
B = 2;
tspan = [0 5];
y0 = [0 0.01];
[t,y] = ode45(@(t,y) odefcn(t,y,A,B), tspan, y0);

dlmwrite('Matlab_t.txt ',t,'delimiter',' ','precision','%.13e')
dlmwrite('Matlab_y.txt ',y,'delimiter',' ','precision','%.13e')

function [ret] = eq1(t,y)
    ret = [y(2)+y(3), -y(1)+2*y(2)+y(3), y(1)+y(3)]';
end

function [ret] = sol1(t)
    ret = [-3+5*exp(2*t), -3+2*exp(t)+5*exp(2*t), 3-2*exp(t)+5*exp(2*t)];
end

function [ret] = eq2(t,y)
    ret = -0.5 * y;
end

function [ret] = sol2(t)
    ret = 2*exp((-1/2)*t);
end

function [ret] = eq3(t,y)
    ret = -cos(1/y)/(y^2);
end

function [ret] = sol3(t)
    ret = sin(1/t);
end

function dydt = odefcn(t,y,A,B)
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = (A/B)*t.*y(1);
end

