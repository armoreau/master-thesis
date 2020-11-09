vect_temps = [0 ,0.11488132, 1.26369452, 3.06074656, 4.81637262, 6.57504937, 8.33467262, 10];
[T, Y] = ode45(@eq1, [0, 19] , [2,4,6],'RequestedPoints'); %point de solve_ivp

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

