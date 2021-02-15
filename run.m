tspan = [0 2];
y0 = [2 4 6];
[t,y] = ode45(@(t,y) syst(t,y), tspan, y0);


function [ret] = expo(t,y)
    ret = -0.5 * y;
end

function [ret] = syst(t,y)
    ret = [y(2)+y(3), -y(1)+2*y(2)+y(3), y(1)+y(3)]';
end

