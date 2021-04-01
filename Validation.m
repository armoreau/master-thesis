%% Test 1
tspan = linspace(0,10,25);
y0 = [2];
[t,y] = ode45(@fun1, tspan, y0);

dlmwrite('Test1_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test1_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 2
tspan = linspace(0,10,25);
y0 = [0];
[t,y] = ode45(@fun2, tspan, y0);

dlmwrite('Test2_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test2_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 3
tspan = linspace(0,2,25);
y0 = [2 4 6];
[t,y] = ode45(@fun3, tspan, y0);

dlmwrite('Test3_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test3_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 4
tspan = linspace(0,20,25);
y0 = [2 0];
[t,y] = ode45(@fun4, tspan, y0);

dlmwrite('Test4_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test4_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 5
A = 1;
B = 2;
tspan = linspace(0,5,25);
y0 = [0 0.01];
[t,y] = ode45(@(t,y) fun5(t,y,A,B), tspan, y0);

dlmwrite('Test5_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test5_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 6
tspan = linspace(0,3,25);
y0 = [-3 -2 -1 0 1 2 3];
[t,y] = ode45(@fun6, tspan, y0);

dlmwrite('Test6_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test6_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 7
tspan = linspace(0,10,25);
y0 = [2];

opts = odeset('Refine',10,'NormControl',true,'MaxStep',1,'InitialStep',0.1);
[t,y] = ode45(@fun1, tspan, y0, opts);

dlmwrite('Test7_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test7_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')


%% Test 8
A = 1;
B = 2;
tspan = linspace(0,5,25);
y0 = [0 0.01];

opts = odeset('Refine',1,'NormControl',false,'MaxStep',0.5,'InitialStep',0.1);
[t,y] = ode45(@(t,y) fun5(t,y,A,B), tspan, y0, opts);

dlmwrite('Test8_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test8_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 9
tspan = linspace(0,3,25);
y0 = [0 0];

[t,y] = ode45(@fun9,tspan,y0);

dlmwrite('Test9_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test9_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 10
tspan = linspace(0,4,25);
y0 = [0; 4; 2; 20; -pi/2; 2];
options = odeset('Mass',@mass_fun8);

[t, y] = ode45(@fun8,tspan,y0,options);

dlmwrite('Test10_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test10_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 11
tspan = linspace(10,0,25);
y0 = [2];
[t,y] = ode45(@fun1, tspan, y0);

dlmwrite('Test11_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test11_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 12
tspan = linspace(10,0,25);
y0 = [0];
[t,y] = ode45(@fun2, tspan, y0);

dlmwrite('Test12_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test12_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 13
tspan = linspace(2,0,25);
y0 = [2 4 6];
[t,y] = ode45(@fun3, tspan, y0);

dlmwrite('Test13_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test13_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 14
tspan = linspace(20,0,25);
y0 = [2 0];
[t,y] = ode45(@fun4, tspan, y0);

dlmwrite('Test14_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test14_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 15
A = 1;
B = 2;
tspan = linspace(5,0,25);
y0 = [0 0.01];
[t,y] = ode45(@(t,y) fun5(t,y,A,B), tspan, y0);

dlmwrite('Test15_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test15_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 16
tspan = linspace(3,0,25);
y0 = [-3 -2 -1 0 1 2 3];
[t,y] = ode45(@fun6, tspan, y0);

dlmwrite('Test16_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test16_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 17
tspan = linspace(10,0,25);
y0 = [2];

opts = odeset('Refine',10,'NormControl',true,'MaxStep',1,'InitialStep',0.1);
[t,y] = ode45(@fun1, tspan, y0, opts);

dlmwrite('Test17_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test17_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')


%% Test 18
A = 1;
B = 2;
tspan = linspace(5,0,25);
y0 = [0 0.01];

opts = odeset('Refine',1,'NormControl',false,'MaxStep',0.5,'InitialStep',0.1);
[t,y] = ode45(@(t,y) fun5(t,y,A,B), tspan, y0, opts);

dlmwrite('Test18_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test18_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 19
tspan = linspace(3,0,25);
y0 = [0 0];

[t,y] = ode45(@fun9,tspan,y0);

dlmwrite('Test19_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test19_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 20
tspan = linspace(4,0,25);
y0 = [0, 4, 2, 20, -pi/2, 2];
options = odeset('Mass',@mass_fun8);

[t, y] = ode45(@fun8,tspan,y0,options);

dlmwrite('Test20_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test20_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 21
tspan = linspace(0,10,25);
y0 = [2];
opts = odeset('RelTol',1e-2,'AbsTol',1e-3);
[t,y] = ode45(@fun1, tspan, y0, opts);

dlmwrite('Test21_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test21_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 22
tspan = linspace(0,10,25);
y0 = [0];
opts = odeset('RelTol',1e-2,'AbsTol',1e-3);
[t,y] = ode45(@fun2, tspan, y0, opts);

dlmwrite('Test22_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test22_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 23
tspan = linspace(0,2,25);
y0 = [2 4 6];
opts = odeset('RelTol',1e-2,'AbsTol',1e-3);
[t,y] = ode45(@fun3, tspan, y0, opts);

dlmwrite('Test23_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test23_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 24
tspan = linspace(0,20,25);
y0 = [2 0];
opts = odeset('RelTol',1e-2,'AbsTol',1e-3);
[t,y] = ode45(@fun4, tspan, y0, opts);

dlmwrite('Test24_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test24_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 25
A = 1;
B = 2;
tspan = linspace(0,5,25);
y0 = [0 0.01];
opts = odeset('RelTol',1e-2,'AbsTol',1e-3);
[t,y] = ode45(@(t,y) fun5(t,y,A,B), tspan, y0, opts);

dlmwrite('Test25_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test25_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 26
tspan = linspace(0,3,25);
y0 = [-3 -2 -1 0 1 2 3];
opts = odeset('RelTol',1e-2,'AbsTol',1e-3);
[t,y] = ode45(@fun6, tspan, y0, opts);

dlmwrite('Test26_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test26_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 27
tspan = linspace(0,10,25);
y0 = [2];

opts = odeset('RelTol',1e-2,'AbsTol',1e-3,'Refine',10,'NormControl',true,'MaxStep',1,'InitialStep',0.1);
[t,y] = ode45(@fun1, tspan, y0, opts);

dlmwrite('Test27_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test27_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')


%% Test 28
A = 1;
B = 2;
tspan = linspace(0,5,25);
y0 = [0 0.01];

opts = odeset('RelTol',1e-2,'AbsTol',1e-3,'Refine',1,'NormControl',false,'MaxStep',0.5,'InitialStep',0.1);
[t,y] = ode45(@(t,y) fun5(t,y,A,B), tspan, y0, opts);

dlmwrite('Test28_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test28_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 29

tspan = linspace(0,3,25);
y0 = [0 0];
opts = odeset('RelTol',1e-2,'AbsTol',1e-3);
[t,y] = ode45(@fun9,tspan,y0,opts);

dlmwrite('Test29_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test29_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 30
tspan = linspace(0,4,25);
y0 = [0; 4; 2; 20; -pi/2; 2];
opts = odeset('RelTol',1e-2,'AbsTol',1e-3,'Mass',@mass_fun8);
[t, y] = ode45(@fun8,tspan,y0, opts);

dlmwrite('Test30_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test30_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 31
tspan = linspace(0,10,25);
y0 = [2];
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t,y] = ode45(@fun1, tspan, y0, opts);

dlmwrite('Test31_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test31_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 32
tspan = linspace(0,10,25);
y0 = [0];
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t,y] = ode45(@fun2, tspan, y0, opts);

dlmwrite('Test32_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test32_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 33
tspan = linspace(0,2,25);
y0 = [2 4 6];
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t,y] = ode45(@fun3, tspan, y0, opts);

dlmwrite('Test33_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test33_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 34
tspan = linspace(0,20,25);
y0 = [2 0];
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t,y] = ode45(@fun4, tspan, y0, opts);

dlmwrite('Test34_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test34_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 35
A = 1;
B = 2;
tspan = linspace(0,5,25);
y0 = [0 0.01];
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t,y] = ode45(@(t,y) fun5(t,y,A,B), tspan, y0, opts);

dlmwrite('Test35_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test35_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 36
tspan = linspace(0,3,25);
y0 = [-3 -2 -1 0 1 2 3];
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t,y] = ode45(@fun6, tspan, y0, opts);

dlmwrite('Test36_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test36_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 37
tspan = linspace(0,10,25);
y0 = [2];

opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'Refine',10,'NormControl',true,'MaxStep',1,'InitialStep',0.1);
[t,y] = ode45(@fun1, tspan, y0, opts);

dlmwrite('Test37_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test37_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')


%% Test 38
A = 1;
B = 2;
tspan = linspace(0,5,25);
y0 = [0 0.01];

opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'Refine',1,'NormControl',false,'MaxStep',0.5,'InitialStep',0.1);
[t,y] = ode45(@(t,y) fun5(t,y,A,B), tspan, y0, opts);

dlmwrite('Test38_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test38_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 39
tspan = linspace(0,3,25);
y0 = [0 0];
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
[t,y] = ode45(@fun9,tspan,y0,opts);

dlmwrite('Test39_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test39_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 40
tspan = linspace(0,4,25);
y0 = [0; 4; 2; 20; -pi/2; 2];
opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'Mass',@mass_fun8);
[t, y] = ode45(@fun8,tspan,y0, opts);

dlmwrite('Test40_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test40_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Test 41
tspan = [0,1,2,3,4,5,6,7,8,9,10];
y0 = [2];
opts = odeset('RelTol',1e-12,'AbsTol',1e-15);
[t,y] = ode45(@fun1, tspan, y0, opts);

dlmwrite('Test41_Matlab_t.txt ',t,'delimiter',' ','precision','%.16e')
dlmwrite('Test41_Matlab_y.txt ',y,'delimiter',' ','precision','%.16e')

%% Function test

function [dydt] = fun1(t,y)
    dydt = -0.5 * y;
end

function [dydt] = fun2(t,y)
    dydt = cos(t);
end

function [dydt] = fun3(t,y)
    dydt = [y(2)+y(3), -y(1)+2*y(2)+y(3), y(1)+y(3)]';
end

function [dydt] = fun4(t,y)
    dydt = zeros(2,1);
    dydt(1) = y(2);
    dydt(2) = (1-y(1)^2)*y(2)-y(1);
end

function [dydt] = fun5(t,y,A,B)
    dydt = zeros(2,1);
    dydt(1) = y(2);
    dydt(2) = (A/B)*t.*y(1);
end

function [dydt] = fun6(t,y)
    dydt = -2*y + 2*cos(t)*sin(2*t);
end

 function dydt = fun7(t,y)
    mu = 1 / 82.45;
    mustar = 1 - mu;
    r13 = ((y(1) + mu)^2 + y(2)^2) ^ 1.5;
    r23 = ((y(1) - mustar)^2 + y(2)^2) ^ 1.5;
    dydt = [ y(3)
    y(4)
    2*y(4) + y(1) - mustar*((y(1)+mu)/r13) - mu*((y(1)-mustar)/r23)
    -2*y(3) + y(2) - mustar*(y(2)/r13) - mu*(y(2)/r23) ];
 end

function [value,isterminal,direction] = events_fun7(t,y)
    y0 = [1.2; 0; 0; -1.04935750];
    dDSQdt = 2 * ((y(1:2)-y0(1:2))' * y(3:4));
    value = [dDSQdt; dDSQdt];
    isterminal = [1;  0];
    direction  = [1; -1];
end

function dydt = fun8(t,y)
    m1 = 0.1;
    m2 = 0.1;
	g = 9.81;
    L = 1;
    dydt = [ y(2)
        m2*L*y(6)^2*cos(y(5))
        y(4)
        m2*L*y(6)^2*sin(y(5))-(m1+m2)*g
        y(6)
        -g*L*cos(y(5)) ];
end

function M = mass_fun8(t,y)
    m1 = 0.1;
    m2 = 0.1;
    L = 1;
    M = zeros(6,6);
    M(1,1) = 1;
    M(2,2) = m1 + m2;
    M(2,6) = -m2*L*sin(y(5));
    M(3,3) = 1;
    M(4,4) = m1 + m2;
    M(4,6) = m2*L*cos(y(5));
    M(5,5) = 1;
    M(6,2) = -L*sin(y(5));
    M(6,4) = L*cos(y(5));
    M(6,6) = L^2;
end

function dydt = fun9(t,y)
    dydt = [y(2);-5*y(2)+4*y(1)+sin(10*t)];
end

