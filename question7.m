%Question 7
%Parameter values
c = 0.05; %Fluid flow coefficient, c* (s^-1)
p_infine = 20e6; %Wall rock pressure, p∞ (Pa)
beta = 1e-10; %Compressibility, β (Pa^−1)
v = 1e-6; %Slip rate, v (m s^−1)
d = 10e-6; %Characteristic slip, dc (m)
porosity_ss = 0.03; %Characteristic porosity, φss
porosity_i = 0.025; %Initial porosity, φ0

N=12;
h_value=logspace(-2,1,N);%Generate logarithmically spaced step sizes

% Initialize arrays to store errors for both methods
euler_pressure_error=zeros(1,N);% Initialize array for Euler method pressure error
euler_porosity_error=zeros(1,N);% Initialize array for Euler method porosity error
modify_pressure_error=zeros(1,N);% Initialize array for modified midpoint method pressure error
modify_porosity_error=zeros(1,N);% Initialize array for modified midpoint method porosity error

for k=1:N
% discretise
dh=h_value(k);% step size
t=0:dh:100;% Time vector with current stepsize
S=length(t);

% initialise solution vector with euler method
y_euler=zeros(2,S);% Initialize solution vector for Euler method
y_euler(1,1)=20e6;% Initial pressure
y_euler(2,1)=0.025;% Initial porosity

% initialise solution vector with modify method
y_modify=zeros(2,S);% Initialize solution vector for modified midpoint method
y_modify(1,1)=20e6;% Initial pressure
y_modify(2,1)=0.025;% Initial porosity

% Main loop
for i=1:S-1
y_euler(:,i+1) = y_euler(:,i) + dh*f(t(i),y_euler(:,i),c ,p_infine, beta, v, d, porosity_ss, porosity_i);
y_modify(:,i+1) = y_modify(:,i) + dh*f(t(i) + dh/2, y_modify(:,i) + (dh/2)* f(t(i), y_modify(:,i),c ,p_infine, beta, v, d, porosity_ss, porosity_i),c ,p_infine, beta, v, d, porosity_ss, porosity_i);
end

%Exact Solutions
%Exact solution of Pressure in Pa
y1_closed = (((porosity_i - porosity_ss)/beta)/(1- c*d/v)) *(exp(-c*t) - exp(-v*t/d)) + p_infine ;
%Exact solution of Porosity
y2_closed = (porosity_i - porosity_ss)*exp(-v*t/d) + porosity_ss;

% Calculate errors for both methods
euler_pressure_error(k)=mean(abs(y1_closed-y_euler(1,:)));%Pressure error for Euler method
euler_porosity_error(k)=mean(abs(y2_closed-y_euler(2,:)));%Porosity error for Euler method
modify_pressure_error(k)=mean(abs(y1_closed-y_modify(1,:)));%Pressure error for modified midpoint method
modify_porosity_error(k)=mean(abs(y2_closed-y_modify(2,:)));%Porosity error for modified midpoint method
end

% Plotting
loglog(h_value,euler_pressure_error, 'kv', 'linewidth',2)% Plot Euler method pressure error
hold on
loglog(h_value,modify_pressure_error, 'bo', 'linewidth',2)% Plot modified midpoint method pressure error
xlabel('step size')
ylabel('average error')
title('Average Error for Pressure')
legend('euler error', 'mid-point error')


figure()
loglog(h_value,euler_porosity_error, 'y^', 'linewidth',2)% Plot Euler method porosity error
hold on
loglog(h_value,modify_porosity_error, 'go', 'linewidth',2)% Plot modified midpoint method porosity error
xlabel('step size')
ylabel('average error')
title('Average Error for Porosity')
legend('euler error', 'mid-point error')