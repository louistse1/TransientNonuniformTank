%%%%%System Model of a Supercritical Thermal Energy Storage System
%%%%%Written by Dr. Adrienne Lavine and Louis Tse
%%%%%Dept. of Mechanical & Aerospace Engineering

clear all; clc; close all;

%%%Set system specs%%%
alpha_HTF = 1e-7;   %m^2/s
rho_HTF = 1000;     %kg/m^3
cp_HTF = 3000;      %J/kgK

alpha_stor = 1e-7;  %m^2/s
rho_stor = 1000;    %kg/m^3
cp_stor = 3000;     %J/kgK

r_i = 2.5e-2;       %m
r_o = 2.6e-2;       %m
L = 100;            %m
m_dot_HTF = 500;    %kg/s
N = 2e5;            %number of tubes
U = 5;              %W/m^2K
P = 2*pi*r_i;       %m

%%%Calculated values%%%
Ac_HTF = pi*(r_o^2 - r_i^2)*N;          %m^2
Ac_stor = pi*r_i^2*N;                   %m^2
u_m = m_dot_HTF/(rho_HTF*Ac_HTF);       %m/s
m_prime_stor = rho_stor*Ac_stor;        %kg/m

%%%Set time and grid steps%%%
time_nodes = 1000;  
dt = 9e-3*(r_i^2/alpha_stor/3600);      %hours
t_final = dt*time_nodes;
time = 0:dt:t_final;

length_nodes = 500;
x_final = L;
dx = x_final/length_nodes;
length = 0:dx:x_final;

%%%Set inlet conditions%%%
theta_stor(1,1) = 1;
theta_HTF(1,1) = 0;


%%%%%%%%%%%%%%%MAIN ROUTINE%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:time_nodes
   
    theta_HTF(1,t) = 0;
    
    for x = 1:length_nodes
        Q_stor(x,t) = U*P*(theta_stor(x,t) - theta_HTF(x,t));
        theta_HTF(x+1,t) = theta_HTF(x,t) + N*dx/(m_dot_HTF*cp_HTF)*Q_stor(x,t);
        
        if x == 1
            theta_stor(x+1,1) = 1;
            theta_stor(x,t+1) = theta_stor(x+1,t) - N*3600*dt/(m_prime_stor*cp_stor)*Q_stor(x,t);
            theta_stor(x+1,t+1) = theta_stor(x,t+1);
        else
            theta_stor(x+1,1) = 1;
            theta_stor(x+1,t+1) = theta_stor(x+1,t) - N*3600*dt/(m_prime_stor*cp_stor)*Q_stor(x,t);
            
        end
       
    end

end

for x = 1:length_nodes
    theta_HTF(x,time_nodes+1) = theta_HTF(x,time_nodes);    %This For loop is to copy an extra step in time at the end to make Time and Theta_HTF the same length.
end
    
%theta=(T-T_min)/(T_max-T_min)

%%Plots%%%
% figure
hold on
plot(length, theta_stor(:,65),'-k')
plot(length, theta_stor(:, 193),'--k')
% plot(length, theta_HTF(:,65),'-b')
% plot(length, theta_HTF(:,193),'--b')


legend('Stor (t = 1 hr)', 'Stor (t = 3 hr)', 'HTF (t = 1 hr)', 'HTF (t = 3 hr)')
xlabel('Length (m)')
ylabel('Nondimensional temperature (T/T_s_t_o_r_,_i_n_i_t_i_a_l)')

figure
hold on
plot(time, theta_stor(length_nodes/2,:),'-r')
plot(time, theta_stor(end, :),'--r')
plot(time, theta_HTF(length_nodes/2,:),'-b')
plot(time, theta_HTF(end,:),'--b')

legend('Stor (halfway)', 'Stor (exit)', 'HTF (halfway)', 'HTF (exit)')
xlabel('Time (hours)')
ylabel('Nondimensional temperature (T/T_s_t_o_r_,_i_n_i_t_i_a_l)')

