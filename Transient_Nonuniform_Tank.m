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
L = 60;             %m
m_dot_HTF = 547;    %kg/s
N = 2e5;            %number of tubes
U = 5;              %W/m^2K
P = 2*pi*r_i*N;     %m


%%%Calculated values%%%
Ac_HTF = pi*(r_o^2 - r_i^2)*N;          %m^2
Ac_stor = pi*r_i^2*N;                   %m^2
u_m = m_dot_HTF/(rho_HTF*Ac_HTF);       %m/s
m_prime_stor = rho_stor*Ac_stor;        %kg/m

%%%Set time and grid steps%%%
time_nodes = 1000;  
t_final = 13*3600;                      %seconds
dt = t_final/time_nodes;                %seconds
time = 0:dt/3600:t_final/3600;          %hours

length_nodes = 1000;
dx = L/length_nodes;
length = 0:dx:L;
x_final = max(size(length));

%%%Set inlet conditions%%%
theta_stor_initial = 1;         %dimensionless; theta = (T-T_low)/(T_high-T_low)
theta_HTF_return = 0;           %dimensionless; theta = (T-T_low)/(T_high-T_low)

T_stor_initial = 500;           %deg C (aka T_high)
T_HTF_return = 289;             %deg C (aka T_low)

theta_stor(1,1) = theta_stor_initial;
theta_HTF(1,1) = theta_HTF_return;

%%%Peng-Robinson Lookup Table for Naphthalene%%%
[PR_T] = xlsread('Naphthalene_PengRobinson_LookupTable', 'u lookup table', 'B2:B3001');     %Varying T
[PR_rho] = xlsread('Naphthalene_PengRobinson_LookupTable', 'u lookup table', 'C1:W1');      %Varying rho
[PR_u] = xlsread('Naphthalene_PengRobinson_LookupTable', 'u lookup table', 'C2:W3001');     %Internal energy varying with T and rho
[PR_P] = xlsread('Naphthalene_PengRobinson_LookupTable', 'P lookup table', 'C2:W3001');     %Pressure varying with T and rho

%%%%%%%%%%%%%%%MAIN ROUTINE%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:time_nodes
   
    theta_HTF(1,t) = theta_HTF_return;
    
    for x = 1:length_nodes
        Q_stor(x,t) = U*P*(theta_stor(x,t) - theta_HTF(x,t));
        theta_HTF(x+1,t) = theta_HTF(x,t) + dx/(m_dot_HTF*cp_HTF)*Q_stor(x,t);
        
        if x == 1
            theta_stor(x+1,1) = theta_stor_initial;                                                 %This If statement has identical equations.  Its purpose is to fill the first row in x because this is a backwards scheme.
            theta_stor(x,t+1) = theta_stor(x+1,t) - dt/(m_prime_stor*cp_stor)*Q_stor(x,t);
            theta_stor(x+1,t+1) = theta_stor(x,t+1);
        else
            theta_stor(x+1,1) = theta_stor_initial;
            theta_stor(x+1,t+1) = theta_stor(x+1,t) - dt/(m_prime_stor*cp_stor)*Q_stor(x,t);    
        end
       
    end

    if theta_HTF(x_final,t) < 0.4787                                                                %Bypass loop
        theta_HTF_return = 0.4351*theta_HTF(x_final,t)-0.2083;                                      %Kolb function with non-dimensional temperature
    end
    
end

%%%Clean-up
theta_HTF(:,time_nodes+1) = theta_HTF(:,time_nodes);                    %This is to copy an extra step in time at the end to make Time and Theta_HTF the same length.
theta_HTF(x_final,time_nodes+1) = theta_HTF(x_final-1,time_nodes+1);    %This is to populate one empty cell at the last time and last grid step

%%%Dimensionalization%%%
T_stor = theta_stor*(T_stor_initial - T_HTF_return) + T_HTF_return;
T_HTF = theta_HTF*(T_stor_initial - T_HTF_return) + T_HTF_return;

%%%Power block%%%
Q_dot_turbine = 0.344*T_HTF(x_final,:)-84.16;                           %MW, Kolb function

%%%Plotting variables%%%
for t = 1:time_nodes
    if time(t) >= 1     %For plotting, T1 = 1 hr. mark (can be changed)
        T1 = t;
    break
    end
end

for t = 1:time_nodes
    if time(t) >= 3     %%For plotting, T2 = 3 hr. mark (can be changed)
        T2 = t;
    break
    end
end

for t = 1:time_nodes
    if time(t) >= 6     %For plotting, T3 = 6 hr. mark (can be changed)
        T3 = t;
    break
    end
end

for t = 1:time_nodes
    if time(t) >= 12     %For plotting, T4 = 12 hr. mark (can be changed)
        T4 = t;
    break
    end
end

x_mid = round(x_final/2);

%%%Plots%%%
figure (1)
hold on
plot(length, T_stor(:,T1),'-r',length,T_stor(:,T2),'--r',length,T_stor(:,T3),':r',length,T_stor(:,T4),'xr')
plot(length, T_HTF(:,T1),'-b',length,T_HTF(:,T2),'--b',length,T_HTF(:,T3),':b',length,T_HTF(:,T4),'xb')
plot(length, 390, '-k')
legend('Stor (t = 1 hr)', 'Stor (t = 3 hr)', 'Stor (t = 6 hr)', 'HTF (t = 1 hr)', 'HTF (t = 3 hr)', 'HTF (t = 6 hr)', 'Drop-off temp.')
xlabel('Length (m)')
ylabel('Temperature (^oC)')
axis([0 L 0 600])

figure (2)
hold on
plot(time, T_stor(x_mid,:),'-r',time, T_stor(x_final, :),'--r')         %Plot storage temp. for all times @ halfway and exit
plot(time, T_HTF(x_mid,:),'-b',time, T_HTF(x_final,:),'--b')            %Plot HTF temp. for all times @ halfway and exit
plot(time, 390, '-k')
legend('Stor (halfway)', 'Stor (exit)', 'HTF (halfway)', 'HTF (exit)', 'Drop-off temp.')
xlabel('Time (hours)')
ylabel('Temperature (^oC)')
axis([0 t_final/3600 0 600])

figure (3)
plot(time, Q_dot_turbine,'-r')
xlabel('Time (hours)')
ylabel('Energy output (MW)')
axis([0 12 0 100])

figure (4)
plot(PR_rho, PR_u(1,:),PR_rho, PR_u(500,:),PR_rho, PR_u(2000,:),PR_rho, PR_u(3000,:))
