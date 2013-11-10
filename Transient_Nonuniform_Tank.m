%%%%%System Model of a Supercritical Thermal Energy Storage System
%%%%%Written by Dr. Adrienne Lavine and Louis Tse
%%%%%Dept. of Mechanical & Aerospace Engineering

clear all; clc; close all;

%%%Set system specs%%%

E_stor = 1621*1000; %kWh  Energy storage capacity of system
alpha_HTF = 1e-7;   %m^2/s
rho_HTF = 1000;     %kg/m^3
cp_HTF = 3000;      %J/kgK

alpha_stor = 1e-7;  %m^2/s
rho_stor = 600;    %kg/m^3
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
time_nodes = 100;  
discharge_time = 12*3600;                      %seconds
dt = discharge_time/time_nodes;                %seconds
time = 0:dt/3600:discharge_time/3600;          %hours
t_final = max(size(time));

length_nodes = 1000;
dx = L/length_nodes;
length = 0:dx:L;
x_final = max(size(length));

%%%Set inlet conditions%%%
theta_stor_initial = 1;         %dimensionless; theta = (T-T_low)/(T_high-T_low)
theta_eo = 0;           %dimensionless; theta = (T-T_low)/(T_high-T_low)

T_stor_initial = 500;           %deg C (aka T_high)
T_eo = 289;             %deg C (aka T_low)

theta_stor(1,1) = theta_stor_initial;
theta_HTF(1,1) = theta_eo;

T_stor(1,1) = T_stor_initial;
T_HTF(1,1) = T_eo;

%%%Peng-Robinson Lookup Table for Naphthalene%%%
[PR_T] = xlsread('Naphthalene_PengRobinson_LookupTable', 'u lookup table', 'B2:B3001');     %Varying T
[PR_rho] = xlsread('Naphthalene_PengRobinson_LookupTable', 'u lookup table', 'C1:W1');      %Varying rho
[PR_u] = xlsread('Naphthalene_PengRobinson_LookupTable', 'u lookup table', 'C2:W3001');     %Internal energy varying with T and rho
[PR_P] = xlsread('Naphthalene_PengRobinson_LookupTable', 'P lookup table', 'C2:W3001');     %Pressure varying with T and rho

%%%%%%%%%%%%%%%MAIN ROUTINE%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:time_nodes
   
    theta_HTF(1,t) = theta_eo;
    T_HTF(1,t) = T_eo(t);
    
    for x = 1:length_nodes
        Q_stor(x,t) = U*P*(theta_stor(x,t) - theta_HTF(x,t));
        
        theta_HTF(x+1,t) = theta_HTF(x,t) + dx/(m_dot_HTF*cp_HTF)*Q_stor(x,t);
        T_HTF(x+1,t) = theta_HTF(x+1,t)*(T_stor_initial - T_eo(t)) + T_eo(t);
        
        if x == 1
            theta_stor(x+1,1) = theta_stor_initial;                                                 %This If statement has identical equations.  Its purpose is to fill the first row in x because this is a backwards scheme.
            theta_stor(x,t+1) = theta_stor(x+1,t) - dt/(m_prime_stor*cp_stor)*Q_stor(x,t);
            theta_stor(x+1,t+1) = theta_stor(x,t+1);
            
            T_stor(x+1,1) = T_stor_initial;
            T_stor(x,t+1) = theta_stor(x,t+1)*(T_stor_initial - T_eo(t)) + T_eo(t);
            T_stor(x+1,t+1) = theta_stor(x+1,t+1)*(T_stor_initial - T_eo(t)) + T_eo(t);
        else
            theta_stor(x+1,1) = theta_stor_initial;
            theta_stor(x+1,t+1) = theta_stor(x+1,t) - dt/(m_prime_stor*cp_stor)*Q_stor(x,t);
            
            T_stor(x+1,1) = T_stor_initial;
            T_stor(x+1,t+1) = theta_stor(x+1,t+1)*(T_stor_initial - T_eo(t)) + T_eo(t);
        end
       
         T_stor_avg(t) = mean(T_stor(:,t));   %Will be used to calculate E_stor = m_storage*delu, where delu is calculated from mean value of T_stor @ beginning and final time
    end

    %Bypass loop is on (by default).  This means T_ei is set, and m_dot_b
    %is solved for.
    T_ei(t) = 390;
    T_ei(t+1) = 390; 
    T_eo(t) = 0.433*T_ei(t) + 120.13;
    T_eo(t+1) = 0.433*T_ei(t+1) + 120.13;
    
    m_dot_b(t) = m_dot_HTF*(T_ei(t)-T_HTF(x_final,t))/(T_eo(t)-T_HTF(x_final,t)); 
    
    %Bypass loop is off.  This means m_dot_b is set to zero, and T_ei is
    %equal to the HTF temp. exiting the tank.
    if m_dot_b(t) <= 0
        m_dot_b(t) = 0;
        T_ei(t) = T_HTF(x_final,t);
        T_eo(t) = 0.433*T_ei(t) + 120.13;
    end
    
end

%%%Power block%%%
Q_dot_turbine = 0.344*T_ei-84.16;                           %MW, Kolb function

%%%Clean-up
theta_HTF(:,time_nodes+1) = theta_HTF(:,time_nodes);                    %This is to copy an extra step in time at the end to make Time and Theta_HTF the same length.
theta_HTF(x_final,time_nodes+1) = theta_HTF(x_final-1,time_nodes+1);    %This is to populate one empty cell at the last time and last grid step

T_HTF(:,time_nodes+1) = T_HTF(:,time_nodes);                    
T_HTF(x_final,time_nodes+1) = T_HTF(x_final-1,time_nodes+1);  
T_stor_avg(:,time_nodes+1) = T_stor_avg(:,time_nodes); 

m_dot_b(time_nodes+1) = m_dot_b(time_nodes);                            %This is to copy an extra step in time to m_dot_b to make it the same length as Time vector.
Q_dot_turbine(t_final) = Q_dot_turbine(t_final-1);

%%%%%%%%%%%%PENG-ROBINSON LOOKUP TABLE%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Naphthalene Lookup Table (generated using EES [2].  If using a different fluid, an entirely new lookup table needs to be generated.)
[PR_T] = xlsread('Naphthalene_PengRobinson_LookupTable', 'u lookup table', 'B2:B3001');     %Varying T [K]
[PR_rho] = xlsread('Naphthalene_PengRobinson_LookupTable', 'u lookup table', 'C1:W1');      %Varying rho [kg/m^3]
[PR_u] = xlsread('Naphthalene_PengRobinson_LookupTable', 'u lookup table', 'C2:W3001');     %Internal energy varying with T and rho [kJ/kg]
[PR_P] = xlsread('Naphthalene_PengRobinson_LookupTable', 'P lookup table', 'C2:W3001');     %Pressure varying with T and rho [kPa]

%%%Find the user-defined rho in the lookup table
diff_rho = PR_rho - rho_stor;  %[kg/m^3]                 %Find the difference between all the rho values in lookup table and user-defined rho
mindiff_rho = min(abs(diff_rho));                   %Minimize all differences to find closest value in lookup table to user-defined rho
index_rho = find(diff_rho == mindiff_rho);          %Find the index of the smallest difference, to pin down a value in the lookup table

%%%Find the user-defined T_initial in the lookup table
diff_Tinitial = PR_T - (T_stor_initial+273);          %[deg C]
mindiff_Tinitial = min(abs(diff_Tinitial));
index_Tinitial = find(abs(diff_Tinitial) == mindiff_Tinitial);

%%%Find the simulated T_final in the lookup table
diff_Tfinal = PR_T - (T_stor_avg(t_final)+273);       %[deg C]
mindiff_Tfinal = min(abs(diff_Tfinal));
index_Tfinal = find(abs(diff_Tfinal) == mindiff_Tfinal);

%%%Calculate delu to find required m_storage
u_final = PR_u(index_Tfinal, index_rho);            %[kJ/kg]
u_initial = PR_u(index_Tinitial, index_rho);        %[kJ/kg]
delu = u_initial - u_final                         %[kJ/kg]
m_storage = E_stor*3600/delu                        %[kg]

%%%Calculate max pressure to find required wall thickness
P_max = PR_P(index_Tinitial, index_rho)/1000;       %[MPa]   %Maximum pressure occurs at initial (highest T_stor)




%%%Plotting variables%%%
T1 = find(time>1,1);        %For plotting, T1 = 1 hr. mark (can be changed)
T2 = find(time>3,1);        %For plotting, T2 = 3 hr. mark (can be changed)
T3 = find(time>6,1);        %For plotting, T3 = 6 hr. mark (can be changed)
T4 = find(time>8,1);        %For plotting, T4 = 8 hr. mark (can be changed)
x_mid = round(x_final/2);   %Index of the middle of tank

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
axis([0 discharge_time/3600 0 600])

figure (3)
plot(time, Q_dot_turbine,'-r')
xlabel('Time (hours)')
ylabel('Energy output (MW)')
axis([0 discharge_time/3600 0 60])


