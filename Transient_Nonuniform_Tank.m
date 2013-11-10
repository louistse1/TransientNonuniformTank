%%%%%System Model of a Supercritical Thermal Energy Storage System
%%%%%Written by Dr. Adrienne Lavine and Louis Tse
%%%%%Dept. of Mechanical & Aerospace Engineering

clear all; clc; close all;

%%%System specifications
E_stor = 1621*1000;                 %[kWh]  Energy storage capacity
discharge_length = 12*3600;         %[s]  Total discharge time
m_stor = 2.1e7;                     %[kg] Total mass of storage fluid
m_dot_HTF = 547;    %[kg/s]  Total mass flow rate of HTF

%%%Heat transfer fluid properties [1]
rho_HTF = 1060;     %[kg/m^3]  Density of HTF
cp_HTF = 2275;      %[J/kgK]  Specific heat of HTF
k_HTF = 0.098;       %[W/mK] Thermal conductivity of HTF
alpha_HTF = k_HTF/(rho_HTF*cp_HTF);   %[m^2/s]  Thermal diffusivity of HTF
price_HTF = 3.96;           %[$/kg]   Price of HTF

%%%Storage fluid properties  [2]
rho_stor = 400;     %[kg/m^3]  Storage fluid loading
cp_stor = 2446;     %[J/kgK]  Specific heat of storage fluid
k_stor = 0.13;       %[W/mK]  Thermal conductivity of storage fluid
alpha_stor = k_stor/(rho_stor*cp_stor);  %[m^2/s]  Thermal diffusivity of storage fluid
price_stor = 1.50;  %[$/kg]  Price of storage fluid

%%%Tube dimensions
r_i = 0.025;        %[m]  Radius of a single tube, with zero wall thickness  
r_o = 0.026;        %[m]  Annulus of HTF surrounding a single tube, which implies Gap = (r_o - r_i)
L = 60;             %[m]  Length of a single tube            
U = 3;              %[W/m^2K]  Overall heat transfer coefficient
N = m_stor/(pi*r_i^2*rho_stor*L)    %Number of tubes
P = 2*pi*r_i*N;     %[m]  Total perimeter of all tubes

%%%Tube material properties
F_tu = 505;         %[MPa]  Ultimate tensile strength of SS316L
n = 4;              %Safety factor
d = 0.6;            %Derating factor for SS316L operating up to 500 deg C
sigma = F_tu*d/n;   %[MPa]  Allowable stress
price_SS = 1.40;    %[$/kg]  Price of tube material
cp_SS = 500;        %[J/kgK] Specific heat of tube material
rho_SS = 7990;      %[kg/m^3] Density of tube material

%%%Calculated values%%%
Ac_stor = pi*r_i^2*N;                   %[m^2]    Total cross-sectional area of storage fluid
m_prime_stor = rho_stor*Ac_stor;     %[kg/m]    Mass of storage fluid per length of all tubes
V_stor = pi*(r_i)^2*L*N;                %[kg/m^3]  Volume of storage fluid
Ac_HTF = pi*(r_o^2 - r_i^2)*N;          %m^2
u_m = m_dot_HTF/(rho_HTF*Ac_HTF);       %m/s

%%%Set time and grid steps%%%
time_nodes = 1000;  
dt = discharge_length/time_nodes;                %seconds
time = 0:dt/3600:discharge_length/3600;          %hours
t_final = max(size(time));

length_nodes = 7000;
dx = L/length_nodes;
length = 0:dx:L;
x_final = max(size(length));

%%%Set inlet conditions%%%
T_stor_initial = 500;           %deg C (aka T_high)
T_eo(1) = 289;             %deg C (aka T_low)

T_stor(1,1) = T_stor_initial;
T_HTF(1,1) = T_eo;

%%%%%%%%%%%%%%%MAIN ROUTINE%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:time_nodes
    
    T_HTF(1,t) = T_eo(t);    
    for x = 1:length_nodes
            
        Q_stor(x,t) = U*P*(T_stor(x,t) - T_HTF(x,t))*(T_stor_initial - T_eo(t));
        T_HTF(x+1,t) = T_HTF(x,t) + dx/(m_dot_HTF*cp_HTF)*Q_stor(x,t)/(T_stor_initial - T_eo(t));
        
        if x == 1

            T_stor(x+1,1) = T_stor_initial;
            T_stor(x,t+1) = T_stor(x+1,t) - dt/(m_prime_stor*cp_stor)*Q_stor(x,t)/(T_stor_initial - T_eo(t));
            T_stor(x+1,t+1) = T_stor(x,t+1);
        else

            T_stor(x+1,1) = T_stor_initial;
            T_stor(x+1,t+1) = T_stor(x+1,t) - dt/(m_prime_stor*cp_stor)*Q_stor(x,t)/(T_stor_initial - T_eo(t)); 
        end
       
    end
  
    %Bypass loop is on (by default).  This means T_ei is set, and m_dot_b
    %is solved for.
    T_ei(t) = 390;
    T_ei(t+1) = 390; 
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
Q_dot_turbine = 0.344*T_ei - 84.16;    %MW, Kolb function
Q_total = sum(Q_dot_turbine*dt/3600)   %[MWh]

%%%Clean-up
T_HTF(:,time_nodes+1) = T_HTF(:,time_nodes);                    %This is to copy an extra step in time at the end to make Time and Theta_HTF the same length.
T_HTF(x_final,time_nodes+1) = T_HTF(x_final-1,time_nodes+1);    %T
Q_dot_turbine(t_final) = Q_dot_turbine(t_final-1);


%%%Plotting variables%%%
T1 = find(time>1,1);        %For plotting, T1 = 1 hr. mark (can be changed)
T2 = find(time>3,1);        %For plotting, T2 = 3 hr. mark (can be changed)
T3 = find(time>4,1);        %For plotting, T3 = 6 hr. mark (can be changed)
T4 = find(time>5,1);        %For plotting, T4 = 8 hr. mark (can be changed)
x_mid = round(x_final/2);   %Index of the middle of tank

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

% figure (2)
% hold on
% plot(time, T_stor(x_mid,:),'-r',time, T_stor(x_final, :),'--r')         %Plot storage temp. for all times @ halfway and exit
% plot(time, T_HTF(x_mid,:),'-b',time, T_HTF(x_final,:),'--b')            %Plot HTF temp. for all times @ halfway and exit
% plot(time, 390, '-k')
% legend('Stor (halfway)', 'Stor (exit)', 'HTF (halfway)', 'HTF (exit)', 'Drop-off temp.')
% xlabel('Time (hours)')
% ylabel('Temperature (^oC)')
% axis([0 t_final/3600 0 600])
% 
figure (3)
plot(time, Q_dot_turbine,'-r')
xlabel('Time (hours)')
ylabel('Energy output (MW)')
axis([0 discharge_length/3600 0 60])
% 
% figure (4)
% plot(PR_rho, PR_u(1,:),PR_rho, PR_u(500,:),PR_rho, PR_u(2000,:),PR_rho, PR_u(3000,:))
