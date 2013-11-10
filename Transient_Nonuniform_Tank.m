%%%%%System Model of a Supercritical Thermal Energy Storage System
%%%%%Written by Dr. Adrienne Lavine and Louis Tse
%%%%%Dept. of Mechanical & Aerospace Engineering

clear all; clc; close all;

%%%%%%%%%%%PROCEDURE%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(1) Guess a value for the m_stor, the mass of storage fluid.
%(2) Main routine will calculate T_stor(x,t) and T_HTF(x,t).
%(3) T_stor(x,t) at beginning and final time are pinpointed on Peng-Robinson
     %look-up table to determine u_initial and u_final.
%(4) E_stor = m_storage*delu, where E_stor and delu are known.  m_storage is calculated, and
     %iterated upon to match m_stor and m_storage.

%%%%%%%END PROCEDURE%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
 
rho_stor = 350;                  %[kg/m^3]  Storage fluid loading
z=1;

while rho_stor < 650

%%%Guess a value of m_stor   
m_stor = 3.55e7;                  %[kg] Total mass of storage fluid
m_stor_percent_diff = 1;          %Assume 100% percent difference between m_stor and m_storage to trigger iterative loop

while m_stor_percent_diff > 0.01  %Criterion for iterating...if m_stor and m_storage are more than 3% apart
    
%%%System specifications
E_stor = 1621*1000;              %[kWh]  Energy storage capacity
discharge_time = 12*3600;        %[s]  Total discharge time
m_dot_e = 547;                   %[kg/s] Total mass flow rate to evaporator/power block

%%%Heat transfer fluid properties [1]
rho_HTF = 1060;     %[kg/m^3]  Density of HTF
cp_HTF = 2275;      %[J/kgK]  Specific heat of HTF
k_HTF = 0.098;      %[W/mK] Thermal conductivity of HTF
alpha_HTF = k_HTF/(rho_HTF*cp_HTF);   %[m^2/s]  Thermal diffusivity of HTF
price_HTF = 3.96;           %[$/kg]   Price of HTF

%%%Storage fluid properties  [2]
cp_stor = 2446;     %[J/kgK]  Specific heat of storage fluid
k_stor = 0.13;      %[W/mK]  Thermal conductivity of storage fluid
alpha_stor = k_stor/(rho_stor*cp_stor);  %[m^2/s]  Thermal diffusivity of storage fluid
price_stor = 1.50;  %[$/kg]  Price of storage fluid

%%%Tube dimensions
r_i = 0.025;        %[m]  Radius of a single tube, with zero wall thickness  
r_o = 0.026;        %[m]  Annulus of HTF surrounding a single tube, which implies Gap = (r_o - r_i)
L = 60;            %[m]  Length of a single tube            
U = 2;              %[W/m^2K]  Overall heat transfer coefficient
N = m_stor/(pi*r_i^2*rho_stor*L);    %Number of tubes
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
m_prime_stor = rho_stor*Ac_stor;        %[kg/m]    Mass of storage fluid per length of all tubes
V_stor = pi*(r_i)^2*L*N;                %[kg/m^3]  Volume of storage fluid
Ac_HTF = pi*(r_o^2 - r_i^2)*N;          %[m^2]
u_m = m_dot_e/(rho_HTF*Ac_HTF);         %[m/s]

%%%Set time and grid steps%%%
time_nodes = 500;              
dt = discharge_time/time_nodes;                %[s]
time = 0:dt/3600:discharge_time/3600;          %[h]
t_final = max(size(time));      %index of final time node

length_nodes = 1000;
dx = L/length_nodes;
tank_length = 0:dx:L;
x_final = max(size(tank_length));       %index of final length node

%%%Set inlet conditions%%%
m_dot_s(1) = 256.1;                 %[kg/s]  (total guess, but is needed to start the routine)  
m_dot_b(1) = m_dot_e - m_dot_s(1);  %[kg/s]  Mass flow rate through the bypass (to cool down the hot HTF exiting the tank)

T_stor_initial = 449;      %[deg C] (aka T_high)
T_eo(1) = 289;             %[deg C] (aka T_low)
T_ei(1) = 390;             %[deg C] (design value of turbine)

T_stor(1,1) = T_stor_initial;
T_HTF(1,1) = T_eo(1);
  
%%%%%%%%%%%%%%%MAIN ROUTINE%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:time_nodes
    
    if t > 1
    T_HTF(1,t) = T_eo(t-1);
    end
    
    for x = 1:length_nodes
            
        Q_stor(x,t) = U*P*(T_stor(x,t) - T_HTF(x,t));
            
            T_HTF(x+1,t) = T_HTF(x,t) + dx/(m_dot_s(t)*cp_HTF)*Q_stor(x,t);

            T_stor(x+1,1) = T_stor_initial;
            T_stor(x+1,t+1) = T_stor(x+1,t) - dt/(m_prime_stor*cp_stor)*Q_stor(x,t); 
            
            if x == 1
            T_stor(1,t+1) = T_stor(2,t+1);     %This copies the 2nd row in x to the 1st row (which was empty)
            end
             
            T_stor_avg(t) = mean(T_stor(:,t));
    end
  
    %Bypass loop is on (by default).  This means T_ei is set, and m_dot_b
    %is solved for.
    T_ei(t) = 390; 
    T_eo(t) = 0.433*T_ei(t) + 120.13;    %Kolb function [1]
    
    m_dot_b(t+1) = m_dot_e*(T_ei(t)-T_HTF(x_final,t))/(T_eo(t)-T_HTF(x_final,t));   %Conservation of enthalpy
    m_dot_s(t+1) = m_dot_e - m_dot_b(t);    %Conservation of mass
    
    %Bypass loop is off.  This means m_dot_b is set to zero, and T_ei is
    %equal to the HTF temp. exiting the tank.
    if m_dot_b(t+1) <= 0
        T_ei(t) = T_HTF(x_final,t);
        T_eo(t) = 0.433*T_ei(t) + 120.13;
        
        m_dot_b(t+1) = 0;
        m_dot_s(t+1) = m_dot_e - m_dot_b(t);
        
    end
end

%%%Power block
Q_dot_turbine = 0.344*T_ei-84.16;      %[MW], Kolb function [1]
Q_total = sum(Q_dot_turbine*dt/3600)   %[MWh]

%%%Clean-up

%These are to copy an extra step in time for each vector to make it the same length as the Time vector.
T_HTF(:,time_nodes+1) = T_HTF(:,time_nodes);                    
T_stor_avg(time_nodes+1) = T_stor_avg(time_nodes);
T_ei(t_final) = T_ei(t_final-1);
m_dot_b(time_nodes+1) = m_dot_b(time_nodes);              
Q_dot_turbine(t_final) = Q_dot_turbine(t_final-1);
%These are to copy an extra step in x to make the vector the same length as the Length vector.
T_HTF(x_final,time_nodes+1) = T_HTF(x_final-1,time_nodes+1);  

%%%%%%%%%%%%PENG-ROBINSON LOOKUP TABLE%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Naphthalene Lookup Table (generated using EES [2].  If using a different fluid, an entirely new lookup table needs to be generated.)
[PR_T] = xlsread('Naphthalene_PengRobinson_LookupTable', 'u lookup table', 'B2:B3001');     %Varying T [K]
[PR_rho] = xlsread('Naphthalene_PengRobinson_LookupTable', 'u lookup table', 'C1:W1');      %Varying rho [kg/m^3]
[PR_u] = xlsread('Naphthalene_PengRobinson_LookupTable', 'u lookup table', 'C2:W3001');     %Internal energy varying with T and rho [kJ/kg]
[PR_P] = xlsread('Naphthalene_PengRobinson_LookupTable', 'P lookup table', 'C2:W3001');     %Pressure varying with T and rho [kPa]

%%%Find the user-defined rho in the lookup table
diff_rho = PR_rho - rho_stor;  %[kg/m^3]            %Find the difference between all the rho values in lookup table and user-defined rho
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
delu = u_initial - u_final;                         %[kJ/kg]
m_storage = E_stor*3600/delu;                       %[kg]

%%%Calculate max pressure to find required wall thickness
P_max = PR_P(index_Tinitial, index_rho)/1000;       %[MPa]   %Maximum pressure occurs at initial (highest T_stor)

%%%Calculate the percent diff. between the guessed and calculated m_stor, and iterate using bisection method
m_stor_percent_diff = abs(m_storage - m_stor)/m_storage      
m_stor = (m_storage + m_stor)/2;

end

if z==1
    results(z,:) = [rho_stor, m_stor, m_stor_percent_diff, delu, Q_total, T_stor_avg(t_final), T_stor(x_final,t_final), T_HTF(x_final,t_final), P_max];
else
    results(z,:) = [rho_stor, m_stor, m_stor_percent_diff, delu, Q_total, T_stor_avg(t_final), T_stor(x_final,t_final), T_HTF(x_final,t_final), P_max];
    results = [results; results(z,:)];
end

rho_stor = rho_stor + 100;
z=z+1;
end

%%%%%%%%%%%%%%%%%%%%%COST MODELING%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Calculating tube parameters
t = P_max*(r_i)/(sigma*(1+2*P_max/sigma));  %[m] Required wall thickness
r_o_tube = r_i + t;      %[m]  Outer radius of tube
V_SS = pi*(r_o^2 - r_i^2)*L*N;   %[m^3]  Total volume of tube material
m_SS = rho_SS*V_SS;              %[kg]  Total mass of tube material
E_SS = m_SS*cp_SS*(max(T_stor_avg) - min(T_stor_avg))/(1e3*3600);    %[kWh]

%%%Calculating costs
cost_SS = m_SS*price_SS/(E_stor+E_SS);   %[$/kWh]   Cost of tube material  
cost_stor = m_stor*price_stor/(E_stor+E_SS); %[$/kWh]  Cost of storage fluid
cost_total = cost_SS + cost_stor;  %[$/kWh]  Total cost



%%%%%%%%%%%%%%%%%%%%%%%%PLOTS%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Plotting variables
T1 = find(time>1,1);        %For plotting, T1 = 1 hr. mark (can be changed)
T2 = find(time>3,1);        %For plotting, T2 = 5 hr. mark (can be changed)
T3 = find(time>6,1);        %For plotting, T3 = 10 hr. mark (can be changed)
T4 = find(time>9,1);        %For plotting, T3 = 10 hr. mark (can be changed)
x_mid = round(x_final/2);   %Index of the middle of tank

%%%Plots%%%
figure (1)
hold on
subplot(2,2,1), plot(tank_length, T_stor(:,T1),'-r',tank_length,T_HTF(:,T1),'-b',tank_length, 390, '--k')
title('t = 1 hr')
xlabel('Length (m)')
ylabel('Temperature (^oC)')
legend('Stor', 'HTF')
subplot(2,2,2), plot(tank_length, T_stor(:,T2),'-r',tank_length,T_HTF(:,T2),'-b',tank_length, 390, '--k')
title('t = 3 hr')
xlabel('Length (m)')
ylabel('Temperature (^oC)')
subplot(2,2,3), plot(tank_length, T_stor(:,T3),'-r',tank_length,T_HTF(:,T3),'-b',tank_length, 390, '--k')
title('t = 6 hr')
xlabel('Length (m)')
ylabel('Temperature (^oC)')
subplot(2,2,4), plot(tank_length, T_stor(:,T4),'-r',tank_length,T_HTF(:,T4),'-b',tank_length, 390, '--k')
title('t = 9 hr')
xlabel('Length (m)')
ylabel('Temperature (^oC)')


figure (2)
hold on
subplot(2,1,1), plot(time, T_stor(x_mid,:),'--r',time, T_stor(x_final,:),'-r', time, 390, ':k')             
xlabel('Time (hours)')
ylabel('Temperature (^oC)')
legend('Halfway', 'Exit')
title('Storage fluid transient temperature')
axis([0 discharge_time/3600 200 600])
subplot(2,1,2), plot(time,T_HTF(x_mid,:),'--b', time, T_HTF(x_final,:),'-b', time, 390, ':k')  
xlabel('Time (hours)')
ylabel('Temperature (^oC)')
legend('Halfway', 'Exit')
title('HTF transient temperature')
axis([0 discharge_time/3600 200 600])


figure (3)
plot(time, Q_dot_turbine,'-b')
xlabel('Time (hours)')
ylabel('Energy output (MW)')
title('Energy output during discharge cycle')
axis([0 discharge_time/3600 0 60])


figure (4)
plot(time,T_stor(x_final,:),'-r', time,T_HTF(x_final,:),'-b', time,T_stor_avg,'--r', time,m_dot_b,':m')
xlabel('Time (hours)')
ylabel('Temperature (^oC)')
legend('T_s_t_o_r (exit)', 'T_H_T_F (exit)', 'T_s_t_o_r_,_a_v_g', 'm_b')
title('Temperatures during discharge cycle')
axis([0 discharge_time/3600 0 600])

%%%%%%%%%REFERENCES%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[1]  Kolb, G.J.  “Evaluation of annual performance of 2-tank and thermocline thermal storage systems for trough plants,” Journal of Solar Energy Engineering, Vol. 133, August 2011.
%[2]  Peng, Ding-Yu, Robinson, Donald B., “A new two-constant equation of state.”  Ind. Eng. Chem. Fundamen., Vol. 15, 1, pp. 59-64, 1976.