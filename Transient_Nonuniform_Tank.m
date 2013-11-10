%%%%%System Model of a Supercritical Thermal Energy Storage System
%%%%%Written by Dr. Adrienne Lavine and Louis Tse
%%%%%Dept. of Mechanical & Aerospace Engineering

clear all; clc; close all;

%%%Set system specs%%%
alpha_HTF = 1e-7;
alpha_stor = 1e-7;
rho_cp_HTF = 3e6;
rho_cp_stor = 3e6;
rho_HTF = 1e3;
r_i = 2.5e-2;
r_o = 2.6e-2;
L = 100;
m_dot_HTF = 500;
N = 2e5;
U = 5;

%%%Calculated values%%%
Ac_HTF = pi*(r_o^2 - r_i^2)*N;
u_m = m_dot_HTF/(rho_HTF*Ac_HTF);
Pe = u_m*r_i/alpha_stor;
U_hat = U*2*r_i/(alpha_stor*rho_cp_HTF);
alpha_ratio = alpha_HTF/alpha_stor;
R = r_i^2/(r_o^2-r_i^2)*(rho_cp_stor/rho_cp_HTF);

%%%Set time and grid steps%%%
time_nodes = 300;  
dtau = 9e-3;
tau_final = dtau*time_nodes;
tau = 0:dtau:tau_final;
time = tau*(2*r_i^2/alpha_stor/3600);

length_nodes = 100;
X_final = L/(r_i*Pe);
dX = X_final/length_nodes;
X = 0:dX:X_final;
length = X*r_i^2*u_m/alpha_stor;

%%%Set inlet conditions%%%
theta_stor(1,1) = 1;
theta_HTF(1,1) = 0;


%%%%%%%%%%%%%%%MAIN ROUTINE%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:time_nodes
   
    theta_HTF(1,t) = 0;
    
    for x = 1:length_nodes
        
        theta_HTF(x+1,t) = theta_HTF(x,t) + U_hat*R*dX*(theta_stor(x,t)-theta_HTF(x,t));
        
        if x == 1
            theta_stor(x+1,1) = 1;
            theta_stor(x,t+1) = theta_stor(x+1,t) + dtau*(U_hat*(theta_HTF(x+1,t)-theta_stor(x+1,t)) + 1/(Pe*dX)^2*(theta_stor(x,t)-2*theta_stor(x+1,t)+theta_stor(x+1,t)));
            theta_stor(x+1,t+1) = theta_stor(x,t+1);
        else
            theta_stor(x+1,1) = 1;
            theta_stor(x+1,t+1) = theta_stor(x+1,t) + dtau*(U_hat*(theta_HTF(x+1,t)-theta_stor(x+1,t)) + 1/(Pe*dX)^2*(theta_stor(x,t)-2*theta_stor(x+1,t)+theta_stor(x+1,t)));
        end
       
    end
    
end

figure
hold on
plot(length, theta_HTF(:,1),'-b')
plot(length, theta_HTF(:,289),'--b')
plot(length, theta_stor(:,1),'-r')
plot(length, theta_stor(:, 289),'--r')
legend('HTF (t = 0)', 'HTF (t = 9 hours)', 'Stor (t = 0)', 'Stor (t = 9 hours)')
xlabel('Length (m)')
ylabel('Nondimensional temperature (T/T_s_t_o_r_,_i_n_i_t_i_a_l)')








