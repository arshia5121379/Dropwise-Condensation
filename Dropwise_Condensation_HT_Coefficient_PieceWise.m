clc;
clear;
addpath('C:\Program Files (x86)\REFPROP');

%% Definition of Dimensionless Numbers : 
% Re = (rho * u * L) / mu;
% Sh = (h * L) / D;
% Sc = nu / D;
% Le = alpha / D;
% Pr = nu / alpha;
% Nu = (h * L) / k;

%% INPUTS 
T_W   = 5.0 + 273.15;   % Temp. of wall [k]
T_B   = 25.0 + 273.15;  % Temp. of bulk of inflow [k]
T_exit= T_B - 2;        % Temp. of exit [k] 
RH    = 1.0;            % relative humidty [-]
D_h   = 0.0158;         % hydrolic diameter [m]
P_amb = 1e5;            % pressure of atmospehre 
theta = 150*pi/180;     % contact angle [rad]
theta_A= 155*pi/180;    % advancing contact angle [rad] (theta_max)
Diff_droplet= 2.9482e-5;% Diff. Coef. of droplet  [m^2/s]
g = 9.81;               % grvity acceleration [m/s^2]
r_d_mean = 0.5e-3;      % mean droplet's radious [m]
x_total = 0.1;          % total condensation length [m]
dx = 0.005;              % elemental length [m]
Alpha = 90*pi/180;      % condensation surface inclination angle [rad]
total_HT_area = 0.062368; % total area of heat transfer [m^2]
A1 = pi*(3.5e-2)^2;     % [m^2]
A2 = 0.0024261;         % [m^2] cross sectional area
U1 = 2.5;               % [m/s]
U2 = A1*U1/A2;  
U_in = U2;              

%% Properties at Bulk Temp.

h_f = refpropm('H','T',T_B,'Q',0,'WATER')*1e-3;
h_g = refpropm('H','T',T_B,'Q',1,'WATER')*1e-3;
h_fg = h_g - h_f;

Alpha_AIR = refpropm('%','T',T_B,'P',P_amb/1000,'AIR.MIX')*1e-4; % thermal Diff. of air 
Alpha_WATER = refpropm('%','T',T_B,'Q',0,'WATER')*1e-4;     % thermal Diff. of water


cp_WATER = refpropm('C','T',T_B,'Q',0,'WATER')*1e-3;
cp_VAPOR = refpropm('C','T',T_B,'Q',1,'WATER')*1e-3;
cp_AIR   = refpropm('C','T',T_B,'P',P_amb/1000,'AIR.MIX')*1e-3; 

k_WATER = refpropm('L','T',T_B,'Q',0,'WATER'); % thermal conductivity of water
k_VAPOR = refpropm('L','T',T_B,'Q',1,'WATER'); % thermal conductivity of steam
k_AIR   = refpropm('L','T',T_B,'P',P_amb/1000,'AIR.MIX');

sigma_WATER = refpropm('I','T',T_B,'Q',0,'WATER'); % surface tension

Rho_WATER   = refpropm('D','T',T_B,'Q',0,'WATER'); % Density of Water
Rho_VAPOR   = refpropm('D','T',T_B,'Q',1,'WATER'); % Density of Steam 
Rho_AIR     = refpropm('D','T',T_B,'P',P_amb/1000,'AIR.MIX');  % Densisty of Air

Nu_WATER = refpropm('$','T',T_B,'Q',0,'WATER')*1e-4; % kinematic viscosity of water
Nu_VAPOR = refpropm('$','T',T_B,'Q',1,'WATER')*1e-4; % kinematic viscosity of water
Nu_AIR   = refpropm('$','T',T_B,'P',P_amb/1000,'AIR.MIX')*1e-4;

Mu_AIR   = Nu_AIR*Rho_AIR;                           % dynamic viscosity of Air
Mu_VAPOR = Nu_VAPOR*Rho_VAPOR;                       % dynamic viscosity of vapor 

M_WATER = refpropm('M','T',T_B,'Q',0,'WATER');             % molar mass of water [g/mol]
M_AIR   = refpropm('M','T',T_B,'P',P_amb/1000,'AIR.MIX');  % molar mass of air

%% calculation of initial parameters 

Pg_Tsat = refpropm('P','T',T_B,'Q',1,'WATER')*1e3; % Pg at T saturation in bulk flow
P_vapor = RH * Pg_Tsat; 
P_air   = P_amb-P_vapor;   

Y_AIR = P_air/P_amb; % mol fraction of air in bulk
Y_H2O = 1 - Y_AIR;   % mol fraction of water in bulk

M_total = Y_H2O*M_WATER + Y_AIR*M_AIR; 

X_AIR = Y_AIR*M_AIR/M_total; % mass fraction of air in bulk 
X_H2O = 1 - X_AIR;           % mass fraction of water in bulk 

%% Initialize variables for piecewise calculation
n_segments = round(x_total / dx);
T_B_current = T_B; % Start with initial bulk temperature
q_i = zeros(n_segments, 1); % Heat transfer for each segment
T_out_i = zeros(n_segments, 1); % Outlet temperature for each segment

% Calculate elemental area
dA = total_HT_area / n_segments;

Rho_mix = Y_AIR*Rho_AIR + Y_H2O*Rho_VAPOR;
Mu_mix  = Y_AIR*Mu_AIR  + Y_H2O*Mu_VAPOR;


fprintf('Piecewise Heat Transfer Calculation:\n');
fprintf('Number of segments: %d\n', n_segments);
fprintf('Elemental length dx = %.4f m\n', dx);
fprintf('Elemental area dA = %.6f m²\n\n', dA);

%% Main loop for piecewise calculation
for i = 1:n_segments
    x_current = i * dx; % Current position along the channel
    
    % (step 1) initializing T_i
    T_init = (T_B_current + T_W)/2;
    eps = 0.00001;
    iter = 0;
    T_i_new = T_init; 
    T_ik = 280; 
    R_max = r_d_mean;

    Re_B = Rho_mix * U_in * D_h / Mu_mix; 

    while (abs(T_ik - T_i_new) > eps && iter < 100)
    % update T interface for next iteration 
    T_i_new = T_ik;

    % (step 2) finding mass fraction of air in interface layer (saturation line): 
    P_vapor_i = refpropm('P','T',T_i_new,'Q',1,'WATER')*1e3; % pressure of vapor at T_i 

    Y_vapor_i = P_vapor_i/P_amb; % mol fraction of vapor at interface flow; 
    Y_AIR_i   = 1 - Y_vapor_i;   % mol fraction of air at interface flow; 
    X_vapor_i = Y_vapor_i * M_WATER/M_total; % mass fraction of vapor at bulk flow; 
    
    X_AIR_i = 1 - X_vapor_i; % mass fraction of air in interface layer

    % (step 3) 
    % A : finding h_d using Eq. (2)
    B_i = log(X_AIR_i/X_AIR);
    h_d = (0.347*Re_B)/(1+1.75*exp(-330.5*B_i)); % h_d = (K/L)_droplet
    
    % B : finding h_g using Eq. (18) first finding Nu_g = C*Nu_film
    % first finding R_max for calculating C (surface enhancement factor)
    % R_max : 

    h_f_i = refpropm('H','T',T_i_new,'Q',0,'WATER')*1e-3;
    h_g_i = refpropm('H','T',T_i_new,'Q',1,'WATER')*1e-3;
    h_fg_i= h_g_i - h_f_i; 

    Nu_VAPOR_i = refpropm('$','T',T_i_new,'Q',1,'WATER')*1e-4;
    Nu_AIR_i   = refpropm('$','T',T_i_new,'P',P_amb/1000,'AIR.MIX')*1e-4;
    
    Nu_mix     = Y_AIR_i*Nu_AIR_i + Y_vapor_i*Nu_VAPOR_i; 

    Alpha_AIR_i   = refpropm('%','T',T_i_new,'P',(P_amb-P_vapor_i)/1000,'AIR.MIX')*1e-4; % thermal Diff. of air within interface layer
    Alpha_VAPOR_i = refpropm('%','T',T_i_new,'Q',1,'WATER')*1e-4; % thermal Diff. of vapor within interface layer

    Alpha_mix     = Y_AIR_i*Alpha_AIR_i + Y_vapor_i*Alpha_VAPOR_i; % thermal Diff. of mixture whithin interface layer     
    cp_mix        = Y_AIR_i*cp_AIR + Y_vapor_i*cp_VAPOR;

    Bo     = (Rho_WATER*g*(r_d_mean*2)^2*sin(Alpha))/sigma_WATER; % Bond number duo to Eq.(24)
    beta   = 1 + 0.096*Bo; % assuming sperichal droplet and D/w = 1 then beta = 1
    Re_D   = Rho_mix*U_in*(r_d_mean*2)/Mu_mix; % Reynolds number based on the droplet radius
    C_drag = 0.28 + 6/sqrt(Re_D) + 21/Re_D; % Drag Coeff. over the droplet

    theta_min = theta_A * (0.01*Bo^2-0.155*Bo+0.97);
    theta_max = theta_A;
    theta_S = pi - theta_min; % theta star for calculating Aproj

    A1 = 2*(cos(theta_max)-cos(theta_min))/pi^3;
    A2 = -3*(cos(theta_max)-cos(theta_min))/pi^2;
    A3 = cos(theta_max);        

    Dd = r_d_mean;

    Lf = (sin(theta_max)*(1-cos(theta_min)))/(sin(theta_min)*(1-cos(theta_max)));
    L1 = 2*Dd*Lf/(1+Lf);
    L2 = 2*Dd/(1+Lf);
    
    % A_project = 0.5*L1^2/(sin(theta_S))^2*theta_min + ...
    %     0.5*L1^2/tan(theta_S) + 0.5*L2^2/(sin(theta_min))^2*theta_min - ...
    %     0.5*L2^2/tan(theta_min) ; %  A_proj is the drag coefficient and the projected area of the droplet.

    a = 2/3*pi*Rho_WATER*g;
    b = 0.5*Rho_WATER*U_in^2*C_drag*1/(1+Lf)^2* ...
        ((Lf^2)*theta_min/(sin(beta)^2) + Lf^2/sin(beta) + ...
        theta_min/(sin(theta_min))^2 - 1/(tan(theta_min))); 
    c = -sigma_WATER*(12*A1*pi^2 + 4*A2*pi);

    coeff = [a^2 0 0 0 b^2 0 0 0 -c^2];
    solutions = roots(coeff);
        
    R_max = min(abs(solutions));

    % calculating Correction Coeff. in table 2 
    %Alpha_mix = Y_AIR_i*Alpha_AIR_i + Y_vapor_i*Alpha_VAPOR_i
      
    % dimensionless numbers for mixture properties in interface layer
    ksi = (1.8*log10(Re_B)-1.5)^(-2);
    Sc = Nu_AIR / Diff_droplet;    
    Pr = Nu_AIR/Alpha_AIR;
    Ja = cp_AIR*(T_B_current-T_W)/h_fg;

    
    Sh0 = (ksi/8 *Re_B*Sc )/(1+12.7*sqrt(ksi/8)*(Sc^(2/3)-1))*(1+(D_h/x_current)^(2/3));   %Debug    
        
    Nu0 = (ksi/8 *Re_B*Pr )/(1+12.7*sqrt(ksi/8)*(Pr^(2/3)-1))*(1+(D_h/x_current)^(2/3));   %Debug
        
    %Sh0 = 47.8;
    %Nu0 = 57.5;

    Le = (Nu0/Sh0)^(1/3) ;         %  Alpha_VAPOR_i/D_droplet; % Lewis number 

    % suction Correction Coeff. for heat transfer
    R_t = (cp_VAPOR*Sh0)/(cp_mix*Le*Nu0)*log((1-X_H2O)/(1-X_vapor_i));
    SIGMA_sq = R_t/(exp(R_t)-1);
    
    % Fog Correction Coeff. for heat transfer
    dP_dT_i = refpropm('E','T',T_i_new,'Q',0,'WATER')*1e3; % pressure gradient along the saturation line  
    dF_dT   = dP_dT_i/P_amb;  % duo to reference "Effect of Fog Formation on
                              % Turbulent Vapor Condensation With
                              % Noncondensable Gases" Table 2

    SIGMA_fq= (1+((h_fg/cp_WATER)*(1/Le)* ( (X_H2O-X_vapor_i)/(T_B_current-T_i_new))) * (Sh0/Nu0)) / ( 1+(h_fg/cp_WATER)*(1/Le)*dF_dT);

    % suction Correction Coeff. for mass transfer
    R_w = (X_H2O-X_vapor_i)/(1-X_vapor_i);
    SIGMA_sm = log( (1-R_w))/(-R_w);
    % Fog Correction Coeff. for mass transfer
    SIGMA_fm= (1+ ((h_fg/cp_mix)*(1/Le)* ( (X_H2O-X_vapor_i)/(T_B_current-T_i_new))* (Sh0/Nu0))^(-1) ) /( 1+((h_fg/cp_mix)*(1/Le)*dF_dT)^(-1));

    % calculating C_factor from Eq. (37)
    % find constants m & n from table 4 
    m = 1.451; 
    n = 0.933; 
    C_Film_Drop =m*R_max+n ;
        
        % calculating Nu_f from Eq. 18 : 
        Nu_g = C_Film_Drop*( Nu0*SIGMA_sq*SIGMA_fq + ...
            Sh0*SIGMA_sm*SIGMA_fm*Pr/Sc*1/Ja*log(X_AIR_i/X_AIR)); 
        k_mix = k_AIR*Y_AIR_i + k_VAPOR*Y_vapor_i;
        hg = Nu_g*k_mix/D_h;
        
        % calculation of convetion (air to dry surface without considering droplets) Nu using eq(8-60) incropera: 
        Nu_convection = 0.0265*Re_B^(4/5)*Pr^(0.3); 
        h_convection  = Nu_convection*k_mix/D_h;  

        %fully developed (hydrodynamicall and thermally) turbulent flow in a smooth circular tube
        h_t = (1/hg + 1/h_d)^(-1) ;  % rectified heat transfer Coeff.
        h_Total = (1/hg + 1/h_d)^(-1) + h_convection;

        q_t     = h_t*(T_B_current-T_W);
        q_Total = h_Total*(T_B_current-T_W);

        q_convetion = h_convection*(T_B_current - T_W);
        q_g = q_Total - q_convetion; 
        q_d = h_d * (T_ik - T_W);
        
        T_ik = T_B_current - (q_g)/hg;
        iter = iter + 1;
    end 

    % Calculate heat transfer for this segment
    q_i(i) = q_Total * dA; % Heat transfer rate [W] for this segment
    
    % Update bulk temperature for next segment
    m_dot_a = Rho_mix * U_in * A2; % Mass flow rate [kg/s] 
    T_out_i(i) = T_B_current - (q_Total * dA) / (m_dot_a * cp_mix * 1000);
    
    fprintf('Segment %d (x = %.3f m):\n', i, x_current);
    fprintf('  Inlet T_B = %.2f K\n', T_B_current);
    fprintf('  Outlet T_B = %.2f K\n', T_out_i(i));
    fprintf('  Heat flux q_Total = %.2f W/m²\n', q_Total);
    fprintf('  Heat transfer Q_i = %.4f W\n\n', q_i(i));
    
    % Update bulk temperature for next segment
    T_B_current = T_out_i(i);
end

%% Display results
total_heat_transfer = sum(q_i);
fprintf('\n=== FINAL RESULTS ===\n');
fprintf('Total heat transfer: %.4f W\n', total_heat_transfer);
fprintf('Final outlet temperature: %.2f K\n', T_out_i(end));
fprintf('Temperature drop: %.2f K\n', T_B - T_out_i(end));

% Display all q_i values
fprintf('\nHeat transfer for each segment:\n');
for i = 1:n_segments
    fprintf('Segment %d: Q = %.4f W\n', i, q_i(i));
end

% Plot results
figure;
subplot(2,1,1);
plot((1:n_segments)*dx, q_i, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('Position along channel (m)');
ylabel('Heat transfer (W)');
title('Heat Transfer vs. Position');
grid on;

subplot(2,1,2);
plot((0:n_segments)*dx, [T_B; T_out_i], 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'r');
xlabel('Position along channel (m)');
ylabel('Bulk Temperature (K)');
title('Bulk Temperature vs. Position');
grid on;