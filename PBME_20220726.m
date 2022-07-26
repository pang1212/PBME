clear
%   Equilibrium Rxns
%                   H   OH  CO2 HCO3     
Chemical(1,:)   = [ 1   1   0   0]; %H20  = H+ + OH-
Chemical(2,:)   = [ 1   0   -1  1]; %CO2 + H2O = H+ + HCO3-
pKa = [14 6.37]';

K_H=29.4; % Henry's law constant for CO2(g)<-> CO2(aq)    [=] atm/M

% Initial Conditions
Water           = [1    1   0   0]*1e-7;
Salts           = [0    0   1/K_H   3];   %Concentrations Added (mol/L)

%Component Properties
Ion_Charge      = [ 1   -1   0  -1];     

%Process Inputs
Rxns        = Chemical;
N_species   = size(Rxns,2);
N_stoic     = size(Chemical,1);
Initial_Conc    = Salts + Water;
Properties      = [Initial_Conc;Ion_Charge];
stoic = Rxns(1:N_stoic,:);
Ci = Properties(1,:);
M = rref(null(stoic)');

%Constants
F = 96485;

% Bulk Electrolyte Concentrations
C_b = Equilibrium_0D(stoic, pKa, Ci, M);

% Determine concentration of H+/OH- addition from HOR at anode or CO2RR/HER at the cathode:
iden = 0.0222491; % Total current density [=] A/cm2, average of 4 single cells in the 4-cell PBME
V = 0.01/60; % Volumetric flow rate of liquid [=] L/s
Ae = 4; % Electrode area [=] cm2
n = 1; % mole e-/mole H+ for HOR at the anode of PBME and mole e-/mole OH- for CO2RR/HER at the cathode
C_H_add = iden*Ae/(n*F)/V; % Concentration of H+ addition from HOR at the anode [=] mol/L
C_OH_add = C_H_add; % Concentration of OH- addition from CO2RR/HER at the cathode [=] mol/L

% Determine concentration of CO addition from CO2RR at the first cathode:
i_CO = 0.012885951; % CO partial current density [=] A/cm2, average of 4 single cells in the 4-cell PBME
n_CO = 2; % mole e-/mole CO produced at the cathode
% C_CO_gen = i_CO*Ae/(n_CO*F)/V

%PBME Loop
Ci = C_b;
pH_i = -log10(Ci(1));
C_CO2g_i = 0;
% C_CO2g_i = 18/1000;
C_CO2g_c = C_CO2g_i;
C_CO_g = 0;
C_aq_out_o = (Ci(3) + Ci(4));
C_aq_out_c = C_aq_out_o;
I0 = [pH_i Ci(3) Ci(4) C_CO2g_i C_CO_g 0]'; % pH C_CO2(aq) C_HCO3- C_CO2(g) CO C_CO2(g)_add
I = [];
I = [I,I0];
CO2_U = [];
CO2_utilization = 0;

number_of_cells = 25;

% Solve PBME
for i = 1:number_of_cells
    
    % Anodes
    Ci(1) = Ci(1) + C_H_add;
    C_a(i,:) = Equilibrium_0D(stoic, pKa, Ci, M);
    if C_a(i,3) > (1/K_H)
        C_a(3) = (1/K_H);
    else
    end
    pH_a = -log10(C_a(i,1));
    C_CO2aq_a = C_a(i,3);
    C_HCO3_a = C_a(i,4);
    C_aq_out_a = C_CO2aq_a + C_HCO3_a;
    C_CO2g_add_a = (C_aq_out_c- C_aq_out_a);
    C_CO2g_a = C_CO2g_c + C_CO2g_add_a;
    I1 = [pH_a C_CO2aq_a C_HCO3_a C_CO2g_a C_CO_g C_CO2g_add_a]';
    I = [I,I1];
    
    % Cathodes
    % Calculate i_CO based on dissolved CO2 concentration
    k_c = 12.886/1000/n_CO/F/(1/K_H); % k_c is mass transfer coefficient
    i_CO = n_CO*F*k_c*C_CO2aq_a;
    C_CO_gen = i_CO*Ae/n_CO/F/V;
    Ci = C_a(i,:);
    Ci(2) = Ci(2) + C_OH_add;
    Ci(3) = Ci(3) - C_CO_gen;
    C_c(i,:) = Equilibrium_0D(stoic, pKa, Ci, M);
    pH_c = -log10(C_c(i,1));
    C_CO2aq_c = C_c(i,3);
    C_HCO3_c = C_c(i,4);
    C_aq_out_c = C_CO2aq_c + C_HCO3_c;
    C_CO2g_add_c = (C_aq_out_a - C_aq_out_c) - C_CO_gen;
    C_CO2g_c = C_CO2g_a + C_CO2g_add_c;
    C_CO_g = C_CO_g + C_CO_gen;
    I2 = [pH_c C_CO2aq_c C_HCO3_c C_CO2g_c C_CO_g C_CO2g_add_c]';
    I = [I,I2];
    CO2_utilization = C_CO_g/(C_CO_g+C_CO2g_c)*100;
    CO2_U = [CO2_U,CO2_utilization];
    
end

figure(1)
plot(1:number_of_cells,CO2_U,'-o','color','k','linewidth',1,'MarkerFaceColor','k')
xlabel('Number of cells')
ylabel('CO_2 utilization / %')
xtickformat('%.0f'); ytickformat('%.0f');
set(gca,'fontsize',15); set(gca,'FontName','Arial');

figure(2)
plot(0:0.5:number_of_cells,I(1,:),'-o','color','k','linewidth',1,'MarkerFaceColor','k')
xlabel('Number of cells')
ylabel('pH')
% ylim([4 9])
xtickformat('%.0f'); ytickformat('%.1f');
set(gca,'fontsize',15); set(gca,'FontName','Arial');

function C = Equilibrium_0D(stoic, pKa, Ci, M)  
    N_Species   = size(stoic,2);
    N_Rxn       = size(stoic,1);
    
    if isempty(pKa)
        C = Ci;
        return
    end
    
    %Define equilibrium function
    fun2 = @(x)React_0D(x, stoic, pKa, Ci, M);

    pC_guess = -log10(Ci);
    
    %Solve System
    options = optimoptions('fsolve','Display','iter-detailed','StepTolerance',1e-10);
    pEquilibrium_Concentrations = fsolve(fun2,pC_guess,options);
    C = 10.^-pEquilibrium_Concentrations;
end


function Outs = React_0D(x, stoic, pKa, Ci, M)

    Potential_Diff_1 = (pKa - stoic*x')';
    Potential_Diff_2 = (M*(10.^-x' - Ci'))';
    Outs = [Potential_Diff_1, Potential_Diff_2];

end