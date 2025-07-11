% ========================================================================
% STOCHASTIC TRANSVERSE CUTTING FIBRINOLYSIS 
% ========================================================================

function [t,ploti] = funHybridTransModelExtra(mu,sigma,rho_0,U_new,finalTime,nt,C_0,N_fibres,lysisLevel)

R_f0 = mu;

% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------

global U_0 check N_old extentLysis X

U_0 = U_new;
check = 0;
E_crit = 0.95;


% Simulation
startTime = 0; %s
delta_t = (finalTime - startTime)/nt;

hybridConstants;
%hybridConstantsExtra;

% Detrmining number of cross-sections
cord = dr;
constraint = 0;
it_check=0;
group = 1;
while constraint < rho_0/1000
    % Update the iteration 
    it_check = it_check + 1;
    % Assign random value 
    %R_f0(it_check) = round(floor(normrnd(mu,sigma))/cord)*cord;
    R_f0(it_check) = round(normrnd(mu,sigma));
    
    % Update constraint
    consStore = constraint;
        
    constraint = group*sum(R_f0.^2)*L_M*rho_fibre*pi()/(U_0*1E9);
   
    
    % Break the loop in case of error
    if it_check > 100000
       
        error('ERROR: Calculated number of cross-sections is too big. Not enough memory');
        
    end
end

N_CS = it_check

% Keeping the constraint of the fibrin density
Rf_end = sqrt((rho_0/1000-consStore)*U_0*1E9/(rho_fibre*pi()*L_M*group));
R_f0(end) = round(Rf_end/cord)*cord;

% Determining binding locations per cross-section
N_0 = zeros(1,N_CS);

for it_out = 1:N_CS

    for it = 1: round(R_f0(it_out)/dr)
        
        N_0(it_out) = N_0(it_out) + pi()/asin(dth/(2*it*dr));
        
    end

end


N_old = N_0(1);


Cs_tPA = C_0(1)*epsilon_0;
Cs_PLG = C_0(2)*epsilon_0;
Cs_PLS = C_0(3)*epsilon_0;

% -------------------------------------------------------------------------
% INITIALISATION
% -------------------------------------------------------------------------

N_species = 5; % n_tPA, n_PLG, n_PLS, n_tot, L_PLS
N_rxns = 9; % tPA (x2), PLG (x2), PLS (x2), M-M, cat, L_PLS dis

% Generating Stochiometry Matrix
%stochiometryMatrix;
stochiometryMatrix2;

% Define Initial State Array

X = zeros(1,N_CS*N_species + 1,'single');
% Order: n_tPA, n_PLG, n_PLS, n_tot, L_PLS

X(1:N_CS)            =  0;
X(N_CS   + 1:2*N_CS) =  0;
X(2*N_CS + 1:3*N_CS) =  0;
%X(3*N_CS + 1:4*N_CS) =  n_0*U_0*1E-6*1E-15*N_AV/N_CS;
X(3*N_CS + 1:4*N_CS) =  N_0;
X(4*N_CS + 1:5*N_CS) =  0;

% -------------------------------------------------------------------------
% ASSIGNING STORAGE
% -------------------------------------------------------------------------

% variables to use for RK gradient calculation
global store

store.ns_tPA = 0;
store.ns_PLG = 0;
store.ns_PLS = 0;
store.ns_tot = n_0;
store.L_PLS = 0;
store.epsilon = epsilon_0

% Free/ bound phases
ploti.CtPA = zeros(1,nt);
ploti.ntPA = zeros(1,nt);
ploti.CPLG = zeros(1,nt);
ploti.nPLG = zeros(1,nt);
ploti.CPLS = zeros(1,nt);
ploti.nPLS = zeros(1,nt);
ploti.LPLS = zeros(1,nt);
ploti.ntot = zeros(1,nt);

% Clot properties
ploti.epsilon = zeros(1,nt);
ploti.extentLysis = zeros(1,nt);
ploti.FDPDist = zeros(N_CS + 1,nt,'single');
CSActive = ones(1,N_CS,'single'); % Binary array that monitors which cross-sections are active

% Molecular Dynamics Properties
ploti.crossSection = zeros(N_CS,nt,'single');


% Storage initialisation
ploti.CtPA(1) = Cs_tPA;
ploti.CPLG(1) = Cs_PLG;
ploti.ntot(1) = n_0;
ploti.epsilon(1) = epsilon_0;

epsilon = epsilon_0;

tic();

I = zeros(1,N_CS,'single');

% -------------------------------------------------------------------------
% SIMULATION ENGINE
% -------------------------------------------------------------------------

for it = 2:nt
    % -------------------------------------------------------------------------
    % RUNGE-KUTTA SOLVER
    % -------------------------------------------------------------------------
    % To solve for the free phase concentrations
    it
    C_input = [Cs_tPA; Cs_PLG; Cs_PLS];
    current_t = (it-1)*delta_t;
    
    % NOTE: only file that needs to be changed is gradient file! The solver
    % is stand-alone. You, need a for loop however to run it for many time
    % steps.
    
    sol = rungeKuttaFourth(C_input, current_t, delta_t);
    
    Cs_tPA = sol(1);
    Cs_PLG = sol(2);
    Cs_PLS = sol(3);
    
    % -------------------------------------------------------------------------
    % GILLESPIE SOLVER
    % -------------------------------------------------------------------------
    
    % NOTE: only file that needs to be changed is the propensity function.
    % Inputs are state array, stochiometry matrix, duration of time
    % and propensity function.
    
    [X_nx,Cs_PLG,I] = gillespieSolver(X,v, delta_t,@propensityFunction, Cs_tPA, Cs_PLG, Cs_PLS, epsilon,I);
    
    X = X_nx;
    
    % -------------------------------------------------------------------------
    % EVALUATING LYSIS PROPERTIES
    % -------------------------------------------------------------------------
    
    [FDP_dist,extentLysis,epsilon] =  clotPropertiesTransExtra(X, CSActive,lysisLevel,N_0);
    
    % -------------------------------------------------------------------------
    % ASSIGNING OUTPUT
    % -------------------------------------------------------------------------
    
    store.ns_tPA = sum(X(1:N_CS))/(U_0*1E-6*1E-15*N_AV);
    store.ns_PLG = sum(X(1*N_CS + 1:2*N_CS))/(U_0*1E-6*1E-15*N_AV);
    store.ns_PLS = sum(X(2*N_CS + 1:3*N_CS))/(U_0*1E-6*1E-15*N_AV);
    store.ns_tot = sum(X(3*N_CS + 1:4*N_CS))/(U_0*1E-6*1E-15*N_AV);
    store.L_PLS  = sum(X(4*N_CS + 1:5*N_CS))/(U_0*1E-6*1E-15*N_AV);
    store.epsilon = epsilon;
    
    % -------------------------------------------------------------------------
    % STORE VALUES
    % -------------------------------------------------------------------------
    
    ploti.CtPA(it) = Cs_tPA;
    ploti.ntPA(it) = store.ns_tPA;
    
    ploti.CPLG(it) = Cs_PLG;
    ploti.nPLG(it) = store.ns_PLG;
    
    ploti.CPLS(it) = Cs_PLS;
    ploti.nPLS(it) = store.ns_PLS;
    ploti.LPLS(it) = store.L_PLS;
    
    ploti.ntot(it) = store.ns_tot;
    ploti.epsilon(it) = epsilon;
    ploti.extentLysis(it) = extentLysis;
    ploti.FDPDist(:,it) = FDP_dist;
    
    ploti.crossSection(:,it) = X(3*N_CS + 1:4*N_CS);
    
    if extentLysis > E_crit
       check = 1; 
    end
    
    
end

t = linspace(startTime, finalTime, nt);

end



