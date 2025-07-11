% =======================================================================
% COMPARING RESULTS, SHRINK MODEL HYBRID V CONTINUOUS
% =======================================================================
% 
clear
clc
global U_0 check N_old extentLysis
global k_a k_r k_2 k_cat K_M epsilon_0 n_0 LtVt N_CS N_AV N_0
global k_a k_r store
global N_CS epsilon_0 N_species
global check
global k_a k_r k_2 k_cat K_M epsilon_0 n_0 LtVt N_CS N_AV N_0 count X


% % Folder has to be in the same directory as this file
% %addpath('\\icnas4.cc.ic.ac.uk\ap4409\MATLAB\stochKinPaperSims\stochasticKineticsShrink');
% %addpath('\\icnas4.cc.ic.ac.uk\ap4409\MATLAB\stochKinPaperSims\continuousKinetics');



% =======================================================================
% INPUTS
% =======================================================================

rho_0 = 3;
U_0 = 100;
finalTime = 10*60;
nt = 1000;
Cs_tPA = 0.04; Cs_PLG = 2.2; Cs_PLS = 0;
lysisLevel = 0.2;
N_fibres = 1;
R_f0 = 100;

% 
sigma_array = 50;
num_sim = 1;

for it_sim = 1:num_sim

    for it = 1 : length(sigma_array)

        sigma = sigma_array(it);

        % Running stochasic approach
        % All variables located in ploti class
        fprintf('Running Hybrid model... \n')
        tic()
       % [t,ploti] = funHybridTransModel(R_f0,rho_0,U_0,finalTime,nt,[Cs_tPA, Cs_PLG, Cs_PLS],N_fibres,lysisLevel);
        [t,ploti] = funHybridTransModelExtra(R_f0,sigma,rho_0,U_0,finalTime,nt,[Cs_tPA, Cs_PLG, Cs_PLS],N_fibres,lysisLevel);
        toc()

        save(strcat(num2str(sigma),'_SigmaComparisonRun',num2str(it_sim)));
        plot(t,ploti.extentLysis,'-k')
        hold on

    end

end

% For the first part (Fibrinolysis as the transverse cutting of fibres)
%plotAllNow;

% For the second part (Fine fibres are lysed slower than coarse clots)
%plotMulti

% For the third part (Changing fibrin distribution)
%plotDist

% For the fourth part (Cellular components)

    % Looking through alpha 
% alpha_list = [1000,100,10,1]; % Real values is divided by 1000
% phi_p0 = 0.4;
% U_0 = 1000;
% 
% for it = 1: length(alpha_list)
% 
%     alpha = alpha_list(it)/1000
%     tic()
%     [t,ploti] = funHybridTransModelExtra2(R_f0,0,rho_0,U_0,finalTime,nt,[Cs_tPA, Cs_PLG, Cs_PLS], N_fibres,lysisLevel,alpha,phi_p0);
%     toc()
%     [k_f,k_t]  = calculate_permeability(ploti.phi_p,ploti.phi_f,1,100);
% 
%     save(strcat('alphaValue_',string(alpha_list(it))));
% %     
% %     plot(t,ploti.phi_p)
% %     hold all
% 
%     %semilogy(t,1./k_t)
%     %hold all
% end
% 
% 
% 
%     % Looking through different retraction
% 
% R_list = [50,60,70,80,90,100];
% phi_f0 = 3/280;
% U_0 = 100;
% 
% 
% for item = 1: length(R_list)
% 
%     R = R_list(item)/100
%     C = (R - (phi_p0 + (1-phi_p0)*phi_f0))/((1-phi_p0)*(1-phi_f0));
%     eps_new = C*(1-phi_f0)/(phi_f0 + C*(1-phi_f0));
%     rho_new = 280*(1-eps_new)
%     finalTime = 60*60;
%     alpha = 0.1;
%     nt = 100000;
% 
%     tic()
%     [t,ploti] = funHybridTransModelExtra2(R_f0,0,rho_new,U_0,finalTime,nt,[Cs_tPA, Cs_PLG, Cs_PLS], N_fibres,lysisLevel,alpha,phi_p0);
%     toc()
%     [k_f,k_t]  = calculate_permeability(ploti.phi_p,ploti.phi_f,1,100);
% 
%     plot(t,1./k_t)
%     hold all
% 
%     save(strcat('retrac_',string(R_list(item))));
% 
% end

%plotCell;





