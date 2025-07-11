% Evaluating epsilon from transverse cutting, and evaluating FDP distribution
% =====================================================================
% Inputs are 
% =====================================================================
% X to determine number of binding sites for every cross-section
% CSActive, binary array that determines which cross-sections are alive
% This information to ease the previous searching algorithm
% =====================================================================
% Outputs are
% ======================================================================
% FDPDist is an array that has the quantity of species of certain lengths
% extentLysis defines extent of lysis
% epsilon defines the voidage

function [FDP_dist,extentLysis,epsilon] =  clotPropertiesTransExtra(X, CSActive, lysisLevel,N_0)

    global epsilon_0 
	
    nCS = length(CSActive);
    threshold = N_0*lysisLevel; 
    
	% =================================================================
	% DETERMINING NUMBER OF CROSS-SECTIONS
	% =================================================================
	for it = 1 : nCS
	
		% If number of binding sites is less than the demanded threshold and cross-section is intact
		if (X(it + 3*nCS) < threshold(it) && CSActive(it) == 1)
		
			% Update the CSActive Array
			CSActive(it) = 0;
		end
	
	end
	
	% Total number of active cross-sections
	nActive = sum(CSActive);
	
	% =================================================================
	% FINDING THE DISTRIBUTION OF FDPs 
	% =================================================================
	
	FDP_dist = zeros(1,nCS + 1);
    
    count = 1;
    
    for it_dist = 1:nCS
        
        if CSActive(it_dist) == 1
            count = count + 1;
            
        else
            
            FDP_dist(count) = FDP_dist(count) + 1;
            count = 1;
            
        end
        
    end
	
	% =================================================================
	% DETERMINING EXTENT OF LYSIS
	% =================================================================

	% Extent of Lysis (defined as the number of fibrin monomers)
	
	extentLysis = 1  -(nCS - FDP_dist(1))/nCS;
    phi_f = (1 - epsilon_0)*(1-extentLysis);
	
	% =================================================================
	% FINDING EPSILON
	% =================================================================
	
	epsilon = 1 - phi_f;



end
