%=========================================================================
%                                                                     
%	TITLE: 
%       CalcPetSinogram.m				
%								
%	DESCRIPTION:						
%	    Calculate sinogram
%
%	INPUT:								
%       projection angles, matrix size, phantom
%
%	OUTPUT:							
%       sinogram
%			
%	VERSION HISTORY:						
%	    191101 INITIAL VERSION 
%
%=========================================================================

%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================
function [sino] = CalcPetSinogram(angles,matrix,phantom)

    sino = zeros(matrix+1,length(angles)); idx = 1;  
        
    % --------------------------------------------------------------------
    % Loop over PET projection angles
    % --------------------------------------------------------------------
    for phi=angles
         
        % ----------------------------------------------------------------
        % Calculate PET line integrals without attenuation
        % ----------------------------------------------------------------
        [r,phi_] = meshgrid(-fix(matrix/2):1:+fix(matrix/2),phi);  
        sino(:,idx) = CalcPetSignals(r,phi_,phantom)/(matrix+1);
        
        % ----------------------------------------------------------------
        % Calculate PET sinogram with attenuation
        % ----------------------------------------------------------------
        attsino = CalcLineIntegrals(r,phi_,phantom)/(matrix+1);
        sino(:,idx) = sino(:,idx).*exp(-attsino');

        idx = idx+1;
    end
end
    
    
%=========================================================================
%=========================================================================     