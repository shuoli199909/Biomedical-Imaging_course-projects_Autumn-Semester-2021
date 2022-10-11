%=========================================================================
%                                                                     
%	TITLE: 
%       CalcCtSinogram.m				
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
function [sino] = CalcCtSinogram(angles,matrix,phantom)

    sino = zeros(matrix+1,length(angles)); idx = 1;  
        
    % --------------------------------------------------------------------
    % Loop over CT projection angles
    % --------------------------------------------------------------------
    for phi=angles
         
        % ----------------------------------------------------------------
        % Calculate CT line integrals 
        % ----------------------------------------------------------------
        [r,phi_] = meshgrid(-fix(matrix/2):1:+fix(matrix/2),phi);  
        sino(:,idx) = CalcLineIntegrals(r,phi_,phantom)/(matrix+1);
        
        idx = idx+1;
    end
end
    
    
%=========================================================================
%=========================================================================     