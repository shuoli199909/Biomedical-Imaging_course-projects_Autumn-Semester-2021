%=========================================================================
%                                                                     
%	TITLE: 
%       CalcFBPRecon.m				
%								
%	DESCRIPTION:						
%	    Calculate FBP reconstruction
%
%	INPUT:								
%       projection angles, matrix size, sino, filter
%
%	OUTPUT:							
%       reconstructed image
%			
%	VERSION HISTORY:						
%	    191101 INITIAL VERSION 
%
%=========================================================================
%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================
function [recon] = CalcFBPRecon(angles,matrix,sino,filter)

    idx    = 1; 
    recon  = zeros(matrix+1,matrix+1);
    
    % --------------------------------------------------------------------
    % Generate indices of static x,y coordinate system
    % --------------------------------------------------------------------
    [x,y]  = meshgrid(-fix(matrix/2):+fix(matrix/2));   
 
    for phi=angles
        
        % ----------------------------------------------------------------
        % Compute indices of rotating r,s coordinate system
        % ----------------------------------------------------------------
        rs = round(x*cos(phi/180*pi)+y*sin(phi/180*pi));
        rs = rs+ceil((matrix+1)/2);
        ix = find((rs>=1)&(rs<=(matrix+1)));
        
        % ----------------------------------------------------------------
        % Convolve projection with filter
        % ----------------------------------------------------------------
        filteredprojection = conv(sino(:,idx),filter,'same');  
        
        % ----------------------------------------------------------------
        % Back-project filtered projection
        % ----------------------------------------------------------------
        recon(ix) = recon(ix)+filteredprojection(rs(ix));      
           
        idx = idx+1;
    end
    
    % ----------------------------------------------------------------
    % Normalize for number of projections
    % ----------------------------------------------------------------
    recon = recon/length(angles);
    
end


%=========================================================================
%=========================================================================     