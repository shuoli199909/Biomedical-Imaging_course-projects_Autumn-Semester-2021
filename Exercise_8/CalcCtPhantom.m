%=========================================================================
%                                                                     
%	TITLE: 
%       CalcCtPhantom.m				
%								
%	DESCRIPTION:						
%	    Calculate discrete CT image from set of ellipses
%
%	INPUT:								
%       image coordinates, phantom struct, organ index		
%
%	OUTPUT:							
%       image
%			
%	VERSION HISTORY:						
%	    191101SK INITIAL VERSION
%
%=========================================================================

%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================f
function [image] = CalcCtPhantom(matrix,phantom,z)
    
    [x,y] = meshgrid(-fix(matrix/2):+fix(matrix/2));   
 
    image = zeros(size(x));
    
    for k = z
        
        theta   = phantom.ellipse(k,5)*pi/180;
        
        X0      = [x(:)'-phantom.ellipse(k,1);y(:)'-phantom.ellipse(k,2)];
        D       = [1/phantom.ellipse(k,3) 0;0 1/phantom.ellipse(k,4)];
        Q       = [cos(theta) sin(theta); -sin(theta) cos(theta)];
         
        % ----------------------------------------------------------------
        % Find inside of ellipse given by X0,D,Q
        % ----------------------------------------------------------------
        equ = sum((D*Q*X0).^2);
        i = find(equ <= 1);
         
        % ----------------------------------------------------------------
        % Assign linear attenuation coefficients
        % ----------------------------------------------------------------
        image(i) = image(i) + phantom.ellipse(k,6);
        
    end
end


%=========================================================================
%=========================================================================     