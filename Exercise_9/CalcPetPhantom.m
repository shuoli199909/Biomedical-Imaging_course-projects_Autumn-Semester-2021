%=========================================================================
%                                                                     
%	TITLE: 
%       CalcPetPhantom.m				
%								
%	DESCRIPTION:						
%	    Calculate discrete PET image from set of ellipses
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
function [image] = CalcPetPhantom(matrix,phantom)
    
    [x,y] = meshgrid(-fix(matrix/2):+fix(matrix/2));   
    image = zeros(size(x));
    
    kk = find(phantom.ellipse(:,7));             % find organs with tracer               
      
    for k = kk'
        
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
        % Assign number of nuclear disintegrations per second 
        % ----------------------------------------------------------------
        image(i) = phantom.ellipse(k,7);                        % [MBq]
        
    end
end


%=========================================================================
%=========================================================================     