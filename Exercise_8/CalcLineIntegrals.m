%=========================================================================
%                                                                     
%	TITLE: 
%       CalcLineIntegrals.m				
%								
%	DESCRIPTION:						
%	    Calculate projection from analytical phantom
%
%	INPUT:								
%       radial coordinate, angle, phantom struct, anode voltage	
%
%	OUTPUT:							
%       projection
%			
%	VERSION HISTORY:						
%	    120216SK INITIAL VERSION 
%	    191020SK UPDATE
%
%=========================================================================

%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================f
function [projection] = CalcLineIntegrals(r,phi,phantom,ua)
    
    projection  = zeros(size(r));
    phi         = phi/180*pi;
    
    sinphi  = sin(phi(:)); 
    cosphi  = cos(phi(:));
    
    rx      = r(:).*cosphi; 
    ry      = r(:).*sinphi;
    
    for k=1:length(phantom.ellipse(:,1))
        
        x0      = phantom.ellipse(k,1); y0 = phantom.ellipse(k,2);
        a       = phantom.ellipse(k,3); b  = phantom.ellipse(k,4);
        
        theta   = phantom.ellipse(k,5)*pi/180; 
        mue     = phantom.ellipse(k,6);
        
        r0      = [rx-x0,ry-y0]';
        
        % ----------------------------------------------------------------
        % Find entry and exit points of each ellipse
        % ----------------------------------------------------------------
        DQ      = [cos(theta)/a sin(theta)/a; -sin(theta)/b cos(theta)/b];
        DQphi   = DQ*[sinphi,-cosphi]'; 
        DQr0    = DQ*r0;
        
        A       = sum(DQphi.^2); 
        B       = 2*sum(DQphi.*DQr0);
        C       = sum(DQr0.^2)-1; 
        equ     = B.^2-4*A.*C;
        
        i       = find(equ>0);
        s1      = 0.5*(-B(i)+sqrt(equ(i)))./A(i);
        s2      = 0.5*(-B(i)-sqrt(equ(i)))./A(i);
        
        % ----------------------------------------------------------------
        % Compute projection i.e. mue * length inside object
        % ----------------------------------------------------------------
        proj    = mue*abs(s1-s2);
        
        % ----------------------------------------------------------------
        % Sum projections
        % ----------------------------------------------------------------
        projection(i) = projection(i)+proj;
    end
    
    projection = reshape(projection,size(r));
    
end


%=========================================================================
%=========================================================================     