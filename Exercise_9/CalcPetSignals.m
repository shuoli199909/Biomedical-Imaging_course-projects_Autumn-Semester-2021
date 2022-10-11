%=========================================================================
%                                                                     
%	TITLE: 
%       CalcPetSignals.m				
%								
%	DESCRIPTION:						
%	    Calculate PET signal of analytical phantom
%
%	INPUT:								
%       radial coordinate, angle, phantom struct
%
%	OUTPUT:							
%       projection
%			
%	VERSION HISTORY:						
%	    1911016SK INITIAL VERSION 
%
%=========================================================================

%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================
function [signal] = CalcPetSignals(r,phi,phantom)
    
    signal  = zeros(size(r));
    phi     = phi/180*pi;
    
    sinphi  = sin(phi(:)); 
    cosphi  = cos(phi(:));
    
    rx      = r(:).*cosphi; 
    ry      = r(:).*sinphi;
    
    kk = find(phantom.ellipse(:,7));             % find organs with tracer               

    for k = kk'
    
        x0      = phantom.ellipse(k,1); y0 = phantom.ellipse(k,2);
        a       = phantom.ellipse(k,3); b  = phantom.ellipse(k,4);

        theta   = phantom.ellipse(k,5)*pi/180; 
        mue     = phantom.ellipse(k,6);

        r0      = [rx-x0,ry-y0]';

        % ----------------------------------------------------------------
        % Find entry and exit points of ellipse k (= tumor)
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
        % Sum projections
        % ----------------------------------------------------------------
        signal(i) = phantom.ellipse(k,7)*abs(s1-s2);               % [MBq]
        signal = reshape(signal,size(r));
    end
end


%=========================================================================
%=========================================================================     