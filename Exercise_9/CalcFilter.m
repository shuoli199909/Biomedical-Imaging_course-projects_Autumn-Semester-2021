%=========================================================================
%                                                                     
%	TITLE: 
%       CalcCtFilter.m				
%								
%	DESCRIPTION:						
%	    Calculate high-pass |u|
%
%	INPUT:								
%       matrix size
%
%	OUTPUT:							
%       filter
%			
%	VERSION HISTORY:						
%	    120216SK INITIAL VERSION 
%	    191020SK UPDATE
%
%=========================================================================

%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================f
function [filter] = CalcCFilter(matrix)

    % --------------------------------------------------------------------
    % Define high-pass filter according to |u|
    % --------------------------------------------------------------------
    filter = abs(-fix(matrix/2):+fix(matrix/2));
       
    % --------------------------------------------------------------------
    % Define low-pass filter according to 0.5+0.5*cos(u)
    % --------------------------------------------------------------------
    weight = 0.5+0.5*cos(filter/matrix*2*pi); 
      
    % --------------------------------------------------------------------
    % Normalize power and transform filter into spatial domain using U2R()
    % --------------------------------------------------------------------
    filter = U2R(3*filter.*weight)+0.003;
        
end


%=========================================================================
%=========================================================================     