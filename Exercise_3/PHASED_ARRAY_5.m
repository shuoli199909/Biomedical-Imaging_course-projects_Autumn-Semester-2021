%=========================================================================
%                                                                     
%       BIOMEDICAL IMAGING
%       ULTRASOUND 2
%
%=========================================================================

%=========================================================================
%	PHASED ARRAY
%=========================================================================

function [] = PHASED_ARRAY()

    clear all; close all; 
    
    fprintf ( '-----------------------------------------\n' );  
    fprintf ( ' PHASED ARRAY                            \n' );  
    fprintf ( '-----------------------------------------\n' );  
    
    
    % Ultrasound specs
    
    f = 1.5E6;                      % frequency [Hz]
    c = 1480;                       % speed of sound in water [m/s]
    k = 2*pi*f/c;                   % wave number [rad/m]
    lambda = c/f;                   % wavelength [m]
    
    
    % Transducer array
    
    % nt = 1;                         % number of transducers
    nt = 40;
    pitch = 0.5*lambda;             % spacing of transducers
    
    
    % Definition of 2D grid for wave calculation
    
    n1 = 4000;                       % grid size in each dimension
    n2 = 1024;
    x_range = [0,0.4];              % covered range in x [m]
    y_range = [-0.2,0.2];           % covered range in y [m]
    
    x_values = x_range(1)+(x_range(2)-x_range(1))*[0:(n1-1)]./(n1-1);    % vector of x values on the grid
    y_values = y_range(1)+(y_range(2)-y_range(1))*[0:(n2-1)]./(n2-1);    % vector of y values on the grid
    
    [x,y] = meshgrid(x_values,y_values);    % the grid, given by 2D arrays x,y   
      
    
    % Loop through transducers and sum up their wave contributions
    
    wave = zeros(n2,n1);                      % container for pressure wave
    
    radi = 0.05;
    
    for t = 1:nt  
        
        x0 = -0.001;                        % x coordinate of transducer, 1 mm outside the container 
        y0 = pitch*(t-(nt+1)/2);            % y coordinate of transducer
        
        r = sqrt((x-x0).^2+(y-y0).^2);      % distance from transducer, on the grid
               
        % phase = 0;                          % phase of transducer input [rad]
        phase = - 2*pi*radi/lambda + 2*pi*sqrt((radi/lambda)^2 - ((((nt+1)/2)-t)/2)^2);
        % disp(phase);
        amplitude = 1;                      % relative amplitude of transducer input
                     
        wave = wave + exp(i*(k*r))./sqrt(r) .* amplitude .* exp(i*phase);   % add partial wave to net wave
        
    end

    figure('position',[100 100 800 500])    
    imagesc(real(wave)), title('real part');  
    % set(gca, 'XTick', []);
    % set(gca, 'YTick', []);
    
    figure('position',[100 100 800 500])
    imagesc(abs(wave)), title('magnitude')
    % set(gca, 'XTick', []);
    % set(gca, 'YTick', []);    
    
end
    