%=========================================================================
%                                                                     
%	TITLE: 
%       XCT - EXERCISE 2
%								
%	DESCRIPTION:						
%       COMPUTE SINOGRAM AND CT RECONSTRUCTIONS (FBP,FT)
%
%	INPUT:								
%       NONE	
%
%	OUTPUT:							
%       DISPLAY
%			
%	VERSION HISTORY:						
%	    150818SK INITIAL VERSION
%       170926SK UPDATE
%       211023SK UPDATE
%
%=========================================================================

%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================
function [] = XCT_EXERCISE2()

    clear; close all; clc;
    
   
    % --------------------------------------------------------------------
    % Display title
    % -------------------------------------------------------------------- 
    fprintf ( '-----------------------------------------\n' );  
    fprintf ( ' BIOMEDICAL IMAGING - XCT-EXERCISE #2\n' );  
    fprintf ( '-----------------------------------------\n' );  
    
    
    % --------------------------------------------------------------------
    % Set imaging parameters
    % --------------------------------------------------------------------
    matrix      = 256;                      % image matrix  [1pixel = 1mm]
    ua          = 50;                       % anode voltage [keV]
    % --------------------------------------------------------------------
    % Set density
    % (see: Table 2 in www.nist.gov/pml/data/xraycoef/)
    % --------------------------------------------------------------------
    rho_blood       = 1.060;                % density blood    [g/cm3]
    rho_bone        = 1.920;                % density bone     [g/cm3]
    rho_lung        = 0.001;                % density lung/air [g/cm3]
    rho_muscle      = 1.050;                % density muscle   [g/cm3]
    
   
    % --------------------------------------------------------------------
    % Set X-ray mass attenuation coefficients for 50 and 150 keV
    % (see: Table 4 in www.nist.gov/pml/data/xraycoef/)
    % --------------------------------------------------------------------
    mac_blood(1)    = 0.228;                % blood @  50 keV   [cm2/g]
    mac_blood(2)    = 0.149;                % blood @ 150 keV   [cm2/g]
    
    mac_bone(1)     = 0.424;                % bone @  50 keV    [cm2/g]
    mac_bone(2)     = 0.148;                % bone @ 150 keV    [cm2/g]
    
    mac_lung(1)     = 0.208;                % lung @  50 keV    [cm2/g]
    mac_lung(2)     = 0.136;                % lung @ 150 keV    [cm2/g]
    
    mac_muscle(1)   = 0.226;                % muscle @  50 keV  [cm2/g]
    mac_muscle(2)   = 0.149;                % muscle @ 150 keV  [cm2/g]
    
    
    % --------------------------------------------------------------------
    % Calculate linear attenuation coefficients
    % --------------------------------------------------------------------
    mue_blood       = rho_blood.*mac_blood(:);    
    mue_bone        = rho_bone.*mac_bone(:); 
    mue_lung        = rho_lung.*mac_lung(:); 
    mue_muscle      = rho_muscle.*mac_muscle(:); 
    
    idx = 1; 
    if ua==150 
        idx = 2; 
    end        % set index
   
    
    % --------------------------------------------------------------------
    % Define analytical phantom using ellipses with [x0 y0 a b phi mue]
    %
    %       x0,y0   - center point [cm] (+x -> left-right, +y -> bottom-up)
    %       a,b     - half axes [cm]
    %       theta   - rotation angle relative to x-axis [deg]
    %       mue     - linear attenuation coefficient [1/cm]
    % --------------------------------------------------------------------
    phantom.ellipse = [   0   0   90  80  0   mue_muscle(idx);                  % thorax
                          0   0   70  60  0   mue_lung(idx)-mue_muscle(idx);    % lung
                       +110   0   15  15  0   mue_muscle(idx);                  % left arm muscle
                       +110   0    5   5  0   mue_bone(idx)-mue_muscle(idx);    % left arm bone
                       -110   0   15  15  0   mue_muscle(idx);                  % right arm muscle
                       -110   0    5   5  0   mue_bone(idx)-mue_muscle(idx);    % right arm bone
                          0   0   10  10  0   mue_blood(idx)-mue_lung(idx);     % aorta 
                        +30 +25   25  20 35   mue_muscle(idx)-mue_lung(idx)];   % heart
    
                    
    % --------------------------------------------------------------------
    % Compute phantom and display
    % --------------------------------------------------------------------
    [x,y]            = meshgrid(-fix(matrix/2):+fix(matrix/2));
    phantom.discrete = CalcDiscretePhantom(x,y,phantom,ua);

    DisplayData(phantom.discrete,[2,3,1]); title('Phantom'); 
    
    
    % --------------------------------------------------------------------
    % Compute and display projections and sinogram
    % --------------------------------------------------------------------
    projection_angles   = linspace(0,179,180); idx = 1;                                  
    phantom.sino        = zeros(matrix+1,length(projection_angles));
    
    for phi=projection_angles
        
        [r,phi_] = meshgrid(-fix(matrix/2):1:+fix(matrix/2),phi);  
        phantom.sino(:,idx) = CalcLineIntegrals(r,phi_,phantom,ua)/(matrix+1);
         
        DisplayData(phantom.sino(:,idx),[2,3,2]); 
        title(sprintf('Projection (phi=%d)',phi)); xlabel('r'); ylabel('mue(r)'); 
        
        DisplayData(phantom.sino,[2,3,3]); 
        title('Sinogram'); xlabel('phi'); ylabel('r'); drawnow;
        
        idx = idx+1;
    end
    
    
    % --------------------------------------------------------------------
    % Reconstruct image using simple back-projection (SBP)
    % --------------------------------------------------------------------
    [x,y]   = meshgrid(-fix(matrix/2):1:+fix(matrix/2));    % x,y coordinates
    [r,s]   = meshgrid(-fix(matrix/2):1:+fix(matrix/2));    % r,s coordinates
    idx     = 1;                                            % counter
    phantom.sbp  = zeros(matrix+1,matrix+1);                % image placeholder
    
    for phi=projection_angles
        
        % TASK 2.2 UNCOMMENT AND FILL IN HERE
        projection = phantom.sino(:,idx);
        projection = repmat(projection',matrix+1,1);
        r_ = x*cos(phi*pi/180) + y*sin(phi*pi/180);
        s_ = x*sin(phi*pi/180) - y*cos(phi*pi/180);
        projection_xy = griddata(r,s,projection,r_,s_);
        phantom.sbp = phantom.sbp+projection_xy;
        
        DisplayData(phantom.sbp,[2,3,4]); 
        title('Backprojection-with gridding'); xlabel('x'); ylabel('y'); drawnow;
        
        idx = idx+1;                                        % counter
    end  
    phantom.sbp(find(isnan(phantom.sbp))) = 0;
 
    % --------------------------------------------------------------------
    % Reconstruct image using Filtered Back-Projection (FBP)
    % --------------------------------------------------------------------
    [x,y]   = meshgrid(-fix(matrix/2):1:+fix(matrix/2));    % x,y coordinates
    [r,s]   = meshgrid(-fix(matrix/2):1:+fix(matrix/2));    % r,s coordinates
    idx     = 1;                                            % counter
    phantom.fbp = zeros(matrix+1,matrix+1);                 % image placeholder
    
    filter      = CalcFilter(matrix+1,0);   % TASK 2.3 EDIT FUNCTION
    
    for phi=projection_angles
        
        % TASK 2.3 UNCOMMENT AND FILL IN HERE
        projection = phantom.sino(:,idx);
        projection_fourier = R2U(projection');
        projection_fourier = projection_fourier.*filter;
        projection = U2R(projection_fourier);
        projection = repmat(projection,matrix+1,1);
        r_ = x*cos(phi*pi/180) + y*sin(phi*pi/180);
        s_ = x*sin(phi*pi/180) - y*cos(phi*pi/180);
        projection_xy = griddata(r,s,projection,r_,s_);
        phantom.fbp = phantom.fbp+projection_xy;
           
        DisplayData(phantom.fbp,[2,3,5]); 
        title('Filtered Backprojection-with gridding'); xlabel('x'); ylabel('y'); drawnow;
        
        idx = idx+1;                                        % counter
    end  
 
    phantom.fbp(find(isnan(phantom.fbp))) = 0;
    % --------------------------------------------------------------------
    % Reconstruct image using Fast Fourier Transform (FFT)
    % --------------------------------------------------------------------
    [p,q]   = meshgrid(-fix(matrix/2):1:+fix(matrix/2));    % p,q coordinates
    [u]     = meshgrid(-fix(matrix/2):1:+fix(matrix/2),1);  % u coordinate
    idx     = 1;                                            % counter
    phantom.fft = zeros(matrix+1,matrix+1);                 % placeholder
    p_ = [];
    q_ = [];
    projection_fourier = [];
    for phi=projection_angles
    
        % TASK 2.4 FILL IN HERE
        projection = phantom.sino(:,idx);
        p_ = [p_, u*cos(phi*pi/180)];
        q_ = [q_, u*sin(phi*pi/180)];
        projection_fourier = [projection_fourier, R2U(projection')];
        idx = idx+1;                                        % counter
    end
    
    % Grid data onto Cartesian matrix
    % TASK 2.4 UNCOMMENT AND FILL IN HERE (use 'help griddata')
    phantom.fft = griddata(p_,q_,projection_fourier,p,q);
    
    % Set NaNs (Not a Number) to zero 
    phantom.fft(find(isnan(phantom.fft))) = 0;
   
    % Inverse Fourier transform
    phantom.fft = U2R(phantom.fft);

    DisplayData(phantom.fft*length(projection_angles),[2,3,6]); 
    title('2D Fourier transform-with gridding'); xlabel('x'); ylabel('y'); drawnow; 
end


%=========================================================================
function [image] = CalcDiscretePhantom(x,y,phantom,~)
    
    image = zeros(size(x));
    
    for k = 1:length(phantom.ellipse(:,1))
        
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
        image(i) = image(i)+phantom.ellipse(k,6);
        
    end
end


%=========================================================================
function [projection] = CalcLineIntegrals(r,phi,phantom,~)
    
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
function [filter] = CalcFilter(matrix,treshold)

    % --------------------------------------------------------------------
    % Define high-pass filter according to |u|
    % --------------------------------------------------------------------
    filter = abs(-fix(matrix/2):+fix(matrix/2));    % filter placeholder
    
    % TASK 2.3 FILL IN HERE
    for i = 0:matrix/2+1
        filter(filter==i) = treshold+(1-treshold)*(i/(matrix/2));
        %filter(filter==i) = 0.5 - 0.5*cos(i*pi/180);
    end
    x = 1:length(filter);
    a = 6.051e-5;
    y = 1-a*(x-0.5*length(filter)).^2;
    filter = filter.*y;
end


%=========================================================================
%=========================================================================