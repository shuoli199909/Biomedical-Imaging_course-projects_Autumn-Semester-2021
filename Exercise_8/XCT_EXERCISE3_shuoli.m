%=========================================================================
%                                                                     
%	TITLE: 
%       XCT - EXERCISE 3
%								
%	DESCRIPTION:						
%       COMPUTE IMAGE SNR/CNR AND IMAGE ARTIFACTS
%
%	INPUT:								
%       NONE	
%
%	OUTPUT:							
%       DISPLAY
%			
%	VERSION HISTORY:						
%	    151001SK INITIAL VERSION
%       170926SK UPDATE
%       211023SK UPDATE
%
%=========================================================================

%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================
function [] = XCT_EXERCISE3()

    clear; close all; clc; 
    
    % --------------------------------------------------------------------
    % Display title
    % -------------------------------------------------------------------- 
    fprintf ( '-----------------------------------------\n' );  
    fprintf ( ' BIOMEDICAL IMAGING - XCT-EXERCISE #3\n' );  
    fprintf ( '-----------------------------------------\n' );  
     
    % --------------------------------------------------------------------
    % Set density
    % (see: Table 2 in www.nist.gov/pml/data/xraycoef/)
    % --------------------------------------------------------------------
    rho_blood       = 1.060;                % density blood     [g/cm3]
    rho_bone        = 1.450;                % density bone      [g/cm3]
    rho_lung        = 0.001;                % density lung/air  [g/cm3]
    rho_muscle      = 1.050;                % density muscle    [g/cm3]
     
    % --------------------------------------------------------------------
    % Set X-ray mass attenuation coefficients for 50 and 100 keV
    % (see: Table 4 in www.nist.gov/pml/data/xraycoef/)
    % --------------------------------------------------------------------
    mac_blood(1)    = 0.228;                % blood  @  50 keV  [cm2/g]
    mac_blood(2)    = 0.149;                % blood  @ 100 keV  [cm2/g]
    mac_bone(1)     = 0.424;                % bone   @  50 keV  [cm2/g]
    mac_bone(2)     = 0.186;                % bone   @ 100 keV  [cm2/g]
    mac_lung(1)     = 0.208;                % lung   @  50 keV  [cm2/g]
    mac_lung(2)     = 0.154;                % lung   @ 100 keV  [cm2/g]
    mac_muscle(1)   = 0.226;                % muscle @  50 keV  [cm2/g]
    mac_muscle(2)   = 0.169;                % muscle @ 100 keV  [cm2/g]
    
    % --------------------------------------------------------------------
    % Calculate linear attenuation coefficients
    % --------------------------------------------------------------------
    mue_blood       = rho_blood.*mac_blood(:);    
    mue_bone        = rho_bone.*mac_bone(:); 
    mue_lung        = rho_lung.*mac_lung(:); 
    mue_muscle      = rho_muscle.*mac_muscle(:);
    
    % --------------------------------------------------------------------
    % Set imaging parameters
    % --------------------------------------------------------------------
    Ua              = 50;                   % anode voltage       [kV]
    Ia              = 100;                  % tube current        [mA]
    Z               = 74;                   % atomic number       []
    yield           = 1e-15*Z;              % anode yield         [1/kV]
    dt              = 0.010;                % burst duration      [s]
    
    % TASK 3.1 UNCOMMENT AND FILL IN HERE
    e               = 1.60217662e-19;       % electron charge     [C]
    Ne              = Ia*0.001/e;           % #electrons          [1/s]   
    N0              = Ne*yield*Ua*dt;       % #photons            []    
    disp(['Number of incident photons: ', num2str(N0)]);
    if Ua==50  
        idx = 1; 
    end                 % set index 
    if Ua==100 
        idx = 2; 
    end
    
    matrix          = 256;                  % image matrix [1pix=1mm]
    R               = 1;                    % radial undersampling factor 
   
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
    phantom.discrete = CalcDiscretePhantom(x,y,phantom);
    figure(1);
    DisplayData(phantom.discrete,[2,3,1]); title('Phantom'); 
    
    % --------------------------------------------------------------------
    % Compute and display projections and sinogram
    % --------------------------------------------------------------------
    projection_angles = 0:R:179;        % R radial undersampling factor
                                   
    phantom.sino = zeros(matrix+1,length(projection_angles)); idx = 1;  
    
    % --------------------------------------------------------------------
    % Keep copies of half axes of heart
    % --------------------------------------------------------------------
    heart_axis = phantom.ellipse(8,4);
    
    % --------------------------------------------------------------------
    % Loop over projection angles
    % --------------------------------------------------------------------
    for phi=projection_angles
        
        % ----------------------------------------------------------------
        % Let the heart beat 
        % ----------------------------------------------------------------
        % TASK 3.4 UNCOMMENT AND FILL IN HERE
        %phantom.ellipse(8,4) = 20 + 4*sin(2*pi*dt*phi/R);
        
        phantom.discrete = CalcDiscretePhantom(x,y,phantom);
        DisplayData(phantom.discrete,[2,3,1]); title('Phantom'); 
        
        % ----------------------------------------------------------------
        % Calculate line integrals
        % ----------------------------------------------------------------
        [r,phi_] = meshgrid(-fix(matrix/2):1:+fix(matrix/2),phi);  
        phantom.sino(:,idx) = CalcLineIntegrals(r,phi_,phantom)/(matrix+1);
       
        % ----------------------------------------------------------------
        % Calculate photons incident to detector using Beer-Lambert's law 
        % ----------------------------------------------------------------
        % TASK 3.1 FILL IN HERE
        projection_mue = phantom.sino(:,idx);
        %N = repmat(N0,length(projection_mue),1);
        N = N0*exp(-projection_mue);
        % ----------------------------------------------------------------
        % Add noise
        % ----------------------------------------------------------------
        % TASK 3.1 UNCOMMENT 
        N = AddNoise(N);
        
        % ----------------------------------------------------------------
        % Calculate projection including noise
        % ----------------------------------------------------------------
        % TASK 3.1 UNCOMMENT AND FILL IN HERE
        phantom.sino(:,idx) = -log(N/N0);
        
        % ----------------------------------------------------------------
        % Display projection and sinogram
        % ----------------------------------------------------------------    
        DisplayData(phantom.sino(:,idx),[2,3,2]); 
        title(sprintf('Projection(with noise) (phi=%d)',phi)); xlabel('r'); ylabel('mue(r)'); 
        
        DisplayData(phantom.sino,[2,3,3]); 
        title('Sinogram(with noise)'); xlabel('phi'); ylabel('r'); 
 
        idx = idx+1;
    end
    
 
    % --------------------------------------------------------------------
    % Reconstruct image using Filtered Back-Projection (FBP)
    % --------------------------------------------------------------------
    [x,y]       = meshgrid(-fix(matrix/2):1:+fix(matrix/2)); idx = 1; 
    phantom.fbp = zeros(matrix+1,matrix+1);
    
    filter      = CalcHighPassFilter(matrix+1);
    
    for phi=projection_angles
        
        rs = round(x*cos(phi/180*pi)+y*sin(phi/180*pi));
        rs = rs+ceil((matrix+1)/2);
        ix = find((rs>=1)&(rs<=(matrix+1)));
        
        filteredprojection = conv(phantom.sino(:,idx),filter,'same');  
        
        phantom.fbp(ix) = phantom.fbp(ix)+filteredprojection(rs(ix));
           
        DisplayData(phantom.fbp,[2,3,4]); 
        title('Filtered Backprojection(with noise)'); xlabel('x'); ylabel('y'); drawnow;
        
        idx = idx+1;
    end
    %phantom.fbp = phantom.fbp/(4*pi^2);
    
    % --------------------------------------------------------------------
    % Reduce noise using Gaussian filter 
    % --------------------------------------------------------------------
    resolfact   = 1;                     
    phantom     = ApplyGaussianFilter(phantom,resolfact);  
    [snr]         = CalcSNR(phantom,8);           % 8 = heart (index in list)   

    DisplayData(phantom.img,[2,3,5]); xlabel('x'); ylabel('y'); title(sprintf('SNR: %3.1f',snr));
    %DisplayData(phantom.img,[2,3,5]); xlabel('x'); ylabel('y'); title(sprintf('resolution reduction factor = 1'));
    
    % --------------------------------------------------------------------
    % Reduce resolution 4-fold using Gaussian filter
    % --------------------------------------------------------------------
    resolfact   = resolfact*4;                     
    phantom     = ApplyGaussianFilter(phantom,resolfact);  
    [snr]         = CalcSNR(phantom,8);           % 8 = heart (index in list)   
    DisplayData(phantom.img,[2,3,6]); xlabel('x'); ylabel('y'); title(sprintf('SNR: %3.1f',snr));    
    %DisplayData(phantom.img,[2,3,6]); xlabel('x'); ylabel('y'); title(sprintf('resolution reduction factor = 4'));
    
end


%=========================================================================
function [projection] = AddNoise(projection)

    % --------------------------------------------------------------------
    % Create Poisson noise 
    % --------------------------------------------------------------------
    % TASK 3.1 FILL IN HERE
    projection = poissrnd(projection);
end

%=========================================================================
function phantom = ApplyGaussianFilter(phantom,factor)

    matrix      = size(phantom.fbp,1);
    [p,q]       = meshgrid(-fix(matrix/2):+fix(matrix/2)); 
    phantom.img = phantom.fbp;
    
    % --------------------------------------------------------------------
    % Define sigma based on full-width-at-half-maximum (FWHM) and the
    % resolution reduction factor 
    % --------------------------------------------------------------------
    % TASK 3.2 FILL IN HERE
    sigma = (128/factor)/sqrt(log(2));
    
    % --------------------------------------------------------------------
    % Calculate Gaussian filter
    % --------------------------------------------------------------------
    % TASK 3.2 FILL IN HERE
    filter = exp(-(p.^2+q.^2)/(sigma^2));
    
    % --------------------------------------------------------------------
    % Filter data and store into phantom.img
    % --------------------------------------------------------------------
    % TASK 3.2 UNCOMMENT HERE
    phantom.img = U2R(R2U(phantom.fbp).*filter);
    
end

%=========================================================================
function [snr_mue] = CalcSNR(phantom,k)

    %snr = 0;
    % --------------------------------------------------------------------
    % Find region-of-interest (ROI)
    % --------------------------------------------------------------------
    [x,y]   = meshgrid(-fix(size(phantom.fbp,1)/2):+fix(size(phantom.fbp,1)/2)); 
    theta   = phantom.ellipse(k,5)*pi/180;

    X0      = [x(:)'-phantom.ellipse(k,1);y(:)'-phantom.ellipse(k,2)];
    D       = [1/(phantom.ellipse(k,3)-3) 0;0 1/(phantom.ellipse(k,4)-3)];
    Q       = [cos(theta) sin(theta); -sin(theta) cos(theta)];

    equ     = sum((D*Q*X0).^2);
    i       = find(equ<=1);         % indices of phantom.img inside ROI

    
    % --------------------------------------------------------------------
    % Calculate signal-to-noise ratio (SNR)
    % --------------------------------------------------------------------    
    % TASK 3.3 FILL IN HERE
    mask    = zeros(size(phantom.img));
    mask(i) = 1;
    mask    = mask.*real(phantom.img);
    mask_vec = [];
    for i1 = 1:length(mask)
        for i2 = 1:width(mask)
            if mask(i1,i2) ~= 0
               mask_vec = [mask_vec,mask(i1,i2)]; 
            end
        end
    end
    disp(mean(mask_vec,'all'));
    snr_mue     = abs(mean(mask_vec,'all')/std(mask_vec,0,'all'));
end

%=========================================================================
%=========================================================================