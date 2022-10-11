%=========================================================================
%                                                                     
%	TITLE: 
%       NUC - EXERCISE 2
%								
%	DESCRIPTION:						
%       QUANTITATIVE PET DATA ANALYSIS
%
%	INPUT:								
%       NONE	
%
%	OUTPUT:							
%       DISPLAY
%			
%	VERSION HISTORY:						
%	    201114SK INITIAL VERSION
%       211023SK UPDATE
%
%=========================================================================

%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================
function [] = NUC_EXERCISE2()

    clear; close all; clc;
    
    % --------------------------------------------------------------------
    % Display title
    % -------------------------------------------------------------------- 
    fprintf ( '-----------------------------------------\n' );  
    fprintf ( ' BIOMEDICAL IMAGING - NUC-EXERCISE #2\n' );  
    fprintf ( '-----------------------------------------\n' );  
     
    % --------------------------------------------------------------------
    % Set rate constants of kinetic model
    % --------------------------------------------------------------------
    k(1) = 0.1;                                         %  [1/min]
    k(2) = 0.3;                                         %  [1/min]
    k(3) = 0.5;                                         %  [1/min]
    
    % --------------------------------------------------------------------
    % Convert from [1/min] to [1/sec]
    % --------------------------------------------------------------------
    k(:) = k(:)/60;                                     % [1/min -> 1/sec]
  
    % --------------------------------------------------------------------
    % Task 2.1 Set up first-order differential equations to obtain impulse
    % response function h(t) (see help dsolve)
    % --------------------------------------------------------------------
    syms k1 k2 k3 ce(t) cm(t)
    
     h = dsolve(diff(ce) == -k2*ce(t)-k3*ce(t),...
                diff(cm) == k3*ce(t),...
                   ce(0) == k1, ...
                   cm(0) == 0 );
  
    % --------------------------------------------------------------------
    % Extract impulse response functions h_ce(t) and h_cm(t) from h(t)
    % --------------------------------------------------------------------
    h_ce  = matlabFunction(h.ce);
    h_cm  = matlabFunction(h.cm);
    k1  = k(1);
    k2  = k(2);
    k3  = k(3);
    
    % --------------------------------------------------------------------
    % Set time span for experiment 
    % --------------------------------------------------------------------
    t   = 0:240;                                        % [sec]
    
    % --------------------------------------------------------------------
    % Set up input blood plasma function cp(t) using gamma variate pdf
    % --------------------------------------------------------------------
    cp  = gampdf(t,2,10);                               % [ml/g] 
    cp  = cp/max(cp(:));      
    
    % --------------------------------------------------------------------
    % Task 2.2. Convolve cp(t) with impulse response h(t) to obtain 
    % extracellular concentration ce(t) and metabolized concentration cm(t)
    % --------------------------------------------------------------------
     ce = filter(h_ce(k1,k2,k3,t),1,cp);
     cm = filter(h_cm(k1,k2,k3,t),1,cp);
     ct = ce+cm;
    
    % --------------------------------------------------------------------
    % Display cp(t), ce(t) and cm(t)
    % --------------------------------------------------------------------
    figure(1); subplot(2,1,1); plot(t,cp,'red'); title('plasma cp(t)'); 
    figure(1); subplot(2,1,2); plot(t,ce,'green'); hold on;
    figure(1); subplot(2,1,2); plot(t,cm,'blue'); hold on;
    figure(1); subplot(2,1,2); plot(t,ct,'black'); 
    legend('extracellular ce(t)','metabolized cm(t)','tissue ct(t)');
    
    % --------------------------------------------------------------------
    % Add noise to cp(t) and ct(t) according to given SNR and fit 
    % concentration-time curves to estimate k1, k2, k3 (repeat 10x)
    % --------------------------------------------------------------------
    SNR = 100;                                          % SNR input 
    
    for i = 1:10
    
        % ----------------------------------------------------------------
        % Task 2.3. Convert concentrations into counts N by assuming
        % a peak SNR of blood plasma activity of 100
        % ----------------------------------------------------------------
         N0      = SNR^2;
         Ncp     = cp*N0;
         Nct     = ct*N0;
        
        % ----------------------------------------------------------------
        % Task 2.3. Add noise to cp(t) and ct(t) = ce(t) + cm(t)
        % ----------------------------------------------------------------
         Ncp_noise = poissrnd(Ncp);
         Nct_noise = poissrnd(Nct);
         
        % ----------------------------------------------------------------
        % Task 2.3.Convert from counts N back to concentrations
        % ----------------------------------------------------------------
         cp_noise  = Ncp_noise/N0;
         ct_noise  = Nct_noise/N0;
        
        % ----------------------------------------------------------------
        % Display noisy concentration-time curves
        % ----------------------------------------------------------------
        figure(2); subplot(2,1,1); plot(t,cp_noise,'red');   title('plasma cp(t)'); hold on;
        figure(2); subplot(2,1,2); plot(t,ct_noise,'black'); title('tissue ct(t)'); hold on;
     
        % ----------------------------------------------------------------
        % Fit noisy data input cp_noise, ct_noise using Matlab's nlinfit
        % to obtain k1, k2, k3
        % ----------------------------------------------------------------
        
        % ----------------------------------------------------------------
        % Task 2.2 Define fit function using h_t(t) = h_ce(t) + h_cm(t)
        % ct_noise(t) = int(cp_noise(tau)*(h_t(t-tau))dtau  
        % use filter() to perform the convolution
        % ----------------------------------------------------------------
        fun = @(k,t)filter(h_ce(k(1),k(2),k(3),t)+h_cm(k(1),k(2),k(3),t),1,cp_noise);
        
        k0 = [0.01 0.01 0.01];                          % Starting values for fit
        
        k_fit  = nlinfit(t,ct_noise,fun,k0);
        
        % ----------------------------------------------------------------
        % Display fitted concentration ct_fit(t) and estimated rates k1 k2 k3
        % ----------------------------------------------------------------
        figure(2); subplot(2,1,2); 
        plot(t,fun(k_fit,t),'yellow*-'); hold on;
        k_result(:,i) = k_fit(:)*60;                    % [1/sec -> 1/min]  
        
        fprintf('Estimated rate constants [k1 k2 k3]: %f %f %f\n',k_result);
        
    end
    
    % --------------------------------------------------------------------
    % Task 2.3. Display mean +/- standard deviation of k1, k2, k3
    % --------------------------------------------------------------------
     fprintf('Mean  [k1 k2 k3]: %f %f %f\n',mean(k_result,2));
     fprintf('StDev [k1 k2 k3]: %f %f %f\n',std(k_result,0,2));
          
end

  
%=========================================================================
%=========================================================================