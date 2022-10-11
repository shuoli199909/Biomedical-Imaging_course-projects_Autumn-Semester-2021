%=========================================================================
%                                                                     
%       BIOMEDICAL IMAGING
%       MRI 3
%
%=========================================================================

%=========================================================================
%	FMRI
%=========================================================================

function [] = FMRI_shuoli()

    clear; close all; clc;
    
    % Load data
    [nx,ny,nt,fmri_data,anatomical_data,paradigm] = load_data();                % Load data
                                                                                % anatomical image of one brain slice 
                                                                                % series of 200 functional images of the same slice
                                                                                % paradigm = time course of visual stimulus
            
    % Show anatomical data
    figure(1);
    subplot(1,2,1);
    imshow(mat2gray(anatomical_data));
    title('Anatomical data');
    
    % Show paradigm
    subplot(1,2,2);
    plot(paradigm),ylim([-2 2]);
    title('Paradigm');
    
    % Show functional data  
    figure(2);
    for i = 1:10
        t = 20*(i-1)+1;
        subplot(2,5,i), imshow(mat2gray(fmri_data(:,:,t))); 
        title(strcat('Frame ',num2str(t)));
    end
    
    
    % TASK: Try to find pixels with temporal fluctuation that resembles the paradigm
    
    figure(3);
    % Use a for loop to search appropriate x and y
    for i1 = 1:1:length(anatomical_data)
        for i2 = 1:1:width(anatomical_data)
            % normalization
            data_temp = squeeze(fmri_data(i1,i2,:));
            data_temp_sort = sort(data_temp);
            data_temp = data_temp - mean([data_temp_sort(end-1),data_temp_sort(2)]); 
            data_temp_sort = sort(data_temp);
            data_temp = data_temp/data_temp_sort(end-1);
            % caculate difference
            if i1+i2 == 2
                diff_min = sum(abs(data_temp-paradigm'));
            else
                diff_temp = sum(abs(data_temp-paradigm'));
                if diff_temp < diff_min
                    diff_min = diff_temp;
                    y = i1;
                    x = i2;
                else
                    continue;
                end
            end
        end
    end
    plot(squeeze(fmri_data(y,x,:)));
    title('Pixel with temporal fluctuation that resembles the paradigm');
    
    % TASK: Calculate an activation map
    
    product = zeros(ny,nx,nt);
    for i = 1:1:nt
        product(:,:,t) = fmri_data(:,:,t).*paradigm(i);
    end
    scalar_product = sum(product,3);
    figure(4);
    imagesc(scalar_product);
    title('Activation map');
        
    % TASK: Calculate map of temporal standard deviation, determine noise level
    
    std_data = std(fmri_data,0,3);
    figure(5);
    imshow(mat2gray(log(std_data)));
    title('Map of temporal standard deviation');
    % We should calculate the standard deviation without activation
    range_bg = 120;
    noise_level = sqrt(mean(mean(std_data(1:1:range_bg,1:1:range_bg).^2)));
         
    % TASK: Display thresholded activation map superimposed on anatomical data
    
    factor = 4.5;
    mask = (abs(scalar_product) > factor*sqrt(nt)*noise_level);
    figure(6);
    imshow(mat2gray(2*anatomical_data.*(~mask)+mask.*scalar_product));
    title('Masked activation superimposed on the anatomical data');
             
end



function [nx,ny,nt,fmri_data,anatomical_data,paradigm] = load_data()

    nx = 382;
    ny = 482;
    nt = 200;
    
    fmri_data = zeros(ny,nx,nt);
    anatomical_data = zeros(ny,nx);
    paradigm = zeros(1,nt);
    
    fileID = fopen('fMRI_data.bin');
    for t = 1:nt
        fmri_data(:,:,t) = fread(fileID,[ny,nx],'float');
    end
    fclose(fileID);
    
    fileID = fopen('anatomical_data.bin');
    anatomical_data(:,:) = fread(fileID,[ny,nx],'float');
    fclose(fileID);
    
    fileID = fopen('paradigm.bin');
    paradigm(:) = fread(fileID,[nt],'float');
    fclose(fileID);
    
end
    
