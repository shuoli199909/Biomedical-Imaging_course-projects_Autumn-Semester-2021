function [filter] = CalcFilter_(matrix)

    % --------------------------------------------------------------------
    % Define high-pass filter according to |u|
    % --------------------------------------------------------------------
    filter = abs(-fix(matrix/2):+fix(matrix/2));    % filter placeholder
      
    % TASK 2.3 FILL IN HERE
    a = 1; b = 0; c = matrix/(1*2.355);
    threshold = 64;
    for u = 0:fix(matrix/2)
        % filter(filter==u) = 1 - a*exp(-(u-b)^2/(2*c^2));
        if u<threshold
            filter(filter==u) = 0; 
        else
            filter(filter==u) = 1;
        end
    end
        
end