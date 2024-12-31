function c = mixSpectrums(y, z, percentage)
    % Ensure percentage is within bounds
    if percentage < 0 || percentage > 100
        error('Percentage must be between 0 and 100.')
    end

    % Get the length of the arrays (they will be the same)
    len = length(y);
    % Get the middle of the arrays in regards to nyquist
    mid = len/2 + 1;  

    % Calculate amount of elements we want for a and b
    A = round((1 - percentage / 100) * (mid - 2));  
    B = (mid - 2) - A;  

    % Create first half of c
    c_first_half = [y(2:2+A-1); z(A+2:A+2+B-1)]; 

    % Create second half of c as the flip of the first half
    c_second_half = flipud(c_first_half);  

    % Just use dc of 1st signal, combine with first half, DC of second
    % (will be the same for both anyway) and the flipped second half
    c = [y(1:1); c_first_half; z(mid); c_second_half];


    % If selection is 100% just use b to avoid errors
     if A == 0
         c = z;
     end


end


