function c = mixSpectrums(y, z, percentage)
    % Ensure percentage is within bounds
    if percentage < 0 || percentage > 100
        error('Percentage must be between 0 and 100.')
    end

    % Get the length of the arrays (they will be the same)
    len = length(y);
    % Get the middle of the arrays in regards to nyquist
    mid = len/2 + 1;  

    % Calculate how many elements we want from the first halfs of signal a
    % and b
    nA = round((1 - percentage / 100) * (mid - 2));  % -1 to exclude the Nyquist component
    nB = (mid - 2) - nA;  % -1 to exclude the Nyquist component

    % Create first half of c
    c_first_half = [y(2:2+nA-1); z(nA+2:nA+2+nB-1)];  % Take Nyquist component from b

    % Create second half of c as the flip of the first half, excluding DC and Nyquist
    c_second_half = flipud(c_first_half);  % Exclude DC and Nyquist when flipping

    % Just use dc of 1st signal, combine with first half, DC of second
    % (will be the same for both anyway) and the flipped second half
    c = [y(1:1); c_first_half; z(mid); c_second_half];


    % If selection is 100% just use b to avoid errors
     if nA == 0
         c = z;
     end


end


