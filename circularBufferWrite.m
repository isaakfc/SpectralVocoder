function [buffer] = circularBufferWrite(in,buffer,n)

% Determine indexes for circular buffer
len = size(buffer, 2);
indexC = mod(n-1,len) + 1; % Current index 
% Store the current output in appropriate index
buffer(:,indexC) = in;

end