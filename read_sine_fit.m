% READ_SINE_FIT
%
% Reads the raw results of the sine fitter code from a file
% 
% Usage: [period, temperature, force, normalised_compliance, ...
%   internal_friction] = read_sine_fit(fileglob)
% 
% Andrew Walker - 20/3/2014

function [period, temperature, force, normalised_compliance, ...
           internal_friction] = read_sine_fit(file)
   
       fh = fopen(file, 'r');
       data = fscanf(fh, '%f %f %f %f %f\n', [5 inf]);
       fclose(fh);

       period = data(1,:);
       temperature = data(2,:);
       force = data(3,:);
       normalised_compliance = data(4,:);
       internal_friction = data(5,:);
end