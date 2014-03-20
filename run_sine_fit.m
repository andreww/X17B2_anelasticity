% RUN_SINE_FIT
%
% Runs the sine fitter code over a list of foil position files
% and stores the results.
% 
% Usage: [period, temperature, force, normalised_compliance, ...
%   internal_friction] = run_sine_fit(fileglob)
% 
% Andrew Walker - 20/3/2014

function [p, t, f, s, phi] = run_sine_fit(fileglob)


    files = dir(fileglob);

    p = zeros(size(files));
    t = zeros(size(files));
    f = zeros(size(files));
    s = zeros(size(files));
    phi = zeros(size(files));

    fh = fopen('SineFits.txt', 'w');

    for i=1:size(files)
        [p(i), t(i), f(i), s(i), phi(i)] = sine_fit_strain(files(i).name);
        %pause
        [~,fn,~] = fileparts(files(i).name);
        fn1 = [fn '_sine_fit.png'];
        fn2 = [fn '_data.png'];
        figure(1)
        print('-dpng', fn1)
        figure(2)
        print('-dpng', fn2)
        close all
        fprintf(fh, '%9.7f %9.7f %9.7f %9.7f %9.7f\n', p(i), t(i), f(i), s(i), phi(i));
    end

    fclose(fh);

end
