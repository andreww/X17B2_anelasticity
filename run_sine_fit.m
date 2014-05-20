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
        status = check_fit_ok;
        while status == 2
            prompt = {'Fraction to remove at start:',...
                'Fraction to remove at end:'};
            dlg_title = 'Cleanup ends';
            num_lines = 1;
            def = {'0.0','0.0'};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            strip_start = str2double(answer{1});
            strip_end = str2double(answer{2});
            close all
            [p(i), t(i), f(i), s(i), phi(i)] = sine_fit_strain(...
                files(i).name, 'trim_data_start', strip_start,...
                'trim_data_end', strip_end);
            status = check_fit_ok;
        end
        if status == 1
            [~,fn,~] = fileparts(files(i).name);
            fn1 = [fn '_sine_fit.png'];
            fn2 = [fn '_data.png'];
            figure(1)
            print('-dpng', fn1)
            figure(2)
            print('-dpng', fn2)
            fprintf(fh, '%9.7f %9.7f %9.7f %9.7f %9.7f\n', p(i), t(i), f(i), s(i), phi(i));
        end
        close all
    end

    fclose(fh);

end


function ok = check_fit_ok

    % Construct a questdlg
    choice = questdlg('Is the fit OK?', 'Fit checker', ...
	    'Yes','Filter','No','No');
    % Handle response
    switch choice
        case 'Yes'
            ok = 1;
        case 'Filter'
            ok = 2;
        case 'No'
            ok = 0;
    end
end