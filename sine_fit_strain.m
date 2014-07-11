% SINE_FIT_STRAIN
%
% Analysis of X17B2 anelasticity experiments. Given a 
% 'position change' file (from Simon's code) containing 
% the length of a standard and a sample, calculate 
% the strains on both and attempt to extract the 
% phase lag and strain amplitude ratio (the internal 
% friction and normalised compliance, respectivly). The 
% approach is to first detrend the data and then fit the 
% phase, period and amplitude in sequence for the sample
% and standard. A global fit (with common period) is then 
% performed before the data is plotted and returned.
% Options to vary the way the data is read and processed
% are avalable.
%
% Arguments: 
%     filename: the filename (or path and filename) or the 
%               'positions change file.
%
% Options:
%     optim_options: set options for the optimiser. Default is 
%                    optimset('Display', 'final', 'TolX', 1e-9, ...
%                    'Tolfun', 1e-9)
%     trim_data_start: remove a fraction of the data from the start of the 
%                time serise. e.g. (..., 'trim_data_start', 0.1) strips out 
%                the first 10% of the data serise.
%     trim_data_end: remove a fraction of the data from the end of the 
%                time serise. e.g. (..., 'trim_data_end', 0.1) strips out 
%                the last 10% of the data serise.
%
%     mode: set a method for measuring strain. Options are 'ES-', '-SE',
%                'ESE', 'E-S', 'S-E', 'SE', 'ES'. Number of characters 
%                must be one larger than the number of foils in the file.
%                E is elastic standard, S is the sample, '-' ignores one
%                element. If two E's are present the lengths are added 
%                before calculating the strain.
% 
% See also: run_sine_fit, plot_sine_fit

% Andrew Walker <a.walker@leeds.ac.uk> - 20/5/2014

function [nom_period, temperature, load,...
    normalised_compliance, internal_friction, normalised_compliance_se, ...
    internal_friction_se] ...
    = sine_fit_strain(filename, varargin)

    % Set defaults for options. 
    minimise_options = optimset('Display', 'final', 'TolX', 1e-9, ...
        'Tolfun', 1e-9);
    trim_data_start = 0.0;
    trim_data_end = 0.0;
    strain_mode = 'default';
    
    ScreenSize = get(0,'ScreenSize');
    fig_width = ScreenSize(3)/2;
    fig_height = ScreenSize(4)/2;
    
    % Process the optional arguments overiding defaults
    iarg = 1 ;
    while iarg <= (length(varargin))
        switch lower(varargin{iarg})
            case 'optim_opts'
                minimise_options = varargin{iarg+1};
                iarg = iarg + 2;
            case 'trim_data_start'
                trim_data_start = varargin{iarg+1};
                iarg = iarg + 2;
            case 'trim_data_end'
                trim_data_end = varargin{iarg+1};
                iarg = iarg + 2;
            case 'mode'
                strain_mode = varargin{iarg+1};
                iarg = iarg + 2;
            otherwise
                error(['Unknown option: ' varargin{iarg}]) ;
        end
    end
    
    
    % Read in the data and process it so we can make use.
    [foil_change_data, time, box_positions, number_boxes, ...
       number_images] = ReadPositionChangeFile(filename);
    [metadata] = read_position_change_metadata(filename);
   
    [~, name, ~] = fileparts(filename);
    name = strrep(name, '_', ' ');
   
    % Time should start at 0
    time_zero = time(1);
    time = time - time_zero;
    

    % Get the nominal period for fitting with
    nom_period = metadata.NominalPeriod;
    assert(nom_period > 0, 'Error, nomimal period must be positive');
              
    temperature = metadata.NominalTemp;
    load = metadata.NominalLoad;
 
    if strcmpi(strain_mode, 'default')
        if number_boxes == 3
            strain_mode = 'SE';
        elseif number_boxes == 4
            strain_mode = 'ESE';
        else
            error('Can only work with 3 or 4 foils!');
        end  
    end
        
    
    switch lower(strain_mode)
        case 'es-'
            std = [1 2];
            sam = [2 3];
        case '-se'
            std = [3 4];
            sam = [2 3];
        case 'ese'
            std = [1 2 3 4];
            sam = [2 3];
        case 'e-s'
            std = [1 2];
            sam = [3 4];
        case 's-e'
            std = [3 4];
            sam = [1 2];
        case 'se'
            std = [2 3];
            sam = [1 2];
        case 'es'
            std = [1 2];
            sam = [2 3];
        otherwise
            error('Unkown mode argument')
    end
    
    if length(std) == 4
        
        error('Not implemented')
        
    elseif length(std) == 2
    
    % pull out individual boxes...
    %box1 = foil_change_data(:,1);
    %box2 = foil_change_data(:,2);
    %box3 = foil_change_data(:,3);
    
    % Switch the data around to agree with the rest of the script.
    time = time';
    %box1 = box1';
    %box2 = box2';
    %box3 = box3';
    
    % Calculate strain of top and bottom block, at each time step.
    % For intresting historical reasons, bot_ is always the elastic 
    % standard and top_ is always the sample.
    top_ref_length = (box_positions(sam(1),1)+box_positions(sam(1),2))/2.0 + ...
        (box_positions(sam(2),1)+box_positions(sam(2),2))/2.0 ;
    bot_ref_length = (box_positions(std(1),1)+box_positions(std(1),2))/2.0 + ...
        (box_positions(std(2),1)+box_positions(std(2),2))/2.0 ;
    
    top_strain = (foil_change_data(:,sam(2)) - foil_change_data(:,sam(1))) / top_ref_length;
    top_strain = top_strain';
    bot_strain = (foil_change_data(:,std(2)) - foil_change_data(:,std(1))) / bot_ref_length;
    bot_strain = bot_strain';
    end
    
    % Optionally throw out data from the start, e.g. where the first
    % cycle is bad. We can plot this up (distinctivly) but we do not
    % use it for fitting.
    if (trim_data_start + trim_data_end ~= 0.0)
        s = 1 + floor(numel(time).*trim_data_start);
        e = numel(time) - floor(numel(time)*trim_data_end);
        unused_time = [time(1:s-1) time(e+1:numel(time))];
        unused_top_strain = [top_strain(1:s-1) top_strain(e+1:numel(time))];
        unused_bot_strain = [bot_strain(1:s-1) bot_strain(e+1:numel(time))];
        time = time(s:e);
        top_strain = top_strain(s:e);
        bot_strain = bot_strain(s:e);
    else
        unused_time = [];
        unused_top_strain = [];
        unused_bot_strain = [];
    end
    
    % Detrend the data and calculate detrended strains
    [detrend_top, bg_top, ~] = detrend_data(time, top_strain, 'verbose',...
        'minimise_options', minimise_options, 'name', 'top data');
    [detrend_bot, bg_bot, ~] = detrend_data(time, bot_strain, 'verbose',...
        'minimise_options', minimise_options, 'name', 'bottom data');
    
    % Guess amplitudes
    top_nom_amp = (sqrt(mean(detrend_top.^2))*sqrt(2));
    bot_nom_amp = (sqrt(mean(detrend_bot.^2))*sqrt(2));
    
    % Guess the phase shift
    top_nom_phase = nominal_phase(time, detrend_top, nom_period, ...
        top_nom_amp);
    bot_nom_phase = nominal_phase(time, detrend_bot, nom_period, ...
        bot_nom_amp);
    
    % Fit sine coeffs for two models allowing different periods...
    [period_top, amplitude_top, phase_top, sine_top] = ...
        fit_single_serise(nom_period, top_nom_amp, top_nom_phase,...
        minimise_options, time, detrend_top, 'top');

    [period_bot, amplitude_bot, phase_bot, sine_bot] = ...
        fit_single_serise(nom_period, bot_nom_amp, bot_nom_phase,...
        minimise_options, time, detrend_bot, 'bottom');
    
    
    % Now use theset to do a joint fit
    sine_coef_both = fminsearch(@model_misfit_both, ...
        [((period_top+period_bot)/2.0) amplitude_top amplitude_bot ...
         phase_top phase_bot], minimise_options, time, ...
         [detrend_top; detrend_bot]);
    
    period = sine_coef_both(1);
    amplitude_top = sine_coef_both(2);
    amplitude_bot = sine_coef_both(3);
    phase_top = sine_coef_both(4);
    phase_bot = sine_coef_both(5);
    fprintf(['Joint fit: period = %6.4f \n' ...
        '    amplitude top = %6.4e phase top = %6.4f\n'...
        '    amplitude bot = %6.4e phase bot = %6.4f\n\n'], ...
        period, amplitude_top, phase_top, amplitude_bot, phase_bot);
  
    period_error = 0;
    if ((abs(nom_period - period)/period) > 0.02)
        fprintf (['Nominal period, %6.4f, and fitted' ...
                  ' period, %6.4f, differ by > 2 percent'], nom_period, period);
        period_error = 1;
    end
    
    
    % Now estimate the errors on the parameters - see 
    % http://www.mathworks.co.uk/matlabcentral/newsreader/view_thread/157530
    % and links therein. Note that this bit needs the Adaptive Robust
    % Derivative toolbox to estimate the Jacobian. See:
    % http://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation#
    % Also note that we do this twice to get errors on the standard and sample, seperatly.
    % Should probably do the calculation on the joint model - but how many 
    % degrees of freedom for this?
    dof = numel(time) - 2; %  Degrees of freedom
    % Sine func
    sine_predict = @(coeff) sine_model(time, coeff(1), coeff(2), coeff(3));
    % standard deviation of the residuals - top
    sdr = sqrt(sum((detrend_top - sine_model(time, period, ...
        amplitude_top, phase_top)).^2)/dof);
    % jacobian matrix
    J = jacobianest(sine_predict,[period, ...
        amplitude_top, phase_top]);
    % Get the errors
    Sigma = sdr^2*inv(J'*J);
    se = sqrt(diag(Sigma))';
    se_period_top = se(1);
    se_amplitude_top = se(2);
    se_phase_top = se(3);
    % standard deviation of the residuals - bottome
    sdr = sqrt(sum((detrend_bot - sine_model(time, period, ...
        amplitude_bot, phase_bot)).^2)/dof);
    % jacobian matrix
    J = jacobianest(sine_predict,[period, ...
        amplitude_bot, phase_bot]);
    % Get the errors
    Sigma = sdr^2*inv(J'*J);
    se = sqrt(diag(Sigma))';
    se_period_bot = se(1);
    se_amplitude_bot = se(2);
    se_phase_bot = se(3);
    % Bit of a hack
    se_period = (se_period_top + se_period_bot) / 2.0;
    fprintf(['Estimated standard errors on parameters: \nperiod = %6.4f \n' ...
        '    amplitude top = %6.4e phase top = %6.4f\n'...
        '    amplitude bot = %6.4e phase bot = %6.4f\n'], ...
        se_period, se_amplitude_top, se_phase_top, se_amplitude_bot, se_phase_bot);
    
    % Draw a graph
    
    fitted_data_top = sine_model(time, period, amplitude_top, phase_top);
    fitted_data_bot = sine_model(time, period, amplitude_bot, phase_bot);
    resid_top = (detrend_top - fitted_data_top)./amplitude_top;
    resid_bot = (detrend_bot - fitted_data_bot)./amplitude_bot;
    
    fig1 = figure;
    set(fig1, 'OuterPosition', [1, fig_height, fig_width, fig_height]);
    subplot(2,1,1)
    plot(time, detrend_top, 'or', time, detrend_bot, 'ob', ...
        time, sine_top, ':r', time, sine_bot, ':b', ...
        time, fitted_data_top, '--r', time, fitted_data_bot, '--b');
    xlim([min(time) max(time)])
    legend('Detrended Zn sample', 'Detrended Al_2O_3 standard');
    xlabel('Timestamp (s)')
    ylabel('Strain')
    title(name);
   
    subplot(2,1,2)
    plot(time, resid_top, '+r', time, resid_bot, '+b', ...
        time, zeros(size(time)), '-k');
    xlim([min(time) max(time)])
    xlabel('Timestamp (s)')
    ylabel('Normalised residual (fraction of modelled amplitude')
    
    fig2 = figure;
    set(fig2, 'OuterPosition', [fig_width, fig_height, fig_width, fig_height]);
    subplot(2,1,1)
    plot(time, top_strain, 'or', time, bg_top, '-r', ...
        unused_time, unused_top_strain, '+r');
    legend('Zn sample data', 'Background');
    xlabel('Timestamp (s)')
    ylabel('Strain')
    
    subplot(2,1,2)
    plot(time, bot_strain, 'ob', time, bg_bot, '-b',...
        unused_time, unused_bot_strain, '+b');
    legend('Al2O3 sample data', 'Background');
    xlabel('Timestamp (s)')
    ylabel('Strain')
    
    if 0
    E_al2o3 = 350; % Youngs mod of elastic standard (Al2O3)
    [~, ~] = fit_maxwell_model(time, ...
    fitted_data_bot, E_al2o3, detrend_top, 'minimise_options',...
    minimise_options);
    end    

    normalised_compliance = amplitude_top / amplitude_bot;
    normalised_compliance_se = (sqrt((se_amplitude_top/amplitude_top)^2 + ...
        (se_amplitude_bot/amplitude_bot)^2)) * normalised_compliance;
    
    phase_bot_frac = (phase_bot/period);
    phase_bot_rad = phase_bot_frac*(2.0*pi);
    phase_bot_frac_se = sqrt ( phase_bot_frac^2 * ( ...
        (se_phase_bot/phase_bot)^2+(se_period/period)^2 ) );
    phase_bot_rad_se = sqrt( phase_bot_frac_se^2 * (2.0*pi)^2);
  
    phase_top_frac = (phase_top/period);
    phase_top_rad = phase_top_frac*(2.0*pi);
    phase_top_frac_se = sqrt ( phase_top_frac^2 * ( ...
        (se_phase_top/phase_top)^2+(se_period/period)^2 ) );
    phase_top_rad_se = sqrt( phase_top_frac_se^2 * (2.0*pi)^2);
    
    internal_friction = abs(phase_bot_rad - phase_top_rad);
    internal_friction_se = sqrt(phase_bot_rad_se^2 + phase_top_rad_se^2);
    
    if (internal_friction < 0)
        % Something went wrong - elastic should not lag viscoelastic
        fprintf(['Calculated phase offset, %6.4f' ...
            ' is negative'], internal_friction);
        internal_friction = NaN;
        normalised_compliance = NaN;
    elseif (internal_friction > pi/2.0)
        % Something went wrong - maximum lag (pure viscosity) is 90 degrees
        fprintf(['Calculated phase offset, %6.4f' ...
            ' is more than 90 degrees'], internal_friction);
        internal_friction = NaN;
        normalised_compliance = NaN;
    elseif (period_error)
        internal_friction = NaN;
        normalised_compliance = NaN;
    end
    
    fprintf(['\nNormalised compliance = %6.4f (s.e. = %6.4f) \n' ...
        'Internal friction = %6.4f (s.e. = %6.4f) rad\n'] , ...
        normalised_compliance, normalised_compliance_se, ...
        internal_friction, internal_friction_se);
    
end

function [best_phase] = nominal_phase(time, data, period, amplitude)

    best_phase = 0;
    best_sum_sq = 10000;
    num_phase = 20;
    for i = 1:num_phase
        phase_guess = (i/num_phase)*period;
        model = sine_model (time, period, amplitude, phase_guess);
        sum_sq = sum((model - data).^2);
        if (sum_sq < best_sum_sq)
            best_phase = phase_guess;
            best_sum_sq = sum_sq;
        end
    end

end

function [sum_sq] = model_misfit_both(coeff, times, data)

    period = coeff(1);
    amplitude_1 = coeff(2);
    amplitude_2 = coeff(3);
    phase_1 = coeff(4);
    phase_2 = coeff(5);
    data_1 = data(1,:);
    data_2 = data(2,:);
    
    model_1 = sine_model(times, period, amplitude_1, phase_1) - data_1;
    model_2 = sine_model(times, period, amplitude_2, phase_2) - data_2;
    
    sum_sq = sum(([model_1 model_2]).^2); 

end

function [sum_sq] = model_misfit(coeff, times, data)

    period = coeff(1);
    amplitude = coeff(2);
    phase = coeff(3);
    a = coeff(4);
    b = coeff(5);
    c = coeff(6);
    
    calc_data = background_sine_model(times, period, amplitude, phase, ...
        a, b, c);
    
    residuals = calc_data - data;
    sum_sq = sum(residuals.^2);

end

function [sum_sq] = sine_misfit(coeff, times, data)

    period = coeff(1);
    amplitude = coeff(2);
    phase = coeff(3);
    
    sum_sq = sum((sine_model(times, period, amplitude, phase) - data).^2);

end


function [background_sine_data] = background_sine_model(times, period, ...
                                    amplitude, phase, a, b, c)
                                
    background_sine_data = sine_model(times, period, amplitude, ...
        phase) + background_model(times, a, b, c);
    
end


function [period, amplitude, phase, fitted_data] = fit_single_serise(...
             nom_period, nom_amp, nom_phase, minimise_options, time, ...
             detrend_data, string)

    sine_coef_top = fminsearch(@sine_misfit, [nom_period nom_amp nom_phase], ...
        minimise_options, time, detrend_data);
    fitted_data = sine_model(time, sine_coef_top(1), sine_coef_top(2),...
        sine_coef_top(3));
    period = sine_coef_top(1);
    amplitude = sine_coef_top(2);
    phase = sine_coef_top(3);
    fprintf(['Sine fit for %s data: period = %6.4f amplitude = %6.4f' ...
        ' phase = %6.4f\n'], string, period, amplitude, phase);
end


function [ sine_data ] = sine_model (times, period, amplitude, phase)

    % phase and period in seconds, amp in px, times in seconds
    sine_data = sin(((times + phase)/period)*(2.0*pi)) * amplitude;
    
end
