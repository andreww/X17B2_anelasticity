

function sine_fit_strain(filename)

    % Options - FIXME: should be optional input arguments. 
    minimise_options = optimset('Display', 'final', 'TolX', 1e-9, 'Tolfun', 1e-9);

    % first read the data file
    [~, ~, time, box1, box2, box3, metadata] ...
        = read_position_change(filename);
    
    time_zero = time(1);
    time = time - time_zero;
    
    [~, name, ~] = fileparts(filename);
    name = strrep(name, '_', ' ');
    
    % Lengths of top and bottom samples in first image 
    % FIXME: ask Simon why this is calculated this way.
    top_ref_length = (metadata.BoxTops(1)+metadata.BoxBotoms(1))/2.0 + ...
                     (metadata.BoxTops(2)+metadata.BoxBotoms(2))/2.0;               
    bot_ref_length = (metadata.BoxTops(2)+metadata.BoxBotoms(2))/2.0 + ...
                     (metadata.BoxTops(3)+metadata.BoxBotoms(3))/2.0;
                 
    % Get the nominal period for fitting with
    nom_period = metadata.NominalPeriod;
              
    % Calculate strain of top and bottom block, at each time step.             
    top_strain = (box2 - box1) / top_ref_length;
    bot_strain = (box3 - box2) / bot_ref_length;

    % Detrend the data and calculate detrended strains
    % FIXME - this should be two calls to a fnction
    fprintf('Detrending the data with a quadratic fit...\n');
    bg_coef_top = fminsearch(@background_misfit, [0.0 0.0 0.0], ...
        minimise_options, time, top_strain);
    a_top = bg_coef_top(1);
    b_top = bg_coef_top(2);
    c_top = bg_coef_top(3);
    fprintf('Best fit background for top strains a = %6.2g b = %6.2g c = %6.2g\n',...
        a_top, b_top, c_top);
    bg_top = background_model(time, a_top, b_top, c_top);
    detrend_top = top_strain-bg_top;
    
    bg_coef_bot = fminsearch(@background_misfit, [0.0 0.0 0.0], ...
       minimise_options, time, bot_strain);
    a_bot = bg_coef_bot(1);
    b_bot = bg_coef_bot(2);
    c_bot = bg_coef_bot(3);
    fprintf('Best fit background for bottom strains a = %6.2g b = %6.2g c = %6.2g\n',...
        a_bot, b_bot, c_bot);
    bg_bot = background_model(time, a_bot, b_bot, c_bot);
    detrend_bot = bot_strain-bg_bot;
    
    % Guess amplitudes
    top_nom_amp = (sqrt(mean(detrend_top.^2))*sqrt(2));
    bot_nom_amp = (sqrt(mean(detrend_bot.^2))*sqrt(2));
    
    % Guess the phase shift
    top_nom_phase = nominal_phase(time, detrend_top, nom_period, ...
        top_nom_amp);
    bot_nom_phase = nominal_phase(time, detrend_bot, nom_period, ...
        bot_nom_amp);
    
    % Fit sine coeffs for two models allowing different periods...
    % FIXME - again this should be a function
    sine_coef_top = fminsearch(@sine_misfit, [nom_period top_nom_amp top_nom_phase], ...
        minimise_options, time, detrend_top);
    sine_top = sine_model(time, sine_coef_top(1), sine_coef_top(2),...
        sine_coef_top(3));
    period_top = sine_coef_top(1);
    amplitude_top = sine_coef_top(2);
    phase_top = sine_coef_top(3);
    fprintf(['Sine fit for top data: period = %6.2g amplitude = %6.2g' ...
        ' phase = %6.2g\n'], period_top, amplitude_top, phase_top);
    
    
    sine_coef_bot = fminsearch(@sine_misfit, [nom_period bot_nom_amp bot_nom_phase], ...
        minimise_options, time, detrend_bot);
    sine_bot = sine_model(time, sine_coef_bot(1), sine_coef_bot(2),...
       sine_coef_bot(3));
    period_bot = sine_coef_bot(1);
    amplitude_bot = sine_coef_bot(2);
    phase_bot = sine_coef_bot(3);
    fprintf(['Sine fit for bottom data: period = %6.2g amplitude = %6.2g' ...
        ' phase = %6.2g\n'], period_bot, amplitude_bot, phase_bot);
    
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
    fprintf(['Joint fit: period = %6.2g \n' ...
        '    amplitude top = %6.2g phase top = %6.2g\n'...
        '    amplitude bot = %6.2g phase bot = %6.2g\n'], ...
        period, amplitude_top, phase_top, amplitude_bot, phase_bot);
    
    % Draw a graph
    
    fitted_data_top = sine_model(time, period, amplitude_top, phase_top);
    fitted_data_bot = sine_model(time, period, amplitude_bot, phase_bot);
    resid_top = (detrend_top - fitted_data_top)./amplitude_top;
    resid_bot = (detrend_bot - fitted_data_bot)./amplitude_bot;
    
    figure
    subplot(2,1,1)
    plot(time, detrend_top, 'or', time, detrend_bot, 'ob', ...
        time, sine_top, ':r', time, sine_bot, ':b', ...
        time, fitted_data_top, '--r', time, fitted_data_bot, '--b');
    legend('Detrended Zn sample', 'Detrended Al_2O_3 standard');
    xlabel('Timestamp (s)')
    ylabel('Displacment (px)')
    title(name);
   
    subplot(2,1,2)
    plot(time, resid_top, '.r', time, resid_bot, '.b', ...
        time, zeros(size(time)), '-k');
    xlabel('Timestamp (s)')
    ylabel('Normalised residual')
    
    figure
    subplot(2,1,1)
    plot(time, top_strain, '.r', time, bg_top, '-r');
    legend('Zn sample data', 'Background');
    xlabel('Timestamp (s)')
    ylabel('Displacment (px)')
    
    subplot(2,1,2)
    plot(time, bot_strain, '.b', time, bg_bot, '-b');
    legend('Al2O3 sample data', 'Background');
    xlabel('Timestamp (s)')
    ylabel('Displacment (px)')
    
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

function [sum_sq] = background_misfit(coeff, times, data)

    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    
    sum_sq = sum((background_model(times, a, b, c) - data).^2);

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

function [background_data] = background_model(times, a, b, c)

    % times in seconds, a, b and c in px, px/sec and px/sec^2
    background_data = a.*times.^2 + b.*times + c;

end

function [ sine_data ] = sine_model (times, period, amplitude, phase)

    % phase and period in degrees, amp in px, times in seconds
    sine_data = sin(((times + phase)/period)*(2.0*pi)) * amplitude;
    
end