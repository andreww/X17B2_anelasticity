% DETREND DATA
%
% Fit a background model to displacment data and remove the long-time
% drift.
% 
% Usage: [detrended_data, trend, params] = detrend_data(time, data, ...)
%
%     time - time for each data point
%     data - the displacments 
%     detrended_data - data with the trend removed
%     trend   - the best fitting long term trend
%     params   - parameters of the long-term trend
% 
% Options:
%     'verbose' - print output to display
%     'function' - bacground function to fit, default is 'quadratic',
%           other functions not implemented.
%     'minimise_options' - optimset output to pass to fitter
%     'name' - a name for the data (to print to display).
% 
% Andrew Walker - 11/3/2014

function [detrended_data, trend, params] = detrend_data(time, data, ...
    varargin)

    % Default options
    verbose = 0;
    method = 'quadratic';
    minimise_options = optimset('Display', 'off');
    name = 'data';
    % process optionsal arguments
    iarg = 1 ;
    while iarg <= (length(varargin))
        switch lower(varargin{iarg})
            case 'verbose'
               verbose = 1;
               iarg = iarg + 1;
            case 'function'
                method = varargin{iarg + 1};
                iarg = iarg + 2;
            case 'minimise_options'
                minimise_options = varargin{iarg + 1};
                iarg = iarg + 2;
            case 'name'
                name = varargin{iarg + 1};
                iarg = iarg + 2;
            otherwise 
               error(['Unknown option in detrend_data: ' varargin{iarg}]);
        end
    end
    
    if verbose
        fprintf('Detrending the %s with a %s fit...\n', name, ...
            lower(method));
    end
    
    switch lower(method)
        case 'quadratic'
            params = fminsearch(@quadratic_misfit, [0.0 0.0 0.0], ...
                minimise_options, time, data);
            % FIXME - check fit worked...
            if verbose
                fprintf(['Best fit to function y = ax^2+bx+c is: \n'...
                    '    a = %6.2g b = %6.2g c = %6.2g\n\n'],...
                    params(1), params(2), params(3));
            end
            trend = quadratic_model(time, params(1), ...
                params(2), params(3));
        otherwise 
            error(['Unknown function in detrend_data: ' varargin{iarg}]);
    end
    detrended_data = data-trend;
    
end

function [sum_sq] = quadratic_misfit(coeff, times, data)

    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    
    sum_sq = sum((quadratic_model(times, a, b, c) - data).^2);

end

function [background_data] = quadratic_model(times, a, b, c)

    % times in seconds, a, b and c in px, px/sec and px/sec^2
    background_data = a.*times.^2 + b.*times + c;

end