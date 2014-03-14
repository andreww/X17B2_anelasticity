% FIT_MAXWELL_MODEL - viscoelastic model
% 
% Fits a visco-elastic maxwell model to strain data for a standard (with
% known Young's modulus) in mechanica, serise with a sample. Returns the
% estimated Young's modulus and viscosity of the sample.
% 
% Usage: [maxwell_E, maxwell_nu] = fit_maxwell_model(time, ...
%    reference_strain, reference_E, sample_strain)
% 
% maxwell_E: the estimated Young's modulus of the sample
% maxwell_nu: the estimated viscosity of the sample
% times: time stamps for the strain time serise
% reference_strain: time serise for strain in the standard
% reference_E: the Young's modulus of the (elastic) standard
% sample_strain: time serise for the (visco-elastic) sample
%
% Options:
%     'minimise_options' - optimset output to pass to fitter
% 
% Andrew Walker 14/3/2014

function [maxwell_E, maxwell_nu] = fit_maxwell_model(time, ...
    reference_strain, reference_E, sample_strain, varargin)

    % Default options
    minimise_options = optimset('Display', 'off');
    % process optionsal arguments
    iarg = 1 ;
    while iarg <= (length(varargin))
        switch lower(varargin{iarg})
            case 'minimise_options'
                minimise_options = varargin{iarg + 1};
                iarg = iarg + 2;
            otherwise 
               error(['Unknown option in fit_maxwell_model: ' ...
                   varargin{iarg}]);
        end
    end
    

    % Calculate stress in reference material
    stresses = reference_strain * reference_E;
    
    % Fit Maxwell model
    maxwell_coef = fminsearch(@maxwell_model_misfit, ...
        [reference_E/2.0, 0.0000001], minimise_options, time, ...
         sample_strain, stresses);
   
    maxwell_E = maxwell_coef(1);
    maxwell_nu = maxwell_coef(2);
   
    % Print the results
    fprintf(['Maxwell model: \n', ...
        '    E sample = %6.3f \n' ...
        '    viscosity = %6.3f \n'], ...
        maxwell_E, maxwell_nu);
    
    % Graph the results
    fitted_maxwell_displacments = maxwell_model(time, stresses, ...
        maxwell_E, maxwell_nu);
    figure
    subplot(1,1,1)
    plot(time, sample_strain, '.g', time, fitted_maxwell_displacments, '-g');
    legend('Measured strain', 'Viscoelastic strain');
    xlabel('Timestamp (s)')
    ylabel('Displacment (px)')
    
end
    
function [sum_sq] = maxwell_model_misfit(coeff, times, data, stress)

    E = coeff(1);
    p = coeff(2);
    model_data = maxwell_model(times, stress, E, p);
    sum_sq = sum((model_data - data).^2);
    
end

function [elastic_displacment, plastic_displacement] = maxwell_model(times, ...
    stress, elasticity, viscosity)

    elastic_displacment = stress / elasticity;
    plastic_rate = stress / viscosity;
    plastic_displacement = zeros(size(times));
    delta_times = zeros(size(times));
    for i = 2:length(times)
        delta_times(i) = times(i) - times(i-1);
        plastic_displacement(i) = plastic_displacement(i-1) + ...
            plastic_rate(i-1)*times(i);
    end

end