% PLOT_YOUNGS_MODULUS
%
% Plot a graph of the idal Young's modulus (from the elasticity) and the 
% viscoelastic Young's modulus as measured at each period. Both as a
% function of temperature.
%
% Andrew Walker 20/3/2014.

function plot_youngs_modulus(file, varargin)

    % Defaults
    force = 27;
    
    % Process the optional arguments overiding defaults
    iarg = 1 ;
    while iarg <= (length(varargin))
        switch lower(varargin{iarg})
            case 'force'
                force = varargin{iarg+1};
                iarg = iarg + 2;
            otherwise
                error(['Unknown option: ' varargin{iarg}]) ;
        end
    end

    if force==27
        pressure = 2.5; % ish
    elseif force==60 
        pressure = 5.0; %Say
    else
        pressure = NaN;
    end
    
    % Get the data from strain fitting
    [p, t, f, s, ~, s_se, ~] = read_sine_fit(file);

    
    % Just work on 27 ton data for now.
    temperature = t(f==force);
    period = p(f==force);
    norm_compl = s(f==force);
    norm_compl_se = s_se(f==force);
    
    % Temperature limits
    t_max = max(temperature);
    t_min = min(temperature);
    Ts = t_min:10:t_max;
    n = length(Ts);

    % Get ideal Young's moduli
    Eiso_Zn = zeros(1,n);
    Emax_Zn = zeros(1,n);
    Emin_Zn = zeros(1,n);
    for i = 1:n
        [C_Zn, Eiso_Zn(i)] = C_zn(pressure, (Ts(i)+273.15));
        [Emax_Zn(i), Emin_Zn(i)] = hex_youngs_mod_lims(C_Zn);
    end
    
    % Calculate E at period and temperature from known Al2O3 E
    n = length(norm_compl);
    E_calc = zeros(1,n);
    E_calc_se = zeros(1,n);
    for i = 1:n
        [~, Eiso_Al2O3] = C_al2o3(pressure, (temperature(i)+273.15));
        E_calc(i) = (Eiso_Al2O3 / norm_compl(i));
        inv_norm_compl_se = sqrt((1/norm_compl(i))^2 * (norm_compl_se(i)/norm_compl(i)^2));
        E_calc_se(i) = sqrt(Eiso_Al2O3^2 * inv_norm_compl_se^2);
    end
    
    % Plot something useful...
    plot(Ts, Eiso_Zn, '-k', Ts, Emax_Zn, '--k', Ts, Emin_Zn, '--k');
    hold on
    errorbar(temperature(period==10), E_calc(period==10), E_calc_se(period==10), 'ok')
    errorbar(temperature(period==30), E_calc(period==30), E_calc_se(period==30), 'ob');
    errorbar(temperature(period==100), E_calc(period==100), E_calc_se(period==100), 'og');
    errorbar(temperature(period==300), E_calc(period==300), E_calc_se(period==300), 'or');
    

end

function [Emax, Emin] = hex_youngs_mod_lims(C)

    % Find the upper and lower bound of Young's modulus for a 
    % hexagonal crystal given the single crystal elastic stiffness
    % matrix C.
    
    MS_checkC(C);
    Emax = 0;
    Emin = 1000;
    
    for i = 1:0.5:90
       Cr = MS_rot3(C, 0.0, i, 0.0);
       Sr = inv(Cr);
       E = 1.0/Sr(1,1);
       if E < Emin
           Emin = E;
       end
       if E > Emax
           Emax = E;
       end
    end

end