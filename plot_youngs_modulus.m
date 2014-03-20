% PLOT_YOUNGS_MODULUS
%
% Plot a graph of the idal Young's modulus (from the elasticity) and the 
% viscoelastic Young's modulus as measured at each period. Both as a
% function of temperature.
%
% Andrew Walker 20/3/2014.

function plot_youngs_modulus(file)

    % Get the data from strain fitting
    [p, t, f, s, ~] = read_sine_fit(file);

    
    % Just work on 27 ton data for now.
    temperature = t(f==27);
    period = p(f==27);
    norm_compl = s(f==27);
    
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
        [C_Zn, Eiso_Zn(i)] = C_zn(2.5, (Ts(i)+273.15));
        [Emax_Zn(i), Emin_Zn(i)] = hex_youngs_mod_lims(C_Zn);
    end
    
    % Calculate E at period and temperature from known Al2O3 E
    n = length(norm_compl);
    E_calc = zeros(1,n);
    for i = 1:n
        [~, Eiso_Al2O3] = C_al2o3(2.5, (temperature(i)+273.15));
        E_calc(i) = Eiso_Al2O3 / norm_compl(i);  
    end
    
    % Plot something useful...
    plot(Ts, Eiso_Zn, '-k', Ts, Emax_Zn, '--k', Ts, Emin_Zn, '--k', ... 
        temperature(period==10), E_calc(period==10), '.k', ...
        temperature(period==30), E_calc(period==30), '.b', ...
        temperature(period==100), E_calc(period==100), '.g', ...
        temperature(period==300), E_calc(period==300), '.r');
    

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