% PLOT_SINE_FIT
%
% Plots the raw results of the sine fitter code.
% 
% Usage: [period, temperature, force, normalised_compliance, ...
%   internal_friction] = run_sine_fit(fileglob)
% 
%    NOTE - A HACK: assumes we have exactly the right data.
%
% Andrew Walker - 20/3/2014

function plot_sine_fit(file)

    [p, t, f, s, phi, s_se, phi_se] = read_sine_fit(file);

    plot_force = unique(f);
    plot_periods = [10 30 100 300];
    period_point = {'ok' 'ob' 'og' 'or'};

    for i = 1:length(plot_force)
        figure
        subplot(2,1,1);
        for j = 1:length(plot_periods)
            plot(t(f==plot_force(i)&p==plot_periods(j)), ...
                s(f==plot_force(i)&p==plot_periods(j)), ...
                period_point{j})
            errorbar(t(f==plot_force(i)&p==plot_periods(j)), ...
                s(f==plot_force(i)&p==plot_periods(j)), ...
                s_se(f==plot_force(i)&p==plot_periods(j)))
            hold on;
        end
        xlabel('Temperature (C)')
        ylabel('Normalised compliance')
        legend('10 s', '30 s', '100 s', '300 s');
        hold off;
        subplot(2,1,2);
        for j = 1:length(plot_periods)
            plot(t(f==plot_force(i)&p==plot_periods(j)), ...
                phi(f==plot_force(i)&p==plot_periods(j)), ...
                period_point{j})
            errorbar(t(f==plot_force(i)&p==plot_periods(j)), ...
                phi(f==plot_force(i)&p==plot_periods(j)), ...
                phi_se(f==plot_force(i)&p==plot_periods(j)))
            hold on
        end
        xlabel('Temperature (C)')
        ylabel('Aparant internal friction (rad)')
        legend('10 s', '30 s', '100 s', '300 s');
        hold off;
    end
    
end
