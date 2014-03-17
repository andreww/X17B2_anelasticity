
function run_sine_fit(fileglob)


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

    figure
    subplot(2,2,1);
    plot(t(f==27&p==10), s(f==27&p==10), '.k', ...
        t(f==27&p==30), s(f==27&p==30), '.b', ...
        t(f==27&p==100), s(f==27&p==100), '.g', ...
        t(f==27&p==300), s(f==27&p==300), '.r')
    xlabel('Temperature (C)')
    ylabel('Normalised compliance')
    legend('10 s', '30 s', '100 s', '300 s');
    subplot(2,2,2);
    plot(t(f==27&p==10), phi(f==27&p==10), '.k', ...
        t(f==27&p==30), phi(f==27&p==30), '.b', ...
        t(f==27&p==100), phi(f==27&p==100), '.g', ...
        t(f==27&p==300), phi(f==27&p==300), '.r')
    xlabel('Temperature (C)')
    ylabel('Aparant internal friction (rad)')
    legend('10 s', '30 s', '100 s', '300 s');
    subplot(2,2,3);
    plot(log(p(f==27&t==25)), s(f==27&t==25), '.k', ...
        log(p(f==27&t==100)), s(f==27&t==100), '.b', ...
        log(p(f==27&t==150)), s(f==27&t==150), '.g', ...
        log(p(f==27&t==200)), s(f==27&t==200), '.r', ...
        log(p(f==27&t==250)), s(f==27&t==250), 'ok', ...
        log(p(f==27&t==300)), s(f==27&t==300), 'ob', ...
        log(p(f==27&t==400)), s(f==27&t==400), 'or')
    xlabel('log Period (s)')
    ylabel('Normalised compliance')
    legend('25 C', '100 C', '150 C', '200 C', '250 C', '300 C', '400 C');
    subplot(2,2,4);
    plot(log(p(f==27&t==25)), phi(f==27&t==25), '.k', ...
        log(p(f==27&t==100)), phi(f==27&t==100), '.b', ...
        log(p(f==27&t==150)), phi(f==27&t==150), '.g', ...
        log(p(f==27&t==200)), phi(f==27&t==200), '.r', ...
        log(p(f==27&t==250)), phi(f==27&t==250), 'ok', ...
        log(p(f==27&t==300)), phi(f==27&t==300), 'ob', ...
        log(p(f==27&t==400)), phi(f==27&t==400), 'or')
    xlabel('log period (s)')
    ylabel('Aparant internal friction (rad)')
    legend('25 C', '100 C', '150 C', '200 C', '250 C', '300 C', '400 C');


    figure
    subplot(2,2,1);
    plot(t(f==60&p==10), s(f==60&p==10), '.k', ...
        t(f==60&p==30), s(f==60&p==30), '.b', ...
        t(f==60&p==100), s(f==60&p==100), '.g', ...
        t(f==60&p==300), s(f==60&p==300), '.r')
    xlabel('Temperature (C)')
    ylabel('Normalised compliance')
    legend('10 s', '30 s', '100 s', '300 s');
    subplot(2,2,2);
    plot(t(f==60&p==10), phi(f==60&p==10), '.k', ...
        t(f==60&p==30), phi(f==60&p==30), '.b', ...
        t(f==60&p==100), phi(f==60&p==100), '.g', ...
        t(f==60&p==300), phi(f==60&p==300), '.r')
    xlabel('Temperature (C)')
    ylabel('Aparant internal friction (rad)')
    legend('10 s', '30 s', '100 s', '300 s');
    subplot(2,2,3);
    plot(log(p(f==60&t==25)), s(f==60&t==25), '.k', ...
        log(p(f==60&t==100)), s(f==60&t==100), '.b', ...
        log(p(f==60&t==200)), s(f==60&t==200), '.r', ...
        log(p(f==60&t==300)), s(f==60&t==300), 'ob', ...
        log(p(f==60&t==400)), s(f==60&t==400), 'or')
    xlabel('log Period (s)')
    ylabel('Normalised compliance')
    legend('25 C', '100 C', '200 C', '300 C',  '400 C');
    subplot(2,2,4);
    plot(log(p(f==60&t==25)), phi(f==60&t==25), '.k', ...
        log(p(f==60&t==100)), phi(f==60&t==100), '.b', ...
        log(p(f==60&t==200)), phi(f==60&t==200), '.r', ...
        log(p(f==60&t==300)), phi(f==60&t==300), 'ob', ...
        log(p(f==60&t==400)), phi(f==60&t==400), 'or')
    xlabel('log period (s)')
    ylabel('Aparant internal friction (rad)')
    legend('25 C', '100 C', '200 C', '300 C',  '400 C');

end
