% C_zn: return the elasticity matrix of zinc at P,T
%
% Usage:
% 
%    [C, Eiso] = C_zn(P, T)
%
%     C - a 6x6 matrix of the elastic cosntants of single crystal zinc
%         interpolated from experimental data at temperature T (in
%         Kelvin) and pressure P (in GPa).
%     E - the Young's modulus for an isotropic polycrystalline aggregate of
%         this material. Both E and C are given in GPa.
%
% The approach is to first find G at 0 GPa by interpolating the results of
% Alers and Neighbours 1958 which should be good between 4 and 670 K. 
% These are then corrected for pressure using the linear dC/dP values 
% from Srinivasan and Ramiji Rao 1971. The elasticiciy of zinc is reviewed
% by Ledbetter 1977.
% 
% References:
% Alers and Neighbours 1958 The elastic constants of zinc between 4.2^\circ
%    and 670^\circK. J. Phys. Chem. Solids 7:58-64
% Ledbetter 1977 Elastic properties of zinc: a compilation and a review
%    J. Phys. Chem. Ref. Data 6:1181-1203
% Turley and Sines 1970 The anisotropy of Young's modulus, shear modulus
%    and Poisson's ratio in cubic materials. J. Phys D. Appl. Phys 
%    4:264-271.
% Srinivasan and Ramiji Rao 1971 Anharmonic properties of the hexagoanl 
%    metals, magnesium zinc and beryllium. I. attice dynamcis and third
%    order elastic constants. J. Phys. Chem. Solids 32:1769-1788.


function [C, Eiso] = C_zn(P, T)

    % Get C at temperature T by interpolation (note T in K here)
    [C11, C33, C44, C13, C66] = interpolate_T(T);
    
    % Update for pressure using pressure derivatives from 
    % of Srinivasan et al (1971). See Table 10 of Ledbetter (1977)
    C11 = C11 + P*7.70;
    C33 = C33 + P*6.50;
    C44 = C44 + P*3.94;
    C13 = C13 + P*4.72;
    C66 = C66 + P*2.10;
    
    % Put C into a matix
    C = [C11 0.0 C13 0.0 0.0 0.0 ; ...
         0.0 0.0 0.0 0.0 0.0 0.0 ; ...
         0.0 0.0 C33 0.0 0.0 0.0 ; ...
         0.0 0.0 0.0 C44 0.0 0.0 ; ...
         0.0 0.0 0.0 0.0 0.0 0.0 ; ...
         0.0 0.0 0.0 0.0 0.0 C66 ];
    C = MS_expand(C, 'hex'); 
    
    % Get the isotropic avarage Young's modulus. This is 
    % Done by first calculating the elastic constants 
    % matrix for an anisotropic poly crystal, then the young's
    % modulus, E, is just 1/S_11 (S is the compliance matrix).
    % See Turley and Sines 1971, for example.
    [ K_vrh, G_vrh ] = MS_polyaverage( C );
    [Ciso]=MS_build_isotropic('K', K_vrh, 'mu', G_vrh);
    Siso = inv(Ciso);
    Eiso = 1.0/Siso(1,1); % Or 2,2 or 3,3

end


function [C11, C33, C44, C13, C66] = interpolate_T(T)

    % Cij data from Table 3 of Alers and Neighbours 2958. Note that C66 is 
    % marked C11-C12/2 and we do not need C12.
    Ts = [...
         4.2    77     150    200    250    295    350    400 ...
         450   500     550    600    670];
    C11s = [...
         179.09 176.77 173.36 170.47 167.30 163.38 160.02 155.90 ...
         151.38 146.48 142.19 138.43 133.95];
    C33s = [...
         68.80  67.66  66.20  25.23  64.28  63.47  62.57  61.68  ...
         60.77  59.82  58.87  57.93  56.61];
    C44s = [...
         45.95  44.79  42.96  41.58  40.18  38.79  37.24  35.73  ...
         34.17  32.61  30.98  29.33  26.67];
    C66s = [...
         70.80  70.00  68.37  66.89  65.33  63.61  61.94  60.06  ...
         57.92  55.67  53.40  50.95  46.53];
    C13s = [...
         54.4   55.2   54.8   54.5   53.7   53.0   52.1   51.5   ...
         51.1   50.8   50.6   50.4   50.3];
    
    
    % Do the interpolation (cubic poly between points)
    C11 = interp1(Ts, C11s, T, 'pchip');
    C33 = interp1(Ts, C33s, T, 'pchip');
    C44 = interp1(Ts, C44s, T, 'pchip');
    C13 = interp1(Ts, C13s, T, 'pchip');
    C66 = interp1(Ts, C66s, T, 'pchip');
      
end