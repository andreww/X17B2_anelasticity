% C_al2o3: return the elasticity matrix of alpha-al2o3 at P,T
%
% Usage:
% 
%    [C, Eiso] = C_al2o3(P, T)
%
%     C - a 6x6 matrix of the elastic cosntants of single crystal alpha-
%         al2o3 interpolated from experimental data at temperature T (in
%         Kelvin) and pressure P (in GPa).
%     E - the Young's modulus for an isotropic polycrystalline aggregate of
%         this material. Both E and C are given in GPa.
%
% The approach is to first find G at 0 GPa by interpolating the results of
% Goto et al 1989 which should be good between 296 and 1825 K. These are
% then corrected for pressure using the linear dC/dP values from Gieske and 
% Barsch 1968 (derived for pressures up to 1 GPa). This is probably the
% biggest error. Duan et al. calculated dC/dP using DFT to 150 GPa.
% Their values are similar to those of Gieske and Barsch and don't show
% much in the way of non-linear behaviour below ~50 GPa.
% 
% References:
% Goto, Anderson, Ohno and Yamamoto 1989 Elastic constants of corundum up 
%    to 1825 K. Journal of Geophysical Research 94:7588-7602.
% Gieske and Barsch 1968 Pressure dependence of the elastic constants of 
%    single crystalline aluminium oxide
% Turley and Sines 1970 The anisotropy of Young's modulus, shear modulus
%    and Poisson's ratio in cubic materials. J. Phys D. Appl. Phys 
%    4:264-271.
% Duan, Karki and Wentzcovitch 1999 High-pressure elasticity of alumina
%    studied by first principles. Am Minner 84:1961 - 1966. 


function [C, Eiso] = C_al2o3(P, T)

    % Get C at temperature T by interpolation (note T in K here)
    [C11, C33, C44, C12, C13, C14] = interpolate_T(T);
    
    % Update for pressure using pressure derivatives from table 2
    % of Gieske and Barsch 1968 assuming this contiues to be linear
    % beyond the 1GPa results they have the P/T cross terms can be 
    % neglected. See table 2 of their paper.
    C11 = C11 + P*6.174;
    C33 = C33 + P*4.998;
    C44 = C44 + P*2.243;
    C12 = C12 + P*3.282;
    C13 = C13 + P*3.653;
    C14 = C14 + P*0.130;
    
    % Put C into a matix - using high trigonal symmetry 
    % (MSAT does not do this, yet!)
    C = [C11  C12 C13  C14 0.0  0.0 ; ...
         C12  C11 C13 -C14 0.0  0.0 ; ...
         C13  C13 C33  0.0 0.0  0.0 ; ...
         C14 -C14 0.0  C44 0.0  0.0 ; ...
         0.0  0.0 0.0  0.0 C44  C14 ; ...
         0.0  0.0 0.0  0.0 C14 (C11-C12)/2.0 ];
    MS_checkC(C); % Check we have a legit tensor
    
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


function [C11, C33, C44, C12, C13, C14] = interpolate_T(T)

    % Cij data from Table 2 of Goto et al. 1989
    Ts = [...
         296   350   400   450   500   550   600   650   700 ... 
         750   800   850   900   950  1000  1050  1100  1150 ...
        1200  1250  1300  1350  1400  1450  1500  1550  1600 ...
        1650  1700  1750  1800  1825  ];
    C11s = [...
        497.3 496.5 494.7 492.7 490.6 488.2 486.0 483.8 481.5 ...
        479.3 476.8 474.5 472.3 469.9 467.4 464.8 462.5 460.0 ...
        457.3 454.5 451.9 449.2 446.7 444.3 442.2 439.5 437.2 ...
        434.8 432.3 429.6 427.2 426.0 ];
    C33s = [...
        500.9 499.0 497.2 495.4 493.6 491.3 489.2 487.2 484.9 ...
        482.6 480.4 478.2 476.0 473.8 471.2 468.9 466.4 464.0 ...
        464.1 458.7 456.2 453.1 450.8 448.3 446.4 443.7 441.3 ...
        439.0 436.5 434.5 432.5 431.4 ];
    C44s = [...
        146.8 145.7 144.4 143.1 141.8 140.5 139.2 137.9 136.5 ...
        135.2 133.9 132.6 131.2 129.9 128.6 127.2 125.8 124.5 ...
        123.2 121.7 120.4 119.0 117.7 116.3 115.1 113.8 112.5 ...
        111.3 110.0 108.7 107.4 106.7 ];
    C12s = [...
        162.8 163.7 163.8 163.8 163.7 163.3 163.1 163.1 162.9 ...
        162.8 162.4 162.3 162.4 162.1 161.8 161.5 161.4 161.1 ...
        160.7 160.3 160.0 159.5 159.5 159.2 159.4 159.0 159.0 ...
        158.9 158.4 158.1 158.0 158.1 ];
    C13s = [...
        116.0 115.6 115.3 114.8 114.4 113.7 113.0 112.5 111.9 ...
        111.2 110.6 110.1 109.6 109.1 108.2 107.6 107.1 106.4 ...
        105.4 104.8 104.1 102.9 102.4 101.6 101.6 100.8 100.5 ...
        100.0  99.4  99.2  99.1  98.9 ];
    C14s = -1.0*[...
        21.90 22.53 22.71 22.85 23.00 23.12 23.26 23.39 23.39 ...
        23.59 23.17 23.82 23.92 23.97 24.10 24.16 24.33 24.28 ...
        24.32 24.37 24.39 24.42 24.46 24.47 24.52 24.54 24.56 ...
        24.53 24.54 24.48 24.46 24.48 ];
    
    % Do the interpolation (cubic poly between points)
    C11 = interp1(Ts, C11s, T, 'pchip');
    C33 = interp1(Ts, C33s, T, 'pchip');
    C44 = interp1(Ts, C44s, T, 'pchip');
    C12 = interp1(Ts, C12s, T, 'pchip');
    C13 = interp1(Ts, C13s, T, 'pchip');
    C14 = interp1(Ts, C14s, T, 'pchip');
      
end