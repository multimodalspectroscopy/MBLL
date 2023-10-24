function conc = MBLL(spectra, wavelengths, optode_dist, DPF)

% Calculate concentrations from raw intensity spectra using UCLn Modified Beer Lambert Law 

% Gemma Bale 2020

% spectra = raw intensity spectra
% wavelengths = of intensity spectra in nm
% optode_dist = separation between optode in cm
% DPF = differential pathlength factor: 4.99 baby head, 6.26 adult head, 4.16 adult arm, 5.51 adult leg (from Duncan 1994)
% conc = change in concentration (mM) of HbO2, HHb and oxCCO

%% calculation

wl = 780:900; % Wavelengths used, nm

% Extinction coefficients
% Wavelength, HbO2, HHb, CCO, wavelength dependency
[specific_extinction_coeffs] = specific_extinction_coeffs_780to900;

% HbO2, HHb, oxCCO
ext_coeffs = specific_extinction_coeffs(:,2:4);
ext_coeffs_inv = pinv(ext_coeffs);

[x,y] = size(spectra);

% Preallocate to improve performance
atten = zeros([x,y]);
atten_int = zeros([length(wl),x]);

% change in attenuation
for i = 1:x 
        atten(i,:) = log10(spectra(1,:)./spectra(i,:)); % log 10
        atten_int(:,i) = interp1(wavelengths, atten(i,:)', wl, 'spline'); % interpolate to wavelengths of interest
end

% wavelength dependency of pathlength
atten_int_wldep = atten_int./ specific_extinction_coeffs(:,5); 

% MBLL
conc =  (ext_coeffs_inv * atten_int_wldep .* 1/(optode_dist*DPF))';

end
