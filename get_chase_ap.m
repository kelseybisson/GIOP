function [ap] = get_chase_ap(chl,wl)

% Chase et al, 2017...Estimation of phytoplankton accessory pigments from
% hyperspectral reflectance spectra: Toward a global algorithm
%
% ap = Ap Chl^Ep

% Kelsey Bisson, Oregon State University, 2019

file = load('chase_ap17.mat'); dat= file.chase_ap17;
ap0   = dat(:,2).* chl'.^dat(:,3);

idx = find(isnan(chl)==1);

for i =1:length(ap0)
    
    if i == idx
        ap(i,:) = NaN;
        continue 
    end
    ap(i,:)   = interp1(dat(:,1),  ap0(:,i),wl,'pchip');
end

end


