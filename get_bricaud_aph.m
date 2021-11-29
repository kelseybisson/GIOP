function [ap,aps,aph,aphs] = get_bricaud_aph(chl,wl,norm)

%
% Bricaud et al., JGR 103, 31033-31044, 1998
%
% ap = Ap Chl^Ep
% ap* = Ap Chl^(Ep-1)
% aph = Aph Chl^Eph
% aph* = Aph Chl^(Eph-1)
% bricaud_1998_aph.txt order of columns: wl,Ap,Ep,Aph,Eph
%
% Jeremy Werdell, NASA Goddard Space Flight Center, July 2013
%

dat = load('bricaud_1998_aph.txt');

ap0   = dat(:,2) .* chl.^dat(:,3);
aps0  = dat(:,2) .* chl.^(dat(:,3) - 1);
aph0  = dat(:,4) .* chl.^dat(:,5);
aphs0 = dat(:,4) .* chl.^(dat(:,5) - 1);


if exist('norm') == 1
	idx = find(dat(:,1) == 442);
	aphs0 = aphs0 * 0.055 / aphs0(idx);
end


ap   = interp1(dat(:,1),  ap0,wl,'pchip');
aps  = interp1(dat(:,1), aps0,wl,'pchip');
aph  = interp1(dat(:,1), aph0,wl,'pchip');
aphs = interp1(dat(:,1),aphs0,wl,'pchip');


