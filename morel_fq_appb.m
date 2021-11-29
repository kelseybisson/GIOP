function [f,q] = morel_fq_appb(chl,solz)

%
% Appendix B from Morel et al. 2002, Appl. Opt.
%
% Jeremy Werdell, NASA Goddard Space Flight Center, July 2013
%

dat = load('morel_fq_appb.txt');

f0 = dat( 1: 6,:);
sf = dat( 7:12,:);
q0 = dat(13:18,:);
sq = dat(19:24,:);

wlut = [412.5, 442.5, 490, 510, 560, 620, 660];
clut = [0.03, 0.1, 0.3, 1, 3, 10];
[wx,cx] = meshgrid(wlut,clut);

wl = 380:1:700;
ilo = find(wl < wlut(1));
wl(ilo) = wlut(1);
ihi = find(wl > wlut(7));
wl(ihi) = wlut(7);

z = 1. - cos(solz * pi / 180);

f1 = f0 + sf * z;
q1 = q0 + sq * z;

f = interp2(wx,cx,f1,wl,chl);
q = interp2(wx,cx,q1,wl,chl);

f = fastsmooth(f,25,1,1);
q = fastsmooth(q,25,1,1);

