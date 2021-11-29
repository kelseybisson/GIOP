function [val] = morel_read(chl,solz,type)

%
% read f, f', and mu_d LUTs provided by Andre Morel
%
% Jeremy Werdell, NASA Goddard Space Flight Center, July 2013
%

if strcmp(type,'mud')
	file = 'morel_mud.txt';
elseif strcmp(type,'f')
	file = 'morel_f.txt';
elseif strcmp(type,'fp')
	file = 'morel_fp.txt';
end

% file columns are solz, wl, chl0.03, chl0.1, chl0.3, chl1, chl3, chl10
%
table = load(file);

sz = table(:,1);

c = [0.03,0.1,0.3,1,3,10];
w = unique(table(:,2));
[cx,wx] = meshgrid(c,w);

wl = 380:1:700;

% isolate values for each solz
%
dat0 = table(find(sz ==  0),3:8);
dat1 = table(find(sz == 15),3:8);
dat2 = table(find(sz == 30),3:8);
dat3 = table(find(sz == 45),3:8);
dat4 = table(find(sz == 60),3:8);
dat5 = table(find(sz == 75),3:8);

% for each solz, interpolate to input chl and wl
%
d00 = interp2(cx,wx,dat0,chl,wl);
d15 = interp2(cx,wx,dat1,chl,wl);
d30 = interp2(cx,wx,dat2,chl,wl);
d45 = interp2(cx,wx,dat3,chl,wl);
d60 = interp2(cx,wx,dat4,chl,wl);
d75 = interp2(cx,wx,dat5,chl,wl);

% interpolate to input solz
%	  
flo = d00;
fhi = d75;
tlo = 0;
thi = 75;

if solz < 15
	fhi = d15;
	thi = 15;
elseif solz >= 15 & solz < 30
	flo = d15;
	fhi = d30;
	tlo = 15;
	thi = 30;
elseif solz >= 30 & solz < 45
	flo = d30;
	fhi = d45;
	tlo = 30;
	thi = 45;
elseif solz >= 45 & solz < 60
	flo = d45;
	fhi = d60;
	tlo = 45;
	thi = 60;
elseif solz >= 60
	flo = d60;
	tlo = 60;
end

[wx,sx] = meshgrid(wl,[tlo,thi]);
val = interp2(wx,sx,[flo';fhi'],wl,solz);

% replace end NaNs
%
inan = find(isnan(val));
val(inan) = val(inan(1)-1);
