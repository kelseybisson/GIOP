function foq = read_fq()

%
% translated from software provided in SeaDAS/l2gen
%
% Jeremy Werdell, NASA Goddard Space Flight Center, July 2013


nw = 7; % wavelengths
ns = 6; % zenith
nc = 6; % chl
nn = 17; % nadir
na = 13; % azimuth

file = 'morel_fq.dat';
data = load(file);

foq = zeros(nw,ns,nc,nn,na);
count = 1;

for i = 1:nw
for j = 1:ns
for k = 1:nc
for l = 1:nn
	
	foq(i,j,k,l,:) = data(count,:);
	count = count + 1;

end
end
end
end


