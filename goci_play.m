
% load GOCI products + do analysis 
clear all
load('/Volumes/data1/bisson/GOCI/files/location1b.mat');
data=B1; clear B1

load('/Volumes/data1/bisson/GOCI/files/location2b.mat');
data=B2; clear B2

load('/Volumes/data1/bisson/GOCI/files/location3b.mat');
data=B3; clear B3

data(1,:)=[];
idx = find(data(:,6)<0);
data(idx,:)=[];

lat = data(:,1);
lon = data(:,2);
year = data(:,3);
doy = data(:,4);
msec = data(:,5);
Rrs412 = data(:,6);
Rrs443 = data(:,7);
Rrs490 = data(:,8);
Rrs555 = data(:,9);
Rrs660 = data(:,10);
Rrs680 = data(:,10);

bad = find(Rrs412 < 8.05e-4 | Rrs443 < 5.49e-4 | Rrs555 < 4.48e-4 | Rrs490 < 2.51e-4 ...
    | Rrs555 < 8.83e-5);

Rrs412(bad)=[]; Rrs443(bad)=[]; Rrs490(bad)=[]; Rrs555(bad)=[];
Rrs660(bad)=[];
lat(bad)=[]; lon(bad)=[]; year(bad)=[]; doy(bad)=[];

% From Rrs, correct for raman scattering of water 
a  = [0.003 0.004 0.011  0.017 ];
B1 = [0.014 0.015 0.010 0.010 ];
B2 = [ -0.022 -0.023 -0.070 -0.080];

% The lee et al 2013 algorithm only 
% (https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrc.20308)
wl =  [412 443  490 555];

% Only use WV - 412nm, 443nm, 488nm, 531nm, 551nm, 667nm
RF = a.*(Rrs443./Rrs555) + B1.*(Rrs555).^B2;
Rrs_t = [Rrs412 Rrs443 Rrs490 Rrs555];
Rrs_r = Rrs_t ./ (1 + RF);
gopt.aw  = get_aw(wl);
gopt.bbw = get_bbw(wl);
cm = get_oc(Rrs443,Rrs490,-1,Rrs555,'oc3m');

for i = 1:length(Rrs_r)
    try
[x,mapg,maph,madg,mbbp,mrrs] = giop(wl,Rrs_r(i,:),cm(i),gopt);
    catch
        continue 
    end
bbp_giop(i,:) = mbbp;
end

bbp443 = real(bbp_giop(:,2)); bbp443(bbp443<0)=NaN;

for j = 1:length(lon)    
[zd] = timezone(lon(j));
offset(j,1)=zd;
end

load('/Volumes/data1/bisson/GOCI/files/timestamps.mat');
time_goci(1)=[];
time_goci(idx)=[];
time_goci(bad)=[];

for i = 1: length(time_goci)
junk = time_goci{i}; 
hr(i,1) = str2num(junk(12:13));
minutes(i,1) = str2num(junk(15:16));
month(i,1) = str2num(junk(6:7));
day(i,1) = str2num(junk(9:10));
end

% get local time in hours 
localh = hr+9;

cd '/Users/bissonk/Documents/MATLAB/scratch_paper'

Time.hour = hr;
Time.year = 2018;
Time.month = month;
Time.day = day;
Time.minute = minutes;
Time.second = 0;
Time.UTCOffset = zeros(length(lat),1);
Location.latitude = lat;
Location.longitude = lon;
Location.altitude = zeros(length(lat),1);

[SunAz, SunEl, ApparentSunEl]= pvl_spa(Time, Location);
% The solar elevation angle is the altitude of the Sun, the angle between 
% the horizon and the centre of the Sun's disc. Since these two angles are
% complementary, the cosine of either one of them equals the sine of the other. 

zenith = asind(cosd(SunEl));


%% now look at hours over the course of a day

[days,ca,ci] = unique([month day],'rows');

data2 = [bbp443 localh month day zenith ci];

% example plot to show diff months 
% cherry pick ref) 

figure 
subplot(2,2,1)
bbp9 = nanmean(data2(localh==9,1)); std9 = nanstd(data2(localh==9,1));
bbp10 = nanmean(data2(localh==10,1)); std10 = nanstd(data2(localh==10,1));
bbp11 = nanmean(data2(localh==11,1)); std11 = nanstd(data2(localh==11,1));
bbp12 = nanmean(data2(localh==12,1)); std12 = nanstd(data2(localh==12,1));
bbp13 = nanmean(data2(localh==13,1)); std13 = nanstd(data2(localh==13,1));
bbp14 = nanmean(data2(localh==14,1)); std14 = nanstd(data2(localh==14,1));
bbp15 = nanmean(data2(localh==15,1)); std15 = nanstd(data2(localh==15,1));
bbp16 = nanmean(data2(localh==16,1)); std16 = nanstd(data2(localh==16,1));
mseb(9:16,[bbp9 bbp10 bbp11 bbp12 bbp13 bbp14 bbp15 bbp16],...
    [std9 std10 std11 std12 std13 std14 std15 std16])

    xlabel('Local hour of day','interpreter','latex')
    ylabel('b$_{bp}$, m$^{-1}$','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','FontSize',14)
    title('Daily changes in b$_{bp}$from GOCI', 'interpreter','latex')
    
subplot(2,2,2)
scatter(bbp443,zenith,30,(data2(:,2)),'filled'); h=colorbar;
colormap(jet(8))
ylabel('Solar Zenith Angle','interpreter','latex')
ylabel(h,'Local hour of day','interpreter','latex')
xlabel('b$_{bp}$ (443nm), m$^{-1}$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex','FontSize',14)
title('Daily changes in b$_{bp}$from GOCI', 'interpreter','latex')
    
refs = [16 57 83 90 105 121 127 129 145]; 
subplot(223)
for i =1:length(refs)
    idx = find(data2(:,end) ==refs(i));
    plot(data2(idx,2),data2(idx,1),'-o','linewidth',2);
    colormap(jet(9))
    hold on
    xlabel('Local hour of day','interpreter','latex')
    ylabel('b$_{bp}$, m$^{-1}$','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','FontSize',14)
    title('Daily changes in b$_{bp}$from GOCI', 'interpreter','latex')
    
end

subplot(224)
for i =1:length(refs)
    idx = find(data2(:,end) ==refs(i));
    plot(data2(idx,2),data2(idx,5),'-o','linewidth',2);
    hold on
    xlabel('Local hour of day','interpreter','latex')
    ylabel('Solar zenith angle','interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','FontSize',14)
    title('Solar zenith angle over daily cycle, GOCI', 'interpreter','latex')
end

