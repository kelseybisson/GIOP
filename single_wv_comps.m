% function to adapt GIOP according to following changes

% 1. Compute chl based on band ratio algorithm
% 2. Assume ap(wl)  = A * Chl ^  b using chase work
% 3. Calculate bbp as D*[ C(wl) ^-eta - ap)wl))
% 4. Let cdom be ag(wl) = E*e^ -s(wl-wl0)

% U = bbp+bsw / (bbp + bsw + asw + ag + ap), and U is obtained from Rrs. 

% INPUTS - raman scattering corrected Rrs 
% Set limits on D, y, E, s. 
% ---------------------------------------------------

% LOAD RRS THAT IS RAMAN CORRECTED AND LOAD WAVELENGTH

% 1. Load float data

load('/Volumes/bissonk/floats/dailymodismatch.mat')
load('/Volumes/bissonk/floats/bbp_median.mat'); days=day; 

% 2. Identify the good points, based on where there is data

good = find(bbp1deg(:,1)>0);
 bbp =bbp(good,:); c= chl(good);
yr = yr(good);
mon = bbp(:,4); day = bbp(:,5); hr = bbp(:,6); mm = bbp(:,7);

% 3. Format dates from matlab into excel format for github code to work
% well for ocean color (nils' work)
DateNum = datenum(yr,mon,day,hr,mm, zeros(length(hr),1));
dates=DateNum-693960;
latlon= [bbp(:,1) bbp(:,2)];

bbp_443=bbp_443(good,1);
bad2 = find(isnan(bbp_443) ==1);

id = 1:14195;
latlon(bad2,:)=[]; dates(bad2)=[]; id(bad2)=[]; bbp_443(bad2)=[]; c(bad2)=[];

% 4. Use this matrix for input input into ocean color and python work!!! 
inpts = [id' dates latlon];

%  Read in matchup Rrs csv files produced through python-- MODIS

% File names are splintered bc you cant run full 4884 amount as a match
filename1 = '/Volumes/bissonk/Python_csv_matchupfiles/floats1_1000_Matchups3hr.csv';
filename2 = '/Volumes/bissonk/Python_csv_matchupfiles/floats_1001_2000Matchups3hr.csv';
filename3 = '/Volumes/bissonk/Python_csv_matchupfiles/floats_2001_3000Matchups3hr.csv';
filename4 = '/Volumes/bissonk/Python_csv_matchupfiles/floats_3000_4000Matchups3hr.csv';

one       = readtable(filename1);
two       = readtable(filename2);
three     = readtable(filename3);
four      = readtable(filename4);

% Pacakge them into one match file 
rrsmatch = [one;two;three; four];

clear one two three four

% unpack variables 
lat= table2array(rrsmatch(2:end,4));
lon = table2array(rrsmatch(2:end,5));

Rrs412 =  table2array(rrsmatch(2:end,9));
Rrs443 =  table2array(rrsmatch(2:end,10));
Rrs469 =  table2array(rrsmatch(2:end,11));
Rrs488 =  table2array(rrsmatch(2:end,12));
Rrs531 =  table2array(rrsmatch(2:end,13));
Rrs547 =  table2array(rrsmatch(2:end,14));
Rrs555 =  table2array(rrsmatch(2:end,15));
Rrs645 =  table2array(rrsmatch(2:end,16));
Rrs667 =  table2array(rrsmatch(2:end,17));
Rrs678 =  table2array(rrsmatch(2:end,18));

% Call SDs from table
Rrs412sd =  table2array(rrsmatch(2:end,29));
Rrs443sd =  table2array(rrsmatch(2:end,30));
Rrs469sd =  table2array(rrsmatch(2:end,31));
Rrs488sd =  table2array(rrsmatch(2:end,32));
Rrs531sd =  table2array(rrsmatch(2:end,33));
Rrs547sd =  table2array(rrsmatch(2:end,34));
Rrs555sd =  table2array(rrsmatch(2:end,35));
Rrs645sd =  table2array(rrsmatch(2:end,36));
Rrs667sd =  table2array(rrsmatch(2:end,37));
Rrs678sd =  table2array(rrsmatch(2:end,38));

par =  table2array(rrsmatch(2:end,25));
ipar = table2array(rrsmatch(2:end,23));

Rrs = [Rrs412 Rrs443 Rrs469 Rrs488 Rrs531 Rrs547 Rrs555 Rrs645 Rrs667 Rrs678];
Rrssd = [Rrs412sd Rrs443sd Rrs469sd Rrs488sd Rrs531sd Rrs547sd Rrs555sd Rrs645sd Rrs667sd Rrs678sd];

id2 =  table2array(rrsmatch(2:end,2));

% identify matches from the subset of float points
for i = 1: length(id2)
matches(i) = find(inpts(:,1) == id2(i));
end

bbpsofar = bbp_443(matches);
c = c(matches);

% From Rrs, correct for raman scattering of water 
a  = [0.003 0.004 0.011 0.015 0.017 0.018];
B1 = [0.014 0.015 0.010 0.010 0.010 0.010];
B2 = [ -0.022 -0.023 -0.051 -0.070 -0.080 -0.081];

% The lee et al 2013 algorithm only 
% (https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrc.20308)

% Only use WV - 412nm, 443nm, 488nm, 531nm, 551nm, 667nm
RF = a.*(Rrs443./Rrs555) + B1.*(Rrs555).^B2;

Rrs_t = [Rrs412 Rrs443 Rrs488 Rrs531 Rrs555 Rrs667];
Rrs_s = [Rrs412sd Rrs443sd Rrs488sd Rrs531sd Rrs555sd Rrs667sd];
Rrs_r = Rrs_t ./ (1 + RF); Rrs_r(Rrs_r<0)=NaN;

cd '/Users/bissonk/Documents/MATLAB/GIOP/GIOP_July2014'
% 1. Calculate chl based on band ratio 
% (from https://oceancolor.gsfc.nasa.gov/atbd/chlor_a/)

% Note: this function expects four Rrs arrays in a particular order:
% blue (~443), blue (~488/490), green (~510/531), green (~547/555/560).
% Use a placeholder for algorithms that don't use four Rrs, such as -1
% in the above examples.  Some knowledge of what bands are used in each 
% algorithm is therefore necessary; details provided here: 

chl = get_oc(Rrs_r(:,2), Rrs_r(:,3), Rrs_r(:,4), Rrs_r(:,5),'oc3m');

% 2. Assume Ap can be calculated as A*chl*b from chase 17
wl = [412,443,488,531,555,667];
ap = get_chase_ap(chl,wl);

% Get seawater contributions to abs and b (update from Zhang!) 
 asw = get_aw(wl); bsw = get_bbw(wl);
 
% take rrs across air-sea interface
 rin = Rrs_r ./ (0.52 + 1.7 .* Rrs_r);
 
rsd = Rrs_s ./ (0.52 + 1.7 .* Rrs_s);
n  = table2array(rrsmatch(2:end,54));

%% perform optimization to get eigenvalues - SIMMULATED ANNEALING
opts = optimoptions('simulannealbnd');
opts = optimoptions(opts,'FunctionTolerance',1e-9,'maxfunevals',1e6,'maxiter',1e5);

lb= [0.005,0,0.6,0.01,0];
ub= [0.01,0.4,1.4,0.02,0.5];
x0 = random('unif',lb,ub);

% if want to remove raman, just consider total rrs
 rin = Rrs_t ./ (0.52 + 1.7 .* Rrs_t);
 
% HERES WHERE YOU COULD PUT DETAILS TO ONLY TAKE LARGE # of Ns 
for i = 1: length(rin)
[x(i,:)] = simulannealbnd(@(x0) bossbisson_fmin(wl(6),x0,rin(i,6),rsd(i,6),asw(6),bsw(6),ap(i,6)),x0,lb,ub,opts)
end

% reconstruct spectra 
modb = x(:,1).*(x(:,2).*   (wl./443).^(-x(:,3))  - ap) + bsw;
moda = x(:,5).*exp(-x(:,4) .* (wl - 443))+ asw + ap;

% forward model (reconstructed) Rrs
if isempty(find(x == -999)) == 1
	modx = modb ./ (modb + moda);
    mrrs = 0.0949.* modx + 0.0794 .*modx.^2;
end

%% figures 
waves = {'412','443','488','531','555','667'};

mrrs(mrrs<0)=NaN;
figure
for i = 1:6
    subplot(3,2,i)
    scatter(rin(n>3,i),mrrs(n>3,i),20,n(n>3),'filled');   colorbar; refline(1,0)
    set(gca,'TickLabelInterpreter','latex','FontSize',16,'yscale','log','xscale','log')
    ylabel('Modeled R$_{rs}$','interpreter','latex')
    xlabel('Observed  R$_{rs}$','interpreter','latex')
    title(waves{i},'interpreter','latex')
 refline(1,0)
end

% plot how std dev changes as f(n) 
 figure
for i = 1:6
    subplot(3,2,i)
    boxplot(rsd(:,i)./mrrs(:,i),n);   
    set(gca,'TickLabelInterpreter','latex','FontSize',12,'yscale','log')
    ylabel('Standard deviation of R$_{rs}$:  R$_{rs}$ ','interpreter','latex')
    xlabel('Number of obs in 5x5','interpreter','latex')
    title(waves{i},'interpreter','latex')
    axis([1 26 0 1]);
end
 

names= {'D','C','$\gamma$','s','E'};
figure
for i = 1:5
    subplot(5,1,i)
    histogram(x(:,i))
    set(gca,'TickLabelInterpreter','latex','FontSize',16)
    ylabel('Frequency','interpreter','latex')
    xlabel('Parameter value','interpreter','latex')
    title(names{i},'interpreter','latex')
end

%% Giop normal 

% quality control; maximum Rrs_measured - Rrs_modeled spectral difference
% ranges from 0 to 1 and considers 400-600 nm
gopt.qc = 0.33; % = 33%

chl(801)=[];
rin(801,:)=[];

for i = 1:length(chl)
[x,apg,aph,adg,bbp,mrrs(i,:)] = giop(wl(6),rin(i,6),chl(i),gopt);
end
 
 