function [oc] = get_oc(r1,r2,r3,r4,type)

%
% function to call empirical OC Chl and Kd algorithms
%
% both algorithms have the form:
% model = 10^(a0 + a1*x + a2*x^2 + a3*x^3 + a4*x^4)
% where a0 ... a4 are polynomial regression coefficients,
% x = log10(MBR), and MBR is maximum Rrs band ratio  
%
% usage:
% chl = get_oc(Rrs443,Rrs490,Rrs510,Rrs555,type)
% where Rrs### is an array of Rrs at wavelength ###
% and type is 'OC4' or 'KD2S' or ... (see below)
%
% examples:
% for SeaWiFS Chl: oc4 = get_oc(Rrs443,Rrs490,Rrs510,Rrs555,'oc4')
% for MODIS Chl: oc3m = get_oc(Rrs443,Rrs488,-1,Rrs547,'oc3m')
% for SeaWiFS Kd: kd2s = get_oc(-1,Rrs490,-1,Rrs555,'kd2s')
%
% Note: this function expects four Rrs arrays in a particular order:
% blue (~443), blue (~488/490), green (~510/531), green (~547/555/560).
% Use a placeholder for algorithms that don't use four Rrs, such as -1
% in the above examples.  Some knowledge of what bands are used in each 
% algorithm is therefore necessary; details provided here: 
%
% http://oceancolor.gsfc.nasa.gov/ANALYSIS/ocv6/
% http://oceancolor.gsfc.nasa.gov/ANALYSIS/kdv4/
%
% Jeremy Werdell, NASA Goddard Space Flight Center, 31 Jan 2012
%

% coefficients for the v6 (operational as of 2009) OC and KD algs

c = [0.3272,-2.9940,2.7218,-1.2259,-0.5683;...  % OC4        1 SeaWiFS operational Chl
     0.3255,-2.7677,2.4409,-1.1288,-0.4990;...  % OC4E       2 MERIS operational Chl
	 0.3325,-2.8278,3.0939,-2.0917,-0.0257;...	% OC4O 		 3 OCTS operational Chl
	 0.2515,-2.3798,1.5823,-0.6372,-0.5692;...	% OC3S		 4 SeaWiFS 3 band Chl
	 0.2424,-2.7423,1.8017, 0.0015,-1.2280;...	% OC3M		 5 MODIS operational Chl
	 0.2521,-2.2146,1.5193,-0.7702,-0.4291;...	% OC3E		 6 MERIS 3 band Chl
	 0.2399,-2.0825,1.6126,-1.0848,-0.2083;...	% OC3O		 7 OCTS 3 band Chl
	 0.3330,-4.3770,7.6267,-7.1457, 1.6673;...	% OC3C		 8 CZCS operational Chl
	 0.2511,-2.0853,1.5035,-3.1747, 0.3383;...	% OC2S		 9 SeaWiFS 2 band Chl
	 0.2389,-1.9369,1.7627,-3.0777,-0.1054;...	% OC2E		10 MERIS 2 band Chl
	 0.2236,-1.8296,1.9094,-2.9481,-0.1718;...	% OC2O		11 OCTS 2 band Chl
	 0.2500,-2.4752,1.4061,-2.8237, 0.5405;...	% OC2M		12 MODIS 2 band Chl
	 0.1464,-1.7953,0.9718,-0.8319,-0.8073;...	% OC2M-HI	13 MODIS high-res band Chl
    -0.8515,-1.8263,1.8714,-2.4414,-1.0690;...	% KD2S		14 SeaWiFS operational Kd
	-0.8813,-2.0584,2.5878,-3.4885,-1.5061;...	% KD2M		15 MODIS operational Kd
	-0.8641,-1.6549,2.0112,-2.5174,-1.1035;...	% KD2E		16 MERIS operational Kd
	-0.8878,-1.5135,2.1459,-2.4943,-1.1043;...	% KD2O		17 OCTS operational Kd
	-1.1358,-2.1146,1.6474,-1.1428,-0.6190;...  % KD2C		18 CZCS operational Kd
	 0.2228,-2.4683,1.5867,-0.4275,-0.7768;...  % OC3V		19 VIIRS operational Chl
	 0.2230,-2.1807,1.4434,-3.1709, 0.5863;...	% OC2V		20 VIIRS 2 band Chl
	-0.8730,-1.8912,1.8021,-2.3865,-1.0453];    % KD2V		21 VIIRS operational Kd



% select coefficients and generate maximum band ratio

switch lower(type)
    case 'oc4e'
        a = c(2,:);
        r = log10(max([r1,r2,r3]')' ./ r4);
    case 'oc4o'
        a = c(3,:);
        r = log10(max([r1,r2,r3]')' ./ r4);
    case 'oc3s'
        a = c(4,:);
        r = log10(max([r1,r2]')' ./ r4);
    case 'oc3m'
        a = c(5,:);
        r = log10(max([r1,r2]')' ./ r4);
    case 'oc3e'
        a = c(6,:);
        r = log10(max([r1,r2]')' ./ r4);
    case 'oc3o'
        a = c(7,:);
        r = log10(max([r1,r2]')' ./ r4);
    case 'oc3c'
        a = c(8,:);
        r = log10(max([r1,r2]')' ./ r4);
    case 'oc2s'
        a = c(9,:);
        r = log10(r2 ./ r4);
    case 'oc2e'
        a = c(10,:);
        r = log10(r2 ./ r4);
    case 'oc2o'
        a = c(11,:);
        r = log10(r2 ./ r4);
    case 'oc2m'
        a = c(12,:);
        r = log10(r2 ./ r4);
    case 'oc2m-hi'
        a = c(13,:);
        r = log10(r2 ./ r4);
    case 'kd2s'
        a = c(14,:);
        r = log10(r2 ./ r4);
    case 'kd2m'
        a = c(15,:);
        r = log10(r2 ./ r4);
    case 'kd2e'
        a = c(16,:);
        r = log10(r2 ./ r4);
    case 'kd2o'
        a = c(17,:);
        r = log10(r2 ./ r4);
    case 'kd2c'
        a = c(18,:);
        r = log10(r2 ./ r3);
    case 'oc3v'
        a = c(19,:);
        r = log10(max([r1,r2]')' ./ r4);
    case 'oc2v'
        a = c(20,:);
        r = log10(r2 ./ r4);
    case 'kd2v'
        a = c(21,:);
        r = log10(r2 ./ r4);
    otherwise 
        a = c(1,:);
        r = log10(max([r1,r2,r3]')' ./ r4);
end

% calculate modeled parameter

oc = 10.^(a(1) + a(2)*r + a(3)*r.^2 + a(4)*r.^3 + a(5)*r.^4);



