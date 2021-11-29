function [x,mapg,maph,madg,mbbp,mrrs] = giop(wl,rrs,chl,gopt)

% GIOP ocean color reflectance inversion model
%
% P.J. Werdell and 18 co-authors, "Generalized ocean color 
% inversion model for retrieving marine inherent optical 
% properties," Appl. Opt. 52, 2019-2037 (2013).
%
% inputs are vectors of wavelength and Rrs, an estimate of
% chlorophyll, plus a structure (gopt) to control the 
% parameterization of the inversion
%
% outputs are a vector of the magnitudes of the eigenvalues
% for adg, bbp, and aph (x), plus modeled spectra of apg,
% aph, adg, bbp, and Rrs
%
% see run_giop.m for details
%
% Jeremy Werdell, NASA Goddard Space Flight Center, July 2013
%

mapg = -999;
maph = -999;
madg = -999;
mbbp = -999;
mrrs = -999;

% require 412, 443, and 547/555 to be present
%
i412 = find(wl >= 410 & wl <= 415);
i443 = find(wl >= 441 & wl <= 445);
i555 = find(wl >= 546 & wl <= 557);
if isempty(i412) | isempty(i443) | isempty(i555)
	'cannot find Rrs(412), Rrs(443), or Rrs(547,555)'
	return
end


% default variables
%
if exist('gopt') == 0; gopt = struct; end
if exist('chl') == 0; chl = 0.2; end
if isfield(gopt,'solz') == 0; gopt.solz = 30; end

% coefficients for pure seawater
%
if isfield(gopt,'aw') == 0; gopt.aw = get_aw(wl); end
if isfield(gopt,'bbw') == 0; gopt.bbw = get_bbw(wl); end


% AOP to IOP relationship
%
% Gordon default
if isfield(gopt,'fq') == 0; gopt.fq = [0.0949,0.0794]; end

% Morel 2002
if strcmp(gopt.fq,'morel')
	if isfield(gopt,'senz') == 0 & isfield(gopt,'relaz') == 0
		[f,q] = morel_fq_appb(chl,gopt.solz);
		f = morel_read(chl,gopt.solz,'fp');
		fq = interp1(380:1:700,f./q,wl,'cubic');
	else
		[fq,fc] = morel_fq(chl,gopt.solz,gopt.senz,gopt.relaz);
		fq = interp1(380:1:700,fq,wl,'cubic');
    end
    gopt.g0 = fq;
    gopt.g1 = 0;
% Gordon (if specifically named)
elseif strcmp(gopt.fq,'gordon')
    gopt.g0 = 0.0949;
    gopt.g1 = 0.0794;
% user defined
else
    gopt.g0 = gopt.fq(1);
    gopt.g1 = gopt.fq(2);
end


% Bricaud 1998 spectra with aph*(443) normalized to 0.055 m2/mg
%
[bap,baps,baph,baphs] = get_bricaud_aph(chl,wl,'norm');


% take Rrs across the air-sea interface (default to Lee) 
%
rin = rrs ./ (0.52 + 1.7 .* rrs);
if isfield(gopt,'trans') == 1
	if strcmp(gopt.trans,'flat')
		% spectrally flat transmission a la Gordon 2007
		rin = rrs ./ 0.529;
	end
end


% define eta (default to QAA)
%
eta = 2 * (1 - 1.2 * exp(-0.9 * rin(i443(1))/rin(i555(1))));
if isfield(gopt,'eta') == 1
	if strcmp(gopt.eta,'qaa'); 
		% from QAA version 5 (Lee et al. 2002, etc.)
		eta = eta;	
	elseif strcmp(gopt.eta,'gsm')
		% from GSM (Maritorena et al. 2002)
		eta = 1.03373;
	else
		% user defined
		eta = gopt.eta;
	end
end


% define Sdg (default to 0.018)
%
sdg = 0.018;
if isfield(gopt,'sdg') == 1
	if strcmp(gopt.sdg,'qaa')
		% from QAA version 5 (Lee et al. 2002, etc.)
		sdg = 0.015 + 0.002 / (0.6 + rin(i443(1))/rin(i555(1)));
	elseif strcmp(gopt.sdg,'obpg')
		% homegrown product from NASA OBPG using NOMAD v2 (unpublished)
		sdg = 0.015 + 0.0038 * log10(rrs(i412(1))/rrs(i555(1)));
		if (sdg <= 0.01); sdg = 0.01; end;
		if (sdg >= 0.02); sdg = 0.02; end;
	elseif strcmp(gopt.sdg,'gsm')
		% from GSM (Maritorena et al. 2002)
		sdg = 0.02061;
	else
		% user defined
		sdg = gopt.sdg;
	end
end


% define aph (default to Bricaud 1998)
%
gopt.aphs = baphs;
if isfield(gopt,'aph') == 1
	if strcmp(gopt.aph,'bricaud')
		% Bricaud 1998 coefficients with aph*(443) normalized to 0.055 m2/mg
		gopt.aphs = baphs;
	elseif strcmp(gopt.aph,'gsm')
		% from GSM (Maritorena et al. 2002)
		gopt.aphs = interp1([412,443,490,510,555,670],[0.006650,0.05582,0.02055,0.0191,0.010150,0.01424],wl,'cubic');
	elseif strcmp(gopt.aph,'ciotti');
		% from Ciotti and Bricaud 2006
		if isfield(gopt,'Sf') ~= 1; gopt.Sf = 0.5; end
		gopt.aphs = get_ciotti_aph(wl,gopt.Sf);
%	elseif strcmp(gopt.aph,'collin');
%		% from Collin Roesler (personal comm.)
%		[dt,dn,ne,nl,tr] = get_collin_aph(wl,'norm');
%		gopt.aphs = dt;
%		gopt.aphse = ne;
%		gopt.aphsl = nl;
	else
		% user defined
		gopt.aphs = gopt.aph;
	end
end


% eigenvectors
%
gopt.bbps = (443 ./ wl).^eta;
gopt.adgs = exp(-sdg .* (wl - 443));


% inversion (default to fminsearch)
%
if isfield(gopt,'inv') == 0; gopt.inv = 'fmin'; end

% nonlinear via fminsearch
%
if strcmp(gopt.inv,'fmin')

	guess = [0.01,0.001,chl];

	opts = optimset('fminsearch');
	opts = optimset(opts,'tolx',1e-6,'tolfun',1e-6,'maxfunevals',1e5,'maxiter',1e3);

	[x,cost,flag] = fminsearch(@(guess) giop_cost(guess,rin,gopt),guess,opts);
	if flag < 1; x = [-999,-999,-999]; end


% linear matrix inversion
%
elseif strcmp(gopt.inv,'lmi')
	
	% positive quadratic root
	%
	q = (-gopt.g0 + sqrt(gopt.g0^2 + 4. * gopt.g1 * rin)) / (2.* gopt.g1);

	% set up matrix
	%
	b = [gopt.bbw .* (1 - q) - gopt.aw .* q]';
	A = [[gopt.adgs .* q]; [gopt.bbps .* (q - 1)]; [gopt.aphs .* q]]';

	% straight up matrix inversion
	%
	%x = A \ b;	

	% QR decomposition
	%
	[Q R] = qr(A);
    x = R \ (R' \ (A' * b));
    r = b - A * x;
    err = R \ (R' \ (A' * r));
    x = x + err;
			
end


% reconstruct spectra
%
madg = x(1) .* gopt.adgs;
mbbp = x(2) .* gopt.bbps;
maph = x(3) .* gopt.aphs;

%if strpos(aph_opt,'collin') ne -1 then begin
%	maph = iop[2]*aphs0 + iop[3]*aphs2
%	mchl = iop[2] + iop[3]
%	caph = [[iop[2]*aphs0],[iop[3]*aphs2]]
%endif

mapg = madg + maph;


% forward model (reconstructed) Rrs
%
if isempty(find(x == -999)) == 1

	moda = gopt.aw + madg + maph;
	modb = gopt.bbw + mbbp;
	modx = modb ./ (modb + moda);
	
    mrrs = gopt.g0 .* modx + gopt.g1 .* modx.^2;

end


% apply QC
% apply delta_Rrs requirement and search for valid retrieval ranges
%
if isfield(gopt,'qc') == 1

	v = find(wl >= 400 & wl <= 600.);
  try
  	rtest = abs(mrrs(v) - rin(v)) ./ rin(v);
  catch E
%     fprintf('rtest failed: maximum number of iterations has been reached');
    mrrs = -999;
		mapg = -999;
		maph = -999;
		madg = -999;
		mbbp = -999;
    return
  end;
    
		
	if isempty(find(rtest > gopt.qc)) ~= 1
% 		fprintf('rtest failed: %f\n', max(rtest));
		mrrs = -999;
		mapg = -999;
		maph = -999;
		madg = -999;
		mbbp = -999;
	else
		u = find(maph < -0.005 | maph > 5);
		if isempty(u) ~= 1; maph(u) = -999; end
		
		u = find(madg < -0.005 | madg > 5);
		if isempty(u) ~= 1; madg(u) = -999; end
		
		u = find(mbbp < -0.005 | mbbp > 5);
		if isempty(u) ~= 1; mbbp(u) = -999; end
		
		u = find(mapg < -0.005 | mapg > 5);
		if isempty(u) ~= 1; mapg(u) = -999; end
	end

end

% print eigenvalues for adg, bbp, and aph
 x;

% clf;
% plot(wl,rin,'r',wl,mrrs,'g');
	


