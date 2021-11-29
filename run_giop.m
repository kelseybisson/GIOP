clear all;

% GIOP ocean color reflectance inversion model
%
% P.J. Werdell and 18 co-authors, "Generalized ocean color 
% inversion model for retrieving marine inherent optical 
% properties," Appl. Opt. 52, 2019-2037 (2013).
%

% GIOP is an ocean reflectance inversion model that
% can be configured at run-time 
%
% requires equally-sized vectors of wavelength and Rrs
% plus an estimate of chlorophyll; all other parameterizations
% are controlled by the structure 'gopt', described below
%
% processing comments:
%
% - defaults to GIOP-DC configuration 
% - currently requires 412, 443, and 547/555 nm to be present
%
% outputs are a vector of the magnitudes of the eigenvalues
% for adg, bbp, and aph (x), plus modeled spectra of apg,
% aph, adg, bbp, and Rrs
%
% Jeremy Werdell, NASA Goddard Space Flight Center, July 2013



% here are example vectors of wavelength and Rrs
% plus an estimate of chlorophyll from the OC4v6 algorithm
%
wl = [412,443,490,510,555,670];
rrs = [0.003478,0.004074,0.004465,0.003588,0.002494,0.000051];
oc = get_oc(rrs(2),rrs(3),rrs(4),rrs(5),'oc4')


% build the parameterization control structure 'gopt'
% will default to GIOP-DC if not passed to giop function
 		  
% air-sea transmission (defaults to Lee et al. 2002)
%gopt.trans='flat'; % the spectrally flat Gordon 2007 relationship ***

% pure seawater absorption and backscattering (defaults to Pope and Fry 1997
% for absorption and a fit to Morel 1974 for backscattering)  	  	  
gopt.aw  = get_aw(wl);
gopt.bbw = get_bbw(wl);

% solar and sensor viewing geometries (defaults to solz=30, senz=0, relaz=0)
% only required when AOP-IOP relationship of Morel et al. 2002 is used
%gopt.solz = 30;
%gopt.senz = 0.001;
%gopt.relaz = 0.001;

% quality control; maximum Rrs_measured - Rrs_modeled spectral difference
% ranges from 0 to 1 and considers 400-600 nm
gopt.qc = 0.33; % = 33%

% AOP-IOP relationship (GIOP-DC default is Gordon et al. 1988)
%gopt.fq = 'morel'; % Morel et al. 2002
%gopt.fq = 'gordon'; % Gordon et al. 1988
%gopt.fq = [A,B]; % user-defined (provide A and B)

% define eta (power-law spectral shape of bbp) (GIOP-DC default QAA)
%gopt.eta = 'qaa'; % Lee et al. 2002 (QAA)
%gopt.eta = 1.03373; % Maritorena et al. 2002 (GSM)
%gopt.eta = A; % user-defined (provide A)

% define Sdg (exponent spectral shape of adg) (GIOP-DC default is 0.018)
%gopt.sdg = 0.018; % GIOP-DC default 
%gopt.sdg = 'qaa'; % Lee et al. 2002 (QAA)
%gopt.sdg = 0.02061; % Maritorena et al. 2002 (GSM)
%gopt.sdg = 'obpg'; % homegrown product from NASA OBPG using NOMAD v2 (unpublished)
%gopt.sdg = A; % user-defined (provide A)

% define aph* (mass-specific absorption of chlorophyll) (GIOP-DC default is Bricaud)
%gopt.aph = 'bricaud'; % Bricaud et al. 1998
%gopt.aph = 'gsm'; % Maritorena et al. 2002 (GSM)
%giop.aph = 'ciotti'; % Ciotti and Bricaud 2006

% define Sf (phytoplankton size fraction from 0-1; use only with Ciotti and Bricaud aph*)
%gopt.Sf = 0.5;

% define statistical inversion method (GIOP-DC default is non-linear 'fmin')
%gopt.inv = 'fmin'; % non-linear inversion using fminsearch Matlab routine
%gopt.inv = 'lmi'; % linear matrix inversion using QR decomposition

% should report 0.0441, 0.0033, and 0.3693 to screen as
% eigenvectors for adg, bbp, and aph
[x,apg,aph,adg,bbp,mrrs] = giop(wl,rrs,oc(1),gopt);
% x is a vector of eigenvalues for adg, bbp, and aph
% apg is modeled adg + aph
% mrrs is modeled (reconstructed) Rrs

figure(1)
plot(wl,rrs,'k-',wl,mrrs,'b--'),ylabel('R_{rs}(\lambda) (sr^{-1})'),xlabel('Wavelength (nm)'),legend('measured','modeled')
figure(2)
plot(wl,aph,'g-',wl,adg,'c-'),ylabel('a(\lambda) (m^{-1})'),xlabel('Wavelength (nm)'),legend('a_{phyt}','a_{cdm}')
figure(3)
plot(wl,gopt.bbw,'g-',wl,bbp,'c-'),ylabel('b_b(\lambda) (m^{-1})'),xlabel('Wavelength (nm)'),legend('b_{bw}','b_{bp}')

%%% Linear Matrix Inversion Option
% this modification should report 0.0414, 0.0022, and 0.1058 to screen as
% eigenvectors for adg, bbp, and aph
gopt.inv = 'lmi';
[t1,apg,aph,adg,bbp,mrrs] = giop(wl,rrs,oc(1),gopt);




