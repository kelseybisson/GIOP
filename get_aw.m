function [aw] = get_aw(wl)

all = load('optics_coef.txt');

aw = interp1(all(:,1),all(:,2),wl,'cubic');

