function [bbw] = get_bbw(wl)

bbw = 0.0038 .* (400 ./ wl).^4.32;


