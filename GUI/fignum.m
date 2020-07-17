function num = fignum(hf)
% FIGNUM Return figure number from figure handle - fix for R2014b
%
% num = fignum(hf)
%

if verLessThan('matlab','8.4.0')
    num = hf ;
else
    num = hf.Number ;
end
