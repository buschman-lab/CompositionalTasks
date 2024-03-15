function out=subtitle2017(Tit,varargin)
%SGTITLE Summary of this function goes here
%   Detailed explanation goes here
fprintf('\nbegan subtitle')
a=ver;
if datetime(a(1).Date)<datetime(2020,1,1)
    out=title(Tit);
else
    out=subtitle(Tit,varargin{:});
end
fprintf('\nend subtitle')

end

