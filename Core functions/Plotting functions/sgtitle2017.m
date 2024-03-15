function out=sgtitle2017(Tit,varargin)
%SGTITLE Summary of this function goes here
%   Detailed explanation goes here
fprintf('\nbegan sgtitle')
a=ver;
if datetime(a(1).Date)<datetime(2018,1,1)
    out=title(Tit);
else
    out=sgtitle(Tit,varargin{:});
end
fprintf('\nend sgtitle')
end

