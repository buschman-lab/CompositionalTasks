function [out] = ConvertAngleDist(in)
%CONVERTANGLEDIST Convert negative angular distance into positive distance

   out=arrayfun(@(x) SwapAngle(x),in,'UniformOutput',1);

end
function ang=SwapAngle(ang)


if ang>=-pi & ang<0
    ang=-ang;
elseif ang<-pi
    ang=ang+2*pi;
end

end