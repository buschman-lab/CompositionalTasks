function [Cos,CosDist,AngSim]=CosineSimilarity(resp,targ)
    Cos =  sum(resp(:).*targ(:))./(sqrt(sum(resp(:).^2))*sqrt(sum(targ(:).^2)));
    CosDist=1-Cos;
    AngSim=1-(2*acos(Cos))/pi;
end