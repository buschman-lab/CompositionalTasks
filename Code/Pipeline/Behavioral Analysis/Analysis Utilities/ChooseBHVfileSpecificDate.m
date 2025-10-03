%%% load files based on their dates
function BHV=ChooseBHVfileSpecificDate(Kidname,opts,Year,Month,Day,Block) 

cd(opts.PATH);
ip=1; 
for DD=[Day] 
    
    A=num2str(padarray(DD,[0 2-length(num2str(DD))],0,'pre'));    HH=A(~isspace(A));  %day
    B=num2str(padarray(Block(ip),[0 2-length(num2str(Block(ip)))],0,'pre'));BB=B(~isspace(B)); %block    
 %   Y=num2str(padarray(Year(ip),[0 2-length(num2str(Year(ip)))],0,'pre')); YY=Y(~isspace(Y));%year
    Y=num2str(Year);YY=Y(end-1:end);
    M=num2str(padarray(Month(ip),[0 2-length(num2str(Month(ip)))],0,'pre')); MM=M(~isspace(M));%month

    FileNames=[Kidname '_' YY MM  num2str(HH) '_' num2str(BB) '_bhv.mat'];
    disp(['Loading ....' FileNames ])
    BHV{ip}=load([ FileNames],'bhv');
    BHV{ip}.bhv.FigNum=ip;
    ip=ip+1;
    close all
end
cd(opts.CodePath);