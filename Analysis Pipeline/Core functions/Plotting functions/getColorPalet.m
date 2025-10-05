function kColors = getColorPalet(K)

%Grabs K colors from predetermined color palet
color_palet = [[0 .6 .3]; [.7 0 .7]; [1 .6 0];  [.1 .3 .9];  [1 .1 .1];...
    [0 .9 .3]; [.4 .2 .7]; [.7 .2 .1]; [.1 .8 1 ]; [1 .3 .7]; [.2 .8 .2];...
    [.7 .4 1]; [.9 .6 .4]; [0 .6 1]; [1 .1 .3]]; 
color_palet = repmat(color_palet, ceil(K/size(color_palet,1)),1); 
kColors = color_palet(1:K,:);

end