if(~exist('fgm', 'dir'))
    disp Please Download FGM package from http://www.f-zhou.com/gm_code.html
    return
end
if(~exist('Cars_and_Motorbikes_Dataset_and_Code','dir'))
    disp Please Download Cars and Motorbikes Dataset from https://sites.google.com/site/graphmatchingmethods/
end
cd fgm
make
cd ../HBPMex 
!make
cd ..