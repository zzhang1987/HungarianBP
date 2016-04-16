if(~exist('HBPMex','dir'))
    !git clone https://github.com/zzhang1987/HungarianBPMex HBPMex
else
    cd HBPMex
    !git pull 
    cd ..
end
if(~exist('fgm', 'dir'))
    disp Please Download FGM package from http://www.f-zhou.com/gm_code.html
    disp Please uncompress FGM package in the same directory as the script.
    return
end
if(~exist('Cars_and_Motorbikes_Dataset_and_Code','dir'))
    disp Please Download Cars and Motorbikes Dataset from https://sites.google.com/site/graphmatchingmethods/
    disp Please uncompress Cars and Motorbikes Dataset in the same directory as the script.
end
cd fgm
make
cd ../HBPMex 
!make
cd ..