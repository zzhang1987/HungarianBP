if(~exist('fgm', 'dir'))
    disp Please Download FGM package from http://www.f-zhou.com/gm_code.html
    disp Please uncompress FGM package in the same directory as the script.
    return
end
if(~exist('Cars_and_Motorbikes_Graph_Matching_Datasets_and_Code','dir'))
    disp Please Download Cars and Motorbikes Dataset from https://sites.google.com/site/graphmatchingmethods/
    disp Please uncompress Cars and Motorbikes Dataset in the same directory as the script.
end
if(~exist('data_chrct','dir'))
    disp Please Download Chinese Character Data set from http://www.escience.cn/system/file?fileId=62549
    disp Please uncompress the downloaded dataset in the same directory as the script.
end
cd fgm
make
cd ..
% cd HBPMex
% mex -v -O -output HungarianBPMex.mexa64 HungarianBPMex.cpp BPSolver.cpp DualSolver.cpp PiSBPSolver.cpp PSBPSolver.cpp HungarianBP.cpp
% mex -v -O -output CluComputeObjMex.mexa64 CluComputeObjMex.cpp
% cd ..
