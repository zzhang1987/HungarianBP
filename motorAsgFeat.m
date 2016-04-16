function wsFeat = motorAsgFeat(wsSrc, parG, parF, varargin)
prex = wsSrc.prex;
Pts = wsSrc.Pts;
mG = length(Pts);
gphs = ps(wsSrc, 'gphs', []);
if isempty(gphs)
    gphs = newGphAs(Pts, parG);
end

XPs{1} = wsSrc.Features{1}(wsSrc.ords{1},9)';
XPs{2} = wsSrc.Features{2}(wsSrc.ords{2},9)';
gphs{1}.XP = XPs{1};
gphs{2}.XP = XPs{2};
wsFeat.gphs = gphs;
wsFeat.XPs = XPs;
wsFeat.Fs = wsSrc.Fs;

