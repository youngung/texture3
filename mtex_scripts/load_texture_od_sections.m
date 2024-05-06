filename='/Users/youngung/repo/evpsc/umat/src/evpsc_incr/matData/lc/texinlc.tex'
filename='/Users/youngung/repo/vpsc8_ReX/TEX_PH2.OUT'

%filename='/Users/youngung/repo/vpsc8_ReX/L_1075C_e0_1.tex'
odf=loadODF_VPSC(filename)
%mtexColorMap('parula')

setMTEXpref('defaultColorMap',parulaColorMap);
%setMTEXpref('defaultColorMap',WhiteJetColorMap);
%setMTEXpref('defaultColorMap',mtexColorMap('parula'));

%plot3d(odf)%,cmap='parula')
%plot3d(odf,'axisAngle','figSize','large')
%plotSection(odf,'phi2','sections',18,'resolution',5*degree)


% set position 1 in a 1x3 matrix ad the current plotting position
t=tiledlayout(4,4)
nexttile
plotSection(odf,'sections',1,'phi2',0*degree )
%f=gcf()
x0=1200; y0=800; dx=250; dy=250
%f.Position=[x0 y0 dx dy];colorbar()
saveFigure('/tmp/phi2=0.png')
nexttile
plotSection(odf,'sections',1,'phi2',45*degree)
%f=gcf()
x0=1200; y0=800; dx=250; dy=250
%f.Position=[x0 y0 dx dy]
saveFigure('/tmp/phi2=45.png')
nexttile
plotSection(odf,'sections',1,'phi2',60*degree)
%f=gcf()
x0=1200; y0=800; dx=250; dy=250
%f.Position=[x0 y0 dx dy]
saveFigure('/tmp/phi2=60.png')