clear; 
close all; 
clc
set(0,'defaultTextInterpreter','latex')

ION="07";
B=importdata("C:\Users\mattb\Documents\Mathematica\Try_3\LevelDensity\Matlab_Input\Sn"+ION+"G.mat");
G=B.G;
Gap=B.Gap;
JPi=B.JPi;
LG=B.L;
dG=B.dG;
V=B.V;
VD=B.VD;
D_WD=B.D_WD;
dD_WD=B.dD_WD;

dVD=VD.*dD_WD./D_WD;

clear B
%%
xHighest=2;
if xHighest==1
else
    for i=1:xHighest-1
        Gap=Gap(Gap~=max(Gap));
    end
end
MAX=max(Gap)
Gmask=Gap==MAX;
GJ=Gap(Gmask);
JP=JPi(Gmask);
JG=G(Gmask);
JdG=dG(Gmask);
JV=V(Gmask);
JVD=VD(Gmask);
JdVD=dVD(Gmask);
JD_WD=D_WD(Gmask);
JdD_WD=dD_WD(Gmask);



i5=0;
i6=0;

for i4=1:length(JP)
    JPI=char(JP(i4));
    if JPI(end)=="o";
        i5=i5+1;
        GJodd(i5)=JG(i4);
        VJodd(i5)=JV(i4);
        VDJodd(i5)=JVD(i4);
        DJodd(i5)=JD_WD(i4);
        if JdG(i4)>20;
            dGJodd(i5)=NaN;
        else
            dGJodd(i5)=JdG(i4);
        end
        if JdVD(i4)>20;
            dVDJodd(i5)=NaN;
        else
            dVDJodd(i5)=JdVD(i4);
        end
        if JdD_WD(i4)>20;
            dDJodd(i5)=NaN;
        else
            dDJodd(i5)=JdD_WD(i4);
        end
        GJoddE(i5)=str2num(JPI(1:end-1));
    else
        i6=i6+1;
        GJeven(i6)=JG(i4);
        VJeven(i6)=JV(i4);
        VDJeven(i6)=JVD(i4);
        DJeven(i6)=JD_WD(i4);
        if JdG(i4)>20;
            dGJeven(i6)=NaN;
        else
            dGJeven(i6)=JdG(i4);
        end
        if JdVD(i4)>20;
            dVDJeven(i6)=NaN;
        else
            dVDJeven(i6)=JdVD(i4);
        end
        if JdD_WD(i4)>20;
            dDJeven(i6)=NaN;
        else
            dDJeven(i6)=JdD_WD(i4);
        end
        GJevenE(i6)=str2num(JPI(1:end-1));
    end
end


MS=10


%%
figure(1)
FONT=12;
xmax=13.99

s(1)=subplot(2,2,2)
ax=gca;
ax.YScale='linear'
grid on
grid minor
hold on
% plot(GJevenE/2,GJeven,'k.','MarkerSize',8)
plot(GJevenE/2,VJeven,'k.','MarkerSize',MS)
plot(GJoddE/2,VJodd,'b.','MarkerSize',MS)
ylabel('V (eV)')
xlabel('J')
hold off
ylim([0,1.5]);
xlim([0,xmax]);
% pbaspect([3 1 1])
set(gcf,'Units','Normalized','Position',[0.1,0.1,0.3,0.3])
ax.XAxis.FontSize = FONT;
ax.YAxis.FontSize = FONT;

s(2)=subplot(2,2,4)
ax=gca;
ax.XAxis.FontSize = FONT;
ax.YAxis.FontSize = FONT;
grid on
grid minor
hold on
% plot(GJevenE/2,GJeven,'k.','MarkerSize',8)
errorbar(GJevenE/2,VDJeven,dVDJeven,'k.','MarkerSize',MS)
errorbar(GJoddE/2,VDJodd,dVDJodd,'b.','MarkerSize',MS)
ylabel('V/D')
xlabel('J')
hold off
ylim([0,30]);
xlim([0,xmax]);

s(3)=subplot(2,2,3)
ax=gca;
ax.YScale='log'
grid on
grid minor
hold on
% plot(GJevenE/2,GJeven,'k.','MarkerSize',8)
errorbar(GJevenE/2,DJeven,dDJeven,'k.','MarkerSize',MS)
errorbar(GJoddE/2,DJodd,dDJodd,'b.','MarkerSize',MS)
ylabel('D (eV)')
xlabel('J')
hold off
ylim([0,1.5]);
xlim([0,xmax]);
% pbaspect([3 1 1])
set(gcf,'Units','Normalized','Position',[0.1,0.1,0.3,0.3])
ax.XAxis.FontSize = FONT;
ax.YAxis.FontSize = FONT;

s(4)=subplot(2,2,1)

ax=gca;
grid on
grid minor
hold on
% plot(GJevenE/2,GJeven,'k.','MarkerSize',8)
errorbar(GJevenE/2,GJeven,dGJeven,'k.','MarkerSize',MS)
errorbar(GJoddE/2,GJodd,dGJodd,'b.','MarkerSize',MS)
ylabel('$\Gamma$ (eV)')
xlabel('J')
hold off
ylim([0,30]);
xlim([0,xmax]);
% pbaspect([3 1 1])
set(gcf,'Units','Normalized','Position',[0.1,0.1,0.3,0.3])
ax.XAxis.FontSize = FONT;
ax.YAxis.FontSize = FONT;


set(gcf,'Units','Normalized','Position',[0.1,0.1,0.5,0.575])
set(s(1),'position',[0.07 0.57 .42 .40])
set(s(2),'position',[0.555 0.57 .42 .40])
set(s(3),'position',[0.07 0.08 .42 .40])
set(s(4),'position',[0.555 0.08 .42 .40])

%%
path="C:/Users/mattb/Documents/Mathematica/Try_3/LevelDensity/Chaos_Param_Images"
saveas(figure(1),[path + '/' + ION + 'All.png'])
