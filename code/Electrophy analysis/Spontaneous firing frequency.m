clear all;clc

fid=dir('SFCh1_*.mat');

%%

xlstitle='Result_BaseF_30s.xls';

for k=1:length(fid)    
%     clearvars -except fid k xlstitle
    
    filename=fid(k).name;
    load(filename)
    disp([num2str(k),'/',num2str(length(fid)),'--',filename])
    
    V=Ch3.values;%target channel
    dtV=Ch3.interval;
     %V=V*20;%%
    duration=30;      
 peaktime=peakfinder(V, -60, -10, 1, 0)*dtV;%peakfinder(x0, sel, thresh, extrema, include_endpoints)
%     peaktimereal=find(diff(peaktime)>0.0005);%>0.0005    peaktime=peaktime([1;peaktimereal(:)+1]);
    peaktimeplot_idx=fix(peaktime/dtV);
    
    instantFreq=[];
for i=2:length(peaktime)
    ISI=peaktime(i)-peaktime(i-1);   
    instantfreq=1/ISI;    
    instantFreq=[instantFreq;instantfreq];
end
APsd=std(instantFreq);

%% spikewaveform parameters
% tspanV=linspace(0,numel(V)*dtV,numel(V));
%  [spkwaveform,ttspanspk,HW, AHP, Threshold, Amp,max_slope,min_slope,...
%                     HWidx,AHPidx, Thresholdidx,max_slopeidx,min_slopeidx] = AP_Statistic(V,peaktime,tspanV,0);
%               
%                 figure(1),clf
%                 plot(ttspanspk,spkwaveform)
%                 pause
%%
FilenameRecocrd{k,1}=filename;
paramRecord(k,:)=[length(peaktime)/duration,APsd,APsd/mean(instantFreq)]

%%
figure(1),clf
set(gcf,'position',[200,100,800,400])

time_V=(1:length(V))*dtV;

h1=subplot(211);title(filename);axis tight;hold on;
yyaxis left
plot(peaktime(2:end),instantFreq,'.','markerfacecolor',[60,176,106]/255);
ylabel('instantFreq (Hz)')

h2=subplot(212);axis tight;hold on;
plot(time_V,V,'color',[220,20,60]/255);
plot(time_V(peaktimeplot_idx),V(peaktimeplot_idx),'o','markerfacecolor',[60,176,106]/255);

pause()

end

xlswrite(xlstitle,{'Name','basepeakFreq','SD','CV','instantFreq'},'para','A1');
xlswrite(xlstitle,FilenameRecocrd,'para','A2');
xlswrite(xlstitle,paramRecord,'para','B2');
