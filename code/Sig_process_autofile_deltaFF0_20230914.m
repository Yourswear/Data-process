%！！！程序运行中由于需要滤波操作，如果峰值位于前5％或者最后5％会被滤除，使用时需注意
clc;
clear;
close all;
tic
%%
    fs=10;%采样频率，单位Hz
    r_um=50;%用来框选的圆的半径,单位um
    pixsize=0.425;%像素尺寸大小，单位um/pix
    thresh=0.3;%colormap中最大值白色对应的ΔF/F0
    thresh_select=0.12;%荧光区域选取阈值
    bp1=1.5;%拟合的断点1
    bp2=5;%拟合的断点2
%%%%%%%%%%%%%%%%%%%%%%%以上为可调参数%%%%%%%%%%%%%%%%%%%%%%%
folderpath=uigetdir;
namelist = dir([folderpath,'\*.tif']);
r_pix=r_um/pixsize;
dt=1/fs;
dtt=1/fs/30;
%%
for n=1:length(namelist)
    %图像读取
    Info=imfinfo([namelist(n).folder,'\',namelist(n).name]);
    format=Info.Format;
    Slice=size(Info,1);                                          %%获取图片z向帧数
    Width=Info.Width;
    Height=Info.Height;
    t=0:dt:(Slice-1)/fs;
    disp('Information');
    toc
    Image = double(tiffreadVolume([namelist(n).folder,'\',namelist(n).name]));
    disp('Read Image ready');
    toc
    
    %%
    %最大荧光变化展示-高精度版
    t_img_bp=bp1/dt:bp2/dt-dt;
    bas=zeros(Height,Width,length(t_img_bp));
    bas_full=zeros(Height,Width,Slice);
    bas=Image(:,:,bp1/dt+1:bp2/dt);
    for i=1:Height
        for j=1:Width
            dif=(bas(i,j,end)-bas(i,j,1))./((bp2-bp1)/dt-1);
            bas_start=bas(i,j,1)-bp1/dt*dif;
            bas_end=bas(i,j,end)+(Slice-bp2/dt-1)*dif;
            temp=zeros(1,1,Slice);
            temp(1,1,:)=linspace(bas_start,bas_end,Slice);  %将bas的长度由两个bp点之间的距离大小延长到整个图像帧数大小
            bas_full(i,j,:)=temp;
        end
    end
    % sumsig=squeeze(sum(Image(:,:,bp1/dt+1:bp2/dt),[1 2]));  %对图像按二维求和得到总的时间信号用于计算荧光猝灭导致的基线衰减的斜率
    % sumsig_trend=sumsig-detrend(sumsig);
    % dif=(sumsig_trend(end)-sumsig_trend(1))./((bp2-bp1)/dt-1)./(Height*Width);  %计算基线远端两点的差值，相当于计算其斜率，用于后续将基线延长到整个图像
    % 
    % 
    % for i=1:Height
    %     for j=1:Width
    %         bas(i,j,:)=squeeze(Image(i,j,bp1/dt+1:bp2/dt))-detrend(squeeze(Image(i,j,bp1/dt+1:bp2/dt)));  %提取每个像素点的两个bp点之间的基线
    %     end
    % end
    % bas_start=bas(:,:,1)-bp1/dt*dif;
    % bas_end=bas(:,:,end)+(Slice-bp2/dt-1)*dif;
    % temp=zeros(1,1,Slice);
    % 
    % for i=1:Height
    %     for j=1:Width
    %         temp(1,1,:)=linspace(bas_start(i,j),bas_end(i,j),Slice);  %将bas的长度由两个bp点之间的距离大小延长到整个图像帧数大小
    %         bas_full(i,j,:)=temp;
    %     end
    % end
    clear temp;
    delta_Image=(Image-bas_full)./bas_full;  %delta_Image为荧光变化量stack
    maximg=max(delta_Image(:,:,5:14/dt),[],3); %看5~14s内荧光变化量的最大值(即最大刺激程度)
    maximg(maximg<thresh_select)=0;
    mask=im2bw(maximg,0);
    pixnum=length(find(mask>0));
    
    %最大荧光变化展示-快速版
    % aver_50=mean(Image(:,:,10:50),3);%取前50张的平均作为用来减去的background
    % subbackimg=Image-aver_50;%减去背景后的stack
    % maximg=max(subbackimg,[],3);
    % maximg=maximg/mean(aver_50,'all');
    % maximg(maximg<thresh_select)=0;
    % mask=im2bw(maximg,0);
    % pixnum=length(find(mask>0));
    
    maximg=maximg./thresh;
    imwrite(maximg,[namelist(n).folder,'\max_',namelist(n).name]);
    figure(666);
    set(figure(666),'position',[100,200,Height,Width]);
    imshow(maximg);
    mymap=[0         0         0
        0.1137         0    0.1765
        0.2275         0    0.3569
        0.3412         0    0.5333
        0.2275         0    0.6902
        0.1137         0    0.8431
             0         0    1.0000
        0.0118    0.2353    0.6863
        0.0196    0.4667    0.3725
        0.0314    0.7020    0.0588
        0.5176    0.8510    0.0314
        1.0000    1.0000         0
        1.0000    0.6667         0
        1.0000    0.3333         0
        1.0000         0         0
        0.9804    0.1059    0.3137
        0.9608    0.2078    0.6275
        0.9412    0.3137    0.9412
        0.9608    0.5412    0.9608
        0.9804    0.7725    0.9804
        1.0000    1.0000    1.0000];
    colormap(mymap);
    labeldis=linspace(0,1,11);
    label=linspace(0,thresh,11);
    bar=colorbar('Ticks',labeldis,'TickLabels',label);
    bar.Label.String='ΔF / F0';
    disp('Image background extraction');
    toc
    
    %%
    %具体数据处理，包含信号提取，拟合曲线，滤波，目标函数计算以及绘图六个部分
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%根据框选的区域进行信号提取%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Signal=zeros(1,Slice);
    Signal=squeeze(sum(Image,[1 2]));%信号提取
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%拟合曲线%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fitob=fit(t',Signal,'smoothingspline','SmoothingParam',0.99);%拟合曲线增大数据量，增加计算结果精确度,fitob为拟合之后的曲线模型
    max_t0=find(Signal==max(Signal));%原始信号的峰值位置
    fitob2=fit(t(max_t0:end)',Signal(max_t0:end),'smoothingspline','SmoothingParam',0.5);
    tt=0:dtt:Slice*dt-dtt;
    tt_2=bp1:dtt:bp2-dtt;
    lent=length(tt);
    yy=feval(fitob,tt');
    max_t=find(yy==max(yy));%拟合后信号的峰值位置
    yy(max_t:end)=feval(fitob2,tt(max_t:end));%峰值靠后的信号用更平滑的曲线拟合
    %用两个断点bp1,bp2间的信号来拟合基线
    bp_Sig=yy(bp1/dtt+1:bp2/dtt);
    Cor_Sig=detrend(bp_Sig);
    bias=yy(bp1/dtt+1:bp2/dtt)-Cor_Sig; %拟合出的衰减曲线
    fitob3=fit(tt_2',bias,'smoothingspline','SmoothingParam',1);
    bias=feval(fitob3,tt');%拟合出来的基线bias
    %还原矫正后的信号：
    Cor_Sig=100.*(yy-bias)./bias;
    maxsig=max(Cor_Sig(5/dtt:20/dtt));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%绘图%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1);
    set(figure(1),'position',[380,270,length(t)*6,450]);%figure窗口位置及大小调整
    subplot(1,2,2);  
    sigmin=min(Cor_Sig,[],'all');
    sigmax=max(Cor_Sig,[],'all');
    ylim([sigmin sigmax]);
    xlabel('Time/s');
    ylabel('ΔF/F0 %');
    plot(tt,Cor_Sig,'Color','r','LineWidth',1.5);hold on;
    text(0.5*tt(end),0.9*max(Cor_Sig),['max ΔF/F0 = ',num2str(maxsig),'%']);
    legend('ΔF/F0');  hold off;
    subplot(1,2,1); 
    plot(t,Signal,'Color',"b");hold on;
    xlabel('Time/s');
    ylabel('Pixvalue');
    legend('raw curve');  hold off;
    res=['max ΔF/F0 = ',num2str(maxsig),'%'];
    toc
    disp(['percentage of pix exceeding the threshold is: ',num2str(100*(pixnum/Height/Width)),'%']);
    disp(res);
    result(n).name=namelist(n).name;
    result(n).Percentage=100*(pixnum/Height/Width);
    result(n).maxdff=maxsig;
end
save([folderpath,'\result.mat'],"result");