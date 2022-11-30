function plotMorphologyParamPM_LPSconc(data)
%% PREPARE DATA

%load('dataAnalysis.mat','data');

%SEPARATE DATA ACCORDING TO CONDITION
%Initialise variables.
LPS0=[];
LPS100=[];
LPS250=[];
LPS500=[];
LPS750=[];
LPS1000=[];
LPS2000=[];

for ii=1:length(data)
    if ~isempty(strfind(data{ii,1},'2000LPS'))
        LPS2000=[LPS2000;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'1000LPS'))
        LPS1000=[LPS1000;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'750LPS'))
        LPS750=[LPS750;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'500LPS'))
        LPS500=[LPS500;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'250LPS'))
        LPS250=[LPS250;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'100LPS'))
        LPS100=[LPS100;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'0LPS'))
        LPS0=[LPS0;data{ii,3}];
    end
end

%EXTRACT PARAMETERS FOR EACH CONDITION
%Initialise variables.
LPS0_Area=zeros(length(LPS0),1);
LPS100_Area=zeros(length(LPS100),1);
LPS250_Area=zeros(length(LPS250),1);
LPS500_Area=zeros(length(LPS500),1);
LPS750_Area=zeros(length(LPS750),1);
LPS1000_Area=zeros(length(LPS1000),1);
LPS2000_Area=zeros(length(LPS2000),1);

LPS0_MajorAxisLength=zeros(length(LPS0),1);
LPS100_MajorAxisLength=zeros(length(LPS100),1);
LPS250_MajorAxisLength=zeros(length(LPS250),1);
LPS500_MajorAxisLength=zeros(length(LPS500),1);
LPS750_MajorAxisLength=zeros(length(LPS750),1);
LPS1000_MajorAxisLength=zeros(length(LPS1000),1);
LPS2000_MajorAxisLength=zeros(length(LPS2000),1);

LPS0_MinorAxisLength=zeros(length(LPS0),1);
LPS100_MinorAxisLength=zeros(length(LPS100),1);
LPS250_MinorAxisLength=zeros(length(LPS250),1);
LPS500_MinorAxisLength=zeros(length(LPS500),1);
LPS750_MinorAxisLength=zeros(length(LPS750),1);
LPS1000_MinorAxisLength=zeros(length(LPS1000),1);
LPS2000_MinorAxisLength=zeros(length(LPS2000),1);

LPS0_Eccentricity=zeros(length(LPS0),1);
LPS100_Eccentricity=zeros(length(LPS100),1);
LPS250_Eccentricity=zeros(length(LPS250),1);
LPS500_Eccentricity=zeros(length(LPS500),1);
LPS750_Eccentricity=zeros(length(LPS750),1);
LPS1000_Eccentricity=zeros(length(LPS1000),1);
LPS2000_Eccentricity=zeros(length(LPS2000),1);

LPS0_Circularity=zeros(length(LPS0),1);
LPS100_Circularity=zeros(length(LPS100),1);
LPS250_Circularity=zeros(length(LPS250),1);
LPS500_Circularity=zeros(length(LPS500),1);
LPS750_Circularity=zeros(length(LPS750),1);
LPS1000_Circularity=zeros(length(LPS1000),1);
LPS2000_Circularity=zeros(length(LPS2000),1);

LPS0_Extent=zeros(length(LPS0),1);
LPS100_Extent=zeros(length(LPS100),1);
LPS250_Extent=zeros(length(LPS250),1);
LPS500_Extent=zeros(length(LPS500),1);
LPS750_Extent=zeros(length(LPS750),1);
LPS1000_Extent=zeros(length(LPS1000),1);
LPS2000_Extent=zeros(length(LPS2000),1);

LPS0_Perimeter=zeros(length(LPS0),1);
LPS100_Perimeter=zeros(length(LPS100),1);
LPS250_Perimeter=zeros(length(LPS250),1);
LPS500_Perimeter=zeros(length(LPS500),1);
LPS750_Perimeter=zeros(length(LPS750),1);
LPS1000_Perimeter=zeros(length(LPS1000),1);
LPS2000_Perimeter=zeros(length(LPS2000),1);

LPS0_Elongation=zeros(length(LPS0),1);
LPS100_Elongation=zeros(length(LPS100),1);
LPS250_Elongation=zeros(length(LPS250),1);
LPS500_Elongation=zeros(length(LPS500),1);
LPS750_Elongation=zeros(length(LPS750),1);
LPS1000_Elongation=zeros(length(LPS1000),1);
LPS2000_Elongation=zeros(length(LPS2000),1);

%CONDITION 1 PARAMETERS
for ii=1:length(LPS0)
    LPS0_Area(ii,1)=LPS0(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS0_MajorAxisLength(ii,1)=LPS0(ii).MajorAxisLength*0.18;
    LPS0_MinorAxisLength(ii,1)=LPS0(ii).MinorAxisLength*0.18;
    LPS0_Eccentricity(ii,1)=LPS0(ii).Eccentricity;
    LPS0_Circularity(ii,1)=LPS0(ii).Circularity;
    LPS0_Extent(ii,1)=LPS0(ii).Extent;
    LPS0_Perimeter(ii,1)=LPS0(ii).Perimeter*0.18;
    LPS0_Elongation(ii,1)=(LPS0(ii).MajorAxisLength*0.18)/(LPS0(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 2 PARAMETERS
for ii=1:length(LPS100)
    LPS100_Area(ii,1)=LPS100(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS100_MajorAxisLength(ii,1)=LPS100(ii).MajorAxisLength*0.18;
    LPS100_MinorAxisLength(ii,1)=LPS100(ii).MinorAxisLength*0.18;
    LPS100_Eccentricity(ii,1)=LPS100(ii).Eccentricity;
    LPS100_Circularity(ii,1)=LPS100(ii).Circularity;
    LPS100_Extent(ii,1)=LPS100(ii).Extent;
    LPS100_Perimeter(ii,1)=LPS100(ii).Perimeter*0.18;
    LPS100_Elongation(ii,1)=(LPS100(ii).MajorAxisLength*0.18)/(LPS100(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 3 PARAMETERS
for ii=1:length(LPS250)
    LPS250_Area(ii,1)=LPS250(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS250_MajorAxisLength(ii,1)=LPS250(ii).MajorAxisLength*0.18;
    LPS250_MinorAxisLength(ii,1)=LPS250(ii).MinorAxisLength*0.18;
    LPS250_Eccentricity(ii,1)=LPS250(ii).Eccentricity;
    LPS250_Circularity(ii,1)=LPS250(ii).Circularity;
    LPS250_Extent(ii,1)=LPS250(ii).Extent;
    LPS250_Perimeter(ii,1)=LPS250(ii).Perimeter*0.18;
    LPS250_Elongation(ii,1)=(LPS250(ii).MajorAxisLength*0.18)/(LPS250(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 4 PARAMETERS
for ii=1:length(LPS500)
    LPS500_Area(ii,1)=LPS500(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS500_MajorAxisLength(ii,1)=LPS500(ii).MajorAxisLength*0.18;
    LPS500_MinorAxisLength(ii,1)=LPS500(ii).MinorAxisLength*0.18;
    LPS500_Eccentricity(ii,1)=LPS500(ii).Eccentricity;
    LPS500_Circularity(ii,1)=LPS500(ii).Circularity;
    LPS500_Extent(ii,1)=LPS500(ii).Extent;
    LPS500_Perimeter(ii,1)=LPS500(ii).Perimeter*0.18;
    LPS500_Elongation(ii,1)=(LPS500(ii).MajorAxisLength*0.18)/(LPS500(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 5 PARAMETERS
for ii=1:length(LPS750)
    LPS750_Area(ii,1)=LPS750(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS750_MajorAxisLength(ii,1)=LPS750(ii).MajorAxisLength*0.18;
    LPS750_MinorAxisLength(ii,1)=LPS750(ii).MinorAxisLength*0.18;
    LPS750_Eccentricity(ii,1)=LPS750(ii).Eccentricity;
    LPS750_Circularity(ii,1)=LPS750(ii).Circularity;
    LPS750_Extent(ii,1)=LPS750(ii).Extent;
    LPS750_Perimeter(ii,1)=LPS750(ii).Perimeter*0.18;
    LPS750_Elongation(ii,1)=(LPS750(ii).MajorAxisLength*0.18)/(LPS750(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 6 PARAMETERS
for ii=1:length(LPS1000)
    LPS1000_Area(ii,1)=LPS1000(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS1000_MajorAxisLength(ii,1)=LPS1000(ii).MajorAxisLength*0.18;
    LPS1000_MinorAxisLength(ii,1)=LPS1000(ii).MinorAxisLength*0.18;
    LPS1000_Eccentricity(ii,1)=LPS1000(ii).Eccentricity;
    LPS1000_Circularity(ii,1)=LPS1000(ii).Circularity;
    LPS1000_Extent(ii,1)=LPS1000(ii).Extent;
    LPS1000_Perimeter(ii,1)=LPS1000(ii).Perimeter*0.18;
    LPS1000_Elongation(ii,1)=(LPS1000(ii).MajorAxisLength*0.18)/(LPS1000(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 7 PARAMETERS
for ii=1:length(LPS2000)
    LPS2000_Area(ii,1)=LPS2000(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS2000_MajorAxisLength(ii,1)=LPS2000(ii).MajorAxisLength*0.18;
    LPS2000_MinorAxisLength(ii,1)=LPS2000(ii).MinorAxisLength*0.18;
    LPS2000_Eccentricity(ii,1)=LPS2000(ii).Eccentricity;
    LPS2000_Circularity(ii,1)=LPS2000(ii).Circularity;
    LPS2000_Extent(ii,1)=LPS2000(ii).Extent;
    LPS2000_Perimeter(ii,1)=LPS2000(ii).Perimeter*0.18;
    LPS2000_Elongation(ii,1)=(LPS2000(ii).MajorAxisLength*0.18)/(LPS2000(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%% REMOVE OUTLIERS

%CONDITION 1
for ii=1:length(LPS0)
    LPS0_Area=rmoutliers(LPS0_Area);
    LPS0_MajorAxisLength=rmoutliers(LPS0_MajorAxisLength);
    LPS0_MinorAxisLength=rmoutliers(LPS0_MinorAxisLength);
    LPS0_Eccentricity=rmoutliers(LPS0_Eccentricity);
    LPS0_Circularity=rmoutliers(LPS0_Circularity);
    LPS0_Extent=rmoutliers(LPS0_Extent);
    LPS0_Perimeter=rmoutliers(LPS0_Perimeter);
    LPS0_Elongation=rmoutliers(LPS0_Elongation);
end

%CONDITION 2
for ii=1:length(LPS100)
    LPS100_Area=rmoutliers(LPS100_Area);
    LPS100_MajorAxisLength=rmoutliers(LPS100_MajorAxisLength);
    LPS100_MinorAxisLength=rmoutliers(LPS100_MinorAxisLength);
    LPS100_Eccentricity=rmoutliers(LPS100_Eccentricity);
    LPS100_Circularity=rmoutliers(LPS100_Circularity);
    LPS100_Extent=rmoutliers(LPS100_Extent);
    LPS100_Perimeter=rmoutliers(LPS100_Perimeter);
    LPS100_Elongation=rmoutliers(LPS100_Elongation);
end

%CONDITION 3
for ii=1:length(LPS250)
    LPS250_Area=rmoutliers(LPS250_Area);
    LPS250_MajorAxisLength=rmoutliers(LPS250_MajorAxisLength);
    LPS250_MinorAxisLength=rmoutliers(LPS250_MinorAxisLength);
    LPS250_Eccentricity=rmoutliers(LPS250_Eccentricity);
    LPS250_Circularity=rmoutliers(LPS250_Circularity);
    LPS250_Extent=rmoutliers(LPS250_Extent);
    LPS250_Perimeter=rmoutliers(LPS250_Perimeter);
    LPS250_Elongation=rmoutliers(LPS250_Elongation);
end

%CONDITION 4
for ii=1:length(LPS500)
    LPS500_Area=rmoutliers(LPS500_Area);
    LPS500_MajorAxisLength=rmoutliers(LPS500_MajorAxisLength);
    LPS500_MinorAxisLength=rmoutliers(LPS500_MinorAxisLength);
    LPS500_Eccentricity=rmoutliers(LPS500_Eccentricity);
    LPS500_Circularity=rmoutliers(LPS500_Circularity);
    LPS500_Extent=rmoutliers(LPS500_Extent);
    LPS500_Perimeter=rmoutliers(LPS500_Perimeter);
    LPS500_Elongation=rmoutliers(LPS500_Elongation);
end

%CONDITION 5
for ii=1:length(LPS750)
    LPS750_Area=rmoutliers(LPS750_Area);
    LPS750_MajorAxisLength=rmoutliers(LPS750_MajorAxisLength);
    LPS750_MinorAxisLength=rmoutliers(LPS750_MinorAxisLength);
    LPS750_Eccentricity=rmoutliers(LPS750_Eccentricity);
    LPS750_Circularity=rmoutliers(LPS750_Circularity);
    LPS750_Extent=rmoutliers(LPS750_Extent);
    LPS750_Perimeter=rmoutliers(LPS750_Perimeter);
    LPS750_Elongation=rmoutliers(LPS750_Elongation);
end

%CONDITION 6
for ii=1:length(LPS1000)
    LPS1000_Area=rmoutliers(LPS1000_Area);
    LPS1000_MajorAxisLength=rmoutliers(LPS1000_MajorAxisLength);
    LPS1000_MinorAxisLength=rmoutliers(LPS1000_MinorAxisLength);
    LPS1000_Eccentricity=rmoutliers(LPS1000_Eccentricity);
    LPS1000_Circularity=rmoutliers(LPS1000_Circularity);
    LPS1000_Extent=rmoutliers(LPS1000_Extent);
    LPS1000_Perimeter=rmoutliers(LPS1000_Perimeter);
    LPS1000_Elongation=rmoutliers(LPS1000_Elongation);
end

%CONDITION 7
for ii=1:length(LPS2000)
    LPS2000_Area=rmoutliers(LPS2000_Area);
    LPS2000_MajorAxisLength=rmoutliers(LPS2000_MajorAxisLength);
    LPS2000_MinorAxisLength=rmoutliers(LPS2000_MinorAxisLength);
    LPS2000_Eccentricity=rmoutliers(LPS2000_Eccentricity);
    LPS2000_Circularity=rmoutliers(LPS2000_Circularity);
    LPS2000_Extent=rmoutliers(LPS2000_Extent);
    LPS2000_Perimeter=rmoutliers(LPS2000_Perimeter);
    LPS2000_Elongation=rmoutliers(LPS2000_Elongation);
end

%% PERFORM STATISTICS FOR EACH PARAMETER

%Prepare data for Kruskall Wallis tests.
g1=ones(length(LPS0_Area),1); %condition 1.
g2=2*ones(length(LPS500_Area),1); %condition 2.
g3=3*ones(length(LPS100_Area),1); %condition 3.
g4=4*ones(length(LPS750_Area),1); %condition 4.
g5=5*ones(length(LPS250_Area),1); %condition 5.
g6=6*ones(length(LPS1000_Area),1); %condition 6.
g7=7*ones(length(LPS2000_Area),1); %condition 7.
group_area=[g1;g2;g3;g4;g5;g6;g7]; %merge conditions in order.

g1=ones(length(LPS0_MajorAxisLength),1); %condition 1.
g2=2*ones(length(LPS500_MajorAxisLength),1); %condition 2.
g3=3*ones(length(LPS100_MajorAxisLength),1); %condition 3.
g4=4*ones(length(LPS750_MajorAxisLength),1); %condition 4.
g5=5*ones(length(LPS250_MajorAxisLength),1); %condition 5.
g6=6*ones(length(LPS1000_MajorAxisLength),1); %condition 6.
g7=7*ones(length(LPS2000_MajorAxisLength),1); %condition 7.
group_majorAxisLength=[g1;g2;g3;g4;g5;g6;g7]; %merge conditions in order.

g1=ones(length(LPS0_MinorAxisLength),1); %condition 1.
g2=2*ones(length(LPS500_MinorAxisLength),1); %condition 2.
g3=3*ones(length(LPS100_MinorAxisLength),1); %condition 3.
g4=4*ones(length(LPS750_MinorAxisLength),1); %condition 4.
g5=5*ones(length(LPS250_MinorAxisLength),1); %condition 5.
g6=6*ones(length(LPS1000_MinorAxisLength),1); %condition 6.
g7=7*ones(length(LPS2000_MinorAxisLength),1); %condition 7.
group_minorAxisLength=[g1;g2;g3;g4;g5;g6;g7]; %merge conditions in order.

g1=ones(length(LPS0_Eccentricity),1); %condition 1.
g2=2*ones(length(LPS500_Eccentricity),1); %condition 2.
g3=3*ones(length(LPS100_Eccentricity),1); %condition 3.
g4=4*ones(length(LPS750_Eccentricity),1); %condition 4.
g5=5*ones(length(LPS250_Eccentricity),1); %condition 5.
g6=6*ones(length(LPS1000_Eccentricity),1); %condition 6.
g7=7*ones(length(LPS2000_Eccentricity),1); %condition 7.
group_eccentricity=[g1;g2;g3;g4;g5;g6;g7]; %merge conditions in order.

g1=ones(length(LPS0_Circularity),1); %condition 1.
g2=2*ones(length(LPS500_Circularity),1); %condition 2.
g3=3*ones(length(LPS100_Circularity),1); %condition 3.
g4=4*ones(length(LPS750_Circularity),1); %condition 4.
g5=5*ones(length(LPS250_Circularity),1); %condition 5.
g6=6*ones(length(LPS1000_Circularity),1); %condition 6.
g7=7*ones(length(LPS2000_Circularity),1); %condition 7.
group_circularity=[g1;g2;g3;g4;g5;g6;g7]; %merge conditions in order.

g1=ones(length(LPS0_Extent),1); %condition 1.
g2=2*ones(length(LPS500_Extent),1); %condition 2.
g3=3*ones(length(LPS100_Extent),1); %condition 3.
g4=4*ones(length(LPS750_Extent),1); %condition 4.
g5=5*ones(length(LPS250_Extent),1); %condition 5.
g6=6*ones(length(LPS1000_Extent),1); %condition 6.
g7=7*ones(length(LPS2000_Extent),1); %condition 7.
group_extent=[g1;g2;g3;g4;g5;g6;g7]; %merge conditions in order.

g1=ones(length(LPS0_Perimeter),1); %condition 1.
g2=2*ones(length(LPS500_Perimeter),1); %condition 2.
g3=3*ones(length(LPS100_Perimeter),1); %condition 3.
g4=4*ones(length(LPS750_Perimeter),1); %condition 4.
g5=5*ones(length(LPS250_Perimeter),1); %condition 5.
g6=6*ones(length(LPS1000_Perimeter),1); %condition 6.
g7=7*ones(length(LPS2000_Perimeter),1); %condition 7.
group_perimeter=[g1;g2;g3;g4;g5;g6;g7]; %merge conditions in order.

g1=ones(length(LPS0_Elongation),1); %condition 1.
g2=2*ones(length(LPS500_Elongation),1); %condition 2.
g3=3*ones(length(LPS100_Elongation),1); %condition 3.
g4=4*ones(length(LPS750_Elongation),1); %condition 4.
g5=5*ones(length(LPS250_Elongation),1); %condition 5.
g6=6*ones(length(LPS1000_Elongation),1); %condition 6.
g7=7*ones(length(LPS2000_Elongation),1); %condition 7.
group_elongation=[g1;g2;g3;g4;g5;g6;g7]; %merge conditions in order.

%Merge data in order.
area=[LPS0_Area;LPS100_Area;LPS250_Area;LPS500_Area;LPS750_Area;LPS1000_Area;LPS2000_Area];
majorAxisLength=[LPS0_MajorAxisLength;LPS100_MajorAxisLength;LPS250_MajorAxisLength;LPS500_MajorAxisLength;LPS750_MajorAxisLength;LPS1000_MajorAxisLength;LPS2000_MajorAxisLength];
minorAxisLength=[LPS0_MinorAxisLength;LPS100_MinorAxisLength;LPS250_MinorAxisLength;LPS500_MinorAxisLength;LPS750_MinorAxisLength;LPS1000_MinorAxisLength;LPS2000_MinorAxisLength];
eccentricity=[LPS0_Eccentricity;LPS100_Eccentricity;LPS250_Eccentricity;LPS500_Eccentricity;LPS750_Eccentricity;LPS1000_Eccentricity;LPS2000_Eccentricity];
circularity=[LPS0_Circularity;LPS100_Circularity;LPS250_Circularity;LPS500_Circularity;LPS750_Circularity;LPS1000_Circularity;LPS2000_Circularity];
extent=[LPS0_Extent;LPS100_Extent;LPS250_Extent;LPS500_Extent;LPS750_Extent;LPS1000_Extent;LPS2000_Extent];
perimeter=[LPS0_Perimeter;LPS100_Perimeter;LPS250_Perimeter;LPS500_Perimeter;LPS750_Perimeter;LPS1000_Perimeter;LPS2000_Perimeter];
elongation=[LPS0_Elongation;LPS100_Elongation;LPS250_Elongation;LPS500_Elongation;LPS750_Elongation;LPS1000_Elongation;LPS2000_Elongation];

%Perform Kruskal-Wallis test.
[p_area,~,stats_area]=kruskalwallis(area,group_area);
[p_majorAxisLength,~,stats_majorAxisLength]=kruskalwallis(majorAxisLength,group_majorAxisLength);
[p_minorAxisLength,~,stats_minorAxisLength]=kruskalwallis(minorAxisLength,group_minorAxisLength);
[p_eccentricity,~,stats_eccentricity]=kruskalwallis(eccentricity,group_eccentricity);
[p_circularity,~,stats_circularity]=kruskalwallis(circularity,group_circularity);
[p_extent,~,stats_extent]=kruskalwallis(extent,group_extent);
[p_perimeter,~,stats_perimeter]=kruskalwallis(perimeter,group_perimeter);
[p_elongation,~,stats_elongation]=kruskalwallis(elongation,group_elongation);
close all

%Perform multiple comparison tests for all parameters.
c_area=multcompare(stats_area,'CType','hsd','Display','off'); %first two columns are the groups being compared. Last column is the p-value.
c_majorAxisLength=multcompare(stats_majorAxisLength,'CType','hsd','Display','off');
c_minorAxisLength=multcompare(stats_minorAxisLength,'CType','hsd','Display','off');
c_eccentricity=multcompare(stats_eccentricity,'CType','hsd','Display','off');
c_circularity=multcompare(stats_circularity,'CType','hsd','Display','off');
c_extent=multcompare(stats_extent,'CType','hsd','Display','off');
c_perimeter=multcompare(stats_perimeter,'CType','hsd','Display','off');
c_elongation=multcompare(stats_elongation,'CType','hsd','Display','off');

save('KruskalWallisResults.mat','p_area','p_majorAxisLength','p_minorAxisLength','p_eccentricity','p_circularity','p_extent','p_perimeter','p_elongation');
save('multCompResults.mat','c_area','c_majorAxisLength','c_minorAxisLength','c_eccentricity','c_circularity','c_extent','c_perimeter','c_elongation');

%% PLOT AREA FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Area));
x2=4*ones(1,length(LPS100_Area));
x3=6*ones(1,length(LPS250_Area));
x4=8*ones(1,length(LPS500_Area));
x5=10*ones(1,length(LPS750_Area));
x6=12*ones(1,length(LPS1000_Area));
x7=14*ones(1,length(LPS2000_Area));

plot(x1,LPS0_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Area),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS100_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS100_Area),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS250_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS250_Area),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS500_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS500_Area),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x5,LPS750_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 5.
plot(x5,mean(LPS750_Area),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS1000_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS1000_Area),'kx','MarkerSize',18); %plot mean for condition 6.
plot(x7,LPS2000_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 7.
plot(x7,mean(LPS2000_Area),'kx','MarkerSize',18); %plot mean for condition 7.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 15]);
xticks(1:15);
ylabel('Cell area (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','0','','100','','250','','500','','750','','1000','','2000',''});

%Indicate significance on plot, if any.
sig=c_area(:,6)<0.05; %extract only significant results.
sig_area=c_area(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_area,1)
    if sig_area(ii,6)<0.05
        if sig_area(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_area(ii,1)==2
            x1=4;
        elseif sig_area(ii,1)==3
            x1=6;
        elseif sig_area(ii,1)==4
            x1=8;
        elseif sig_area(ii,1)==5
            x1=10;
        elseif sig_area(ii,1)==6
            x1=12;
        else
            x1=14;
        end
        if sig_area(ii,2)==2
            x2=4;
        elseif sig_area(ii,2)==3
            x2=6;
        elseif sig_area(ii,2)==4
            x2=8;
        elseif sig_area(ii,2)==5
            x2=10;
        elseif sig_area(ii,2)==6
            x2=12;
        else
            x2=14;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_area(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_area(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_area(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end


%Save plot.
savefig(fig,'CellArea'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'CellArea.svg'); %save figure as an SVG file.
saveas(fig,'CellArea.jpg'); %save figure as an JPEG file.
close

%% PLOT MAJOR AXIS LENGTH FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_MajorAxisLength));
x2=4*ones(1,length(LPS100_MajorAxisLength));
x3=6*ones(1,length(LPS250_MajorAxisLength));
x4=8*ones(1,length(LPS500_MajorAxisLength));
x5=10*ones(1,length(LPS750_MajorAxisLength));
x6=12*ones(1,length(LPS1000_MajorAxisLength));
x7=14*ones(1,length(LPS2000_MajorAxisLength));

plot(x1,LPS0_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS100_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS100_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS250_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS250_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS500_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS500_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x5,LPS750_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 5.
plot(x5,mean(LPS750_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS1000_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS1000_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 6.
plot(x7,LPS2000_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 7.
plot(x7,mean(LPS2000_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 7.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 15]);
xticks(1:15);
ylabel('Cell MajorAxisLength (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','0','','100','','250','','500','','750','','1000','','2000',''});

%Indicate significance on plot, if any.
sig=c_majorAxisLength(:,6)<0.05; %extract only significant results.
sig_MajorAxisLength=c_majorAxisLength(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_MajorAxisLength,1)
    if sig_MajorAxisLength(ii,6)<0.05
        if sig_MajorAxisLength(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_MajorAxisLength(ii,1)==2
            x1=4;
        elseif sig_MajorAxisLength(ii,1)==3
            x1=6;
        elseif sig_MajorAxisLength(ii,1)==4
            x1=8;
        elseif sig_MajorAxisLength(ii,1)==5
            x1=10;
        elseif sig_MajorAxisLength(ii,1)==6
            x1=12;
        else
            x1=14;
        end
        if sig_MajorAxisLength(ii,2)==2
            x2=4;
        elseif sig_MajorAxisLength(ii,2)==3
            x2=6;
        elseif sig_MajorAxisLength(ii,2)==4
            x2=8;
        elseif sig_MajorAxisLength(ii,2)==5
            x2=10;
        elseif sig_MajorAxisLength(ii,2)==6
            x2=12;
        else
            x2=14;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_MajorAxisLength(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_MajorAxisLength(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_MajorAxisLength(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'MajorAxisLength'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'MajorAxisLength.svg'); %save figure as an SVG file.
saveas(fig,'MajorAxisLength.jpg'); %save figure as an JPEG file.
close

%% PLOT MINOR AXIS LENGTH FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_MinorAxisLength));
x2=4*ones(1,length(LPS100_MinorAxisLength));
x3=6*ones(1,length(LPS250_MinorAxisLength));
x4=8*ones(1,length(LPS500_MinorAxisLength));
x5=10*ones(1,length(LPS750_MinorAxisLength));
x6=12*ones(1,length(LPS1000_MinorAxisLength));
x7=14*ones(1,length(LPS2000_MinorAxisLength));

plot(x1,LPS0_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS100_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS100_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS250_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS250_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS500_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS500_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x5,LPS750_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 5.
plot(x5,mean(LPS750_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS1000_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS1000_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 6.
plot(x7,LPS2000_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 7.
plot(x7,mean(LPS2000_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 7.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 15]);
xticks(1:15);
ylabel('Cell MinorAxisLength (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','0','','100','','250','','500','','750','','1000','','2000',''});

%Indicate significance on plot, if any.
sig=c_minorAxisLength(:,6)<0.05; %extract only significant results.
sig_MinorAxisLength=c_minorAxisLength(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_MinorAxisLength,1)
    if sig_MinorAxisLength(ii,6)<0.05
        if sig_MinorAxisLength(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_MinorAxisLength(ii,1)==2
            x1=4;
        elseif sig_MinorAxisLength(ii,1)==3
            x1=6;
        elseif sig_MinorAxisLength(ii,1)==4
            x1=8;
        elseif sig_MinorAxisLength(ii,1)==5
            x1=10;
        elseif sig_MinorAxisLength(ii,1)==6
            x1=12;
        else
            x1=14;
        end
        if sig_MinorAxisLength(ii,2)==2
            x2=4;
        elseif sig_MinorAxisLength(ii,2)==3
            x2=6;
        elseif sig_MinorAxisLength(ii,2)==4
            x2=8;
        elseif sig_MinorAxisLength(ii,2)==5
            x2=10;
        elseif sig_MinorAxisLength(ii,2)==6
            x2=12;
        else
            x2=14;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_MinorAxisLength(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_MinorAxisLength(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_MinorAxisLength(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'MinorAxisLength'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'MinorAxisLength.svg'); %save figure as an SVG file.
saveas(fig,'MinorAxisLength.jpg'); %save figure as an JPEG file.
close

%% PLOT ECCENTRICITY FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Eccentricity));
x2=4*ones(1,length(LPS100_Eccentricity));
x3=6*ones(1,length(LPS250_Eccentricity));
x4=8*ones(1,length(LPS500_Eccentricity));
x5=10*ones(1,length(LPS750_Eccentricity));
x6=12*ones(1,length(LPS1000_Eccentricity));
x7=14*ones(1,length(LPS2000_Eccentricity));

plot(x1,LPS0_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS100_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS100_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS250_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS250_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS500_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS500_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x5,LPS750_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 5.
plot(x5,mean(LPS750_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS1000_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS1000_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 6.
plot(x7,LPS2000_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 7.
plot(x7,mean(LPS2000_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 7.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 15]);
xticks(1:15);
ylabel('Cell Eccentricity (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','0','','100','','250','','500','','750','','1000','','2000',''});

%Indicate significance on plot, if any.
sig=c_eccentricity(:,6)<0.05; %extract only significant results.
sig_eccentricity=c_eccentricity(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_eccentricity,1)
    if sig_eccentricity(ii,6)<0.05
        if sig_eccentricity(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_eccentricity(ii,1)==2
            x1=4;
        elseif sig_eccentricity(ii,1)==3
            x1=6;
        elseif sig_eccentricity(ii,1)==4
            x1=8;
        elseif sig_eccentricity(ii,1)==5
            x1=10;
        elseif sig_eccentricity(ii,1)==6
            x1=12;
        else
            x1=14;
        end
        if sig_eccentricity(ii,2)==2
            x2=4;
        elseif sig_eccentricity(ii,2)==3
            x2=6;
        elseif sig_eccentricity(ii,2)==4
            x2=8;
        elseif sig_eccentricity(ii,2)==5
            x2=10;
        elseif sig_eccentricity(ii,2)==6
            x2=12;
        else
            x2=14;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_eccentricity(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_eccentricity(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_eccentricity(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'Eccentricity'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'Eccentricity.svg'); %save figure as an SVG file.
saveas(fig,'Eccentricity.jpg'); %save figure as an JPEG file.
close

%% PLOT CIRCULARITY FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Circularity));
x2=4*ones(1,length(LPS100_Circularity));
x3=6*ones(1,length(LPS250_Circularity));
x4=8*ones(1,length(LPS500_Circularity));
x5=10*ones(1,length(LPS750_Circularity));
x6=12*ones(1,length(LPS1000_Circularity));
x7=14*ones(1,length(LPS2000_Circularity));

plot(x1,LPS0_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Circularity),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS100_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS100_Circularity),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS250_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS250_Circularity),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS500_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS500_Circularity),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x5,LPS750_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 5.
plot(x5,mean(LPS750_Circularity),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS1000_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS1000_Circularity),'kx','MarkerSize',18); %plot mean for condition 6.
plot(x7,LPS2000_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 7.
plot(x7,mean(LPS2000_Circularity),'kx','MarkerSize',18); %plot mean for condition 7.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 15]);
xticks(1:15);
ylabel('Cell Circularity (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','0','','100','','250','','500','','750','','1000','','2000',''});

%Indicate significance on plot, if any.
sig=c_circularity(:,6)<0.05; %extract only significant results.
sig_circularity=c_circularity(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_circularity,1)
    if sig_circularity(ii,6)<0.05
        if sig_circularity(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_circularity(ii,1)==2
            x1=4;
        elseif sig_circularity(ii,1)==3
            x1=6;
        elseif sig_circularity(ii,1)==4
            x1=8;
        elseif sig_circularity(ii,1)==5
            x1=10;
        elseif sig_circularity(ii,1)==6
            x1=12;
        else
            x1=14;
        end
        if sig_circularity(ii,2)==2
            x2=4;
        elseif sig_circularity(ii,2)==3
            x2=6;
        elseif sig_circularity(ii,2)==4
            x2=8;
        elseif sig_circularity(ii,2)==5
            x2=10;
        elseif sig_circularity(ii,2)==6
            x2=12;
        else
            x2=14;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_circularity(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_circularity(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_circularity(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'Circularity'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'Circularity.svg'); %save figure as an SVG file.
saveas(fig,'Circularity.jpg'); %save figure as an JPEG file.
close

%% PLOT EXTENT FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Extent));
x2=4*ones(1,length(LPS100_Extent));
x3=6*ones(1,length(LPS250_Extent));
x4=8*ones(1,length(LPS500_Extent));
x5=10*ones(1,length(LPS750_Extent));
x6=12*ones(1,length(LPS1000_Extent));
x7=14*ones(1,length(LPS2000_Extent));

plot(x1,LPS0_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Extent),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS100_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS100_Extent),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS250_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS250_Extent),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS500_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS500_Extent),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x5,LPS750_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 5.
plot(x5,mean(LPS750_Extent),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS1000_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS1000_Extent),'kx','MarkerSize',18); %plot mean for condition 6.
plot(x7,LPS2000_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 7.
plot(x7,mean(LPS2000_Extent),'kx','MarkerSize',18); %plot mean for condition 7.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 15]);
xticks(1:15);
ylabel('Cell Extent (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','0','','100','','250','','500','','750','','1000','','2000',''});

%Indicate significance on plot, if any.
sig=c_extent(:,6)<0.05; %extract only significant results.
sig_extent=c_extent(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_extent,1)
    if sig_extent(ii,6)<0.05
        if sig_extent(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_extent(ii,1)==2
            x1=4;
        elseif sig_extent(ii,1)==3
            x1=6;
        elseif sig_extent(ii,1)==4
            x1=8;
        elseif sig_extent(ii,1)==5
            x1=10;
        elseif sig_extent(ii,1)==6
            x1=12;
        else
            x1=14;
        end
        if sig_extent(ii,2)==2
            x2=4;
        elseif sig_extent(ii,2)==3
            x2=6;
        elseif sig_extent(ii,2)==4
            x2=8;
        elseif sig_extent(ii,2)==5
            x2=10;
        elseif sig_extent(ii,2)==6
            x2=12;
        else
            x2=14;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_extent(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_extent(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_extent(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'Extent'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'Extent.svg'); %save figure as an SVG file.
saveas(fig,'Extent.jpg'); %save figure as an JPEG file.
close

%% PLOT PERIMETER FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Perimeter));
x2=4*ones(1,length(LPS100_Perimeter));
x3=6*ones(1,length(LPS250_Perimeter));
x4=8*ones(1,length(LPS500_Perimeter));
x5=10*ones(1,length(LPS750_Perimeter));
x6=12*ones(1,length(LPS1000_Perimeter));
x7=14*ones(1,length(LPS2000_Perimeter));

plot(x1,LPS0_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Perimeter),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS100_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS100_Perimeter),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS250_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS250_Perimeter),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS500_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS500_Perimeter),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x5,LPS750_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 5.
plot(x5,mean(LPS750_Perimeter),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS1000_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS1000_Perimeter),'kx','MarkerSize',18); %plot mean for condition 6.
plot(x7,LPS2000_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 7.
plot(x7,mean(LPS2000_Perimeter),'kx','MarkerSize',18); %plot mean for condition 7.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 15]);
xticks(1:15);
ylabel('Cell Perimeter (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','0','','100','','250','','500','','750','','1000','','2000',''});

%Indicate significance on plot, if any.
sig=c_perimeter(:,6)<0.05; %extract only significant results.
sig_perimeter=c_perimeter(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_perimeter,1)
    if sig_perimeter(ii,6)<0.05
        if sig_perimeter(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_perimeter(ii,1)==2
            x1=4;
        elseif sig_perimeter(ii,1)==3
            x1=6;
        elseif sig_perimeter(ii,1)==4
            x1=8;
        elseif sig_perimeter(ii,1)==5
            x1=10;
        elseif sig_perimeter(ii,1)==6
            x1=12;
        else
            x1=14;
        end
        if sig_perimeter(ii,2)==2
            x2=4;
        elseif sig_perimeter(ii,2)==3
            x2=6;
        elseif sig_perimeter(ii,2)==4
            x2=8;
        elseif sig_perimeter(ii,2)==5
            x2=10;
        elseif sig_perimeter(ii,2)==6
            x2=12;
        else
            x2=14;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_perimeter(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_perimeter(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_perimeter(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'Perimeter'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'Perimeter.svg'); %save figure as an SVG file.
saveas(fig,'Perimeter.jpg'); %save figure as an JPEG file.
close

%% PLOT ELONGATION FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Elongation));
x2=4*ones(1,length(LPS100_Elongation));
x3=6*ones(1,length(LPS250_Elongation));
x4=8*ones(1,length(LPS500_Elongation));
x5=10*ones(1,length(LPS750_Elongation));
x6=12*ones(1,length(LPS1000_Elongation));
x7=14*ones(1,length(LPS2000_Elongation));

plot(x1,LPS0_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Elongation),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x2,LPS100_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x2,mean(LPS100_Elongation),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x3,LPS250_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x3,mean(LPS250_Elongation),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x4,LPS500_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 4.
plot(x4,mean(LPS500_Elongation),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x5,LPS750_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 5.
plot(x5,mean(LPS750_Elongation),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS1000_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS1000_Elongation),'kx','MarkerSize',18); %plot mean for condition 6.
plot(x7,LPS2000_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 7.
plot(x7,mean(LPS2000_Elongation),'kx','MarkerSize',18); %plot mean for condition 7.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([1 15]);
xticks(1:15);
ylabel('Cell Elongation (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','0','','100','','250','','500','','750','','1000','','2000',''});

%Indicate significance on plot, if any.
sig=c_elongation(:,6)<0.05; %extract only significant results.
sig_elongation=c_elongation(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_elongation,1)
    if sig_elongation(ii,6)<0.05
        if sig_elongation(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_elongation(ii,1)==2
            x1=4;
        elseif sig_elongation(ii,1)==3
            x1=6;
        elseif sig_elongation(ii,1)==4
            x1=8;
        elseif sig_elongation(ii,1)==5
            x1=10;
        elseif sig_elongation(ii,1)==6
            x1=12;
        else
            x1=14;
        end
        if sig_elongation(ii,2)==2
            x2=4;
        elseif sig_elongation(ii,2)==3
            x2=6;
        elseif sig_elongation(ii,2)==4
            x2=8;
        elseif sig_elongation(ii,2)==5
            x2=10;
        elseif sig_elongation(ii,2)==6
            x2=12;
        else
            x2=14;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_elongation(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_elongation(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_elongation(ii,6)<0.01
            x=[x(1)+0.25,x(1)+0.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        else
            x=x(1)+0.25; %x-coordinate of asterisk.
            
            %Find y-coordinate.
            yrange=a+a;
            y8=yrange/8;
            y=max(ylim)+y8; %y-coordinate of asterisk.
            
            plot(x,y,'k*'); %plot asterisk.
        end
    end
    ylim([0,yt(end)+ydiff]);
end

%Save plot.
savefig(fig,'Elongation'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'Elongation.svg'); %save figure as an SVG file.
saveas(fig,'Elongation.jpg'); %save figure as an JPEG file.
close
end