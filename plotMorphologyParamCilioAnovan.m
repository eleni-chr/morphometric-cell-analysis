function plotMorphologyParamCilioAnovan(data)
%% PREPARE DATA
%load('dataAnalysis.mat','data');
%data2=[data2([1:36],:);data2([38:57],:);data2([59:end],:)];
%data3=[data([1:7],:);data([9:60],:);data([62:71],:);data([73:end],:)];

%SEPARATE DATA ACCORDING TO CONDITION
%Initialise variables.
LPS0_Cilio0=[];
LPS0_Cilio20=[];
LPS0_Cilio50=[];
LPS10_Cilio0=[];
LPS10_Cilio20=[];
LPS10_Cilio50=[];

for ii=1:length(data)
    if ~isempty(strfind(data{ii,1},'_LPS+DMSO'))
        LPS10_Cilio0=[LPS10_Cilio0;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'_LPS+20uMcilio'))
        LPS10_Cilio20=[LPS10_Cilio20;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'_LPS+50uMcilio'))
        LPS10_Cilio50=[LPS10_Cilio50;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'_HBSS+DMSO'))
        LPS0_Cilio0=[LPS0_Cilio0;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'_HBSS+20uMcilio'))
        LPS0_Cilio20=[LPS0_Cilio20;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'_HBSS+50uMcilio'))
        LPS0_Cilio50=[LPS0_Cilio50;data{ii,3}];
    end
end

%EXTRACT PARAMETERS FOR EACH CONDITION
%Initialise variables.
LPS0_Cilio0_Area=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_Area=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_Area=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_Area=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_Area=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_Area=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_MajorAxisLength=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_MajorAxisLength=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_MajorAxisLength=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_MajorAxisLength=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_MajorAxisLength=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_MajorAxisLength=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_MinorAxisLength=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_MinorAxisLength=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_MinorAxisLength=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_MinorAxisLength=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_MinorAxisLength=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_MinorAxisLength=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_Eccentricity=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_Eccentricity=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_Eccentricity=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_Eccentricity=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_Eccentricity=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_Eccentricity=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_Circularity=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_Circularity=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_Circularity=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_Circularity=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_Circularity=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_Circularity=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_Extent=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_Extent=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_Extent=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_Extent=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_Extent=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_Extent=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_Perimeter=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_Perimeter=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_Perimeter=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_Perimeter=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_Perimeter=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_Perimeter=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_Elongation=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_Elongation=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_Elongation=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_Elongation=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_Elongation=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_Elongation=zeros(length(LPS10_Cilio50),1);

%CONDITION 1 PARAMETERS
for ii=1:length(LPS0_Cilio0)
    LPS0_Cilio0_Area(ii,1)=LPS0_Cilio0(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS0_Cilio0_MajorAxisLength(ii,1)=LPS0_Cilio0(ii).MajorAxisLength*0.18;
    LPS0_Cilio0_MinorAxisLength(ii,1)=LPS0_Cilio0(ii).MinorAxisLength*0.18;
    LPS0_Cilio0_Eccentricity(ii,1)=LPS0_Cilio0(ii).Eccentricity;
    LPS0_Cilio0_Circularity(ii,1)=LPS0_Cilio0(ii).Circularity;
    LPS0_Cilio0_Extent(ii,1)=LPS0_Cilio0(ii).Extent;
    LPS0_Cilio0_Perimeter(ii,1)=LPS0_Cilio0(ii).Perimeter*0.18;
    LPS0_Cilio0_Elongation(ii,1)=(LPS0_Cilio0(ii).MajorAxisLength*0.18)/(LPS0_Cilio0(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 2 PARAMETERS
for ii=1:length(LPS0_Cilio20)
    LPS0_Cilio20_Area(ii,1)=LPS0_Cilio20(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS0_Cilio20_MajorAxisLength(ii,1)=LPS0_Cilio20(ii).MajorAxisLength*0.18;
    LPS0_Cilio20_MinorAxisLength(ii,1)=LPS0_Cilio20(ii).MinorAxisLength*0.18;
    LPS0_Cilio20_Eccentricity(ii,1)=LPS0_Cilio20(ii).Eccentricity;
    LPS0_Cilio20_Circularity(ii,1)=LPS0_Cilio20(ii).Circularity;
    LPS0_Cilio20_Extent(ii,1)=LPS0_Cilio20(ii).Extent;
    LPS0_Cilio20_Perimeter(ii,1)=LPS0_Cilio20(ii).Perimeter*0.18;
    LPS0_Cilio20_Elongation(ii,1)=(LPS0_Cilio20(ii).MajorAxisLength*0.18)/(LPS0_Cilio20(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 3 PARAMETERS
for ii=1:length(LPS0_Cilio50)
    LPS0_Cilio50_Area(ii,1)=LPS0_Cilio50(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS0_Cilio50_MajorAxisLength(ii,1)=LPS0_Cilio50(ii).MajorAxisLength*0.18;
    LPS0_Cilio50_MinorAxisLength(ii,1)=LPS0_Cilio50(ii).MinorAxisLength*0.18;
    LPS0_Cilio50_Eccentricity(ii,1)=LPS0_Cilio50(ii).Eccentricity;
    LPS0_Cilio50_Circularity(ii,1)=LPS0_Cilio50(ii).Circularity;
    LPS0_Cilio50_Extent(ii,1)=LPS0_Cilio50(ii).Extent;
    LPS0_Cilio50_Perimeter(ii,1)=LPS0_Cilio50(ii).Perimeter*0.18;
    LPS0_Cilio50_Elongation(ii,1)=(LPS0_Cilio50(ii).MajorAxisLength*0.18)/(LPS0_Cilio50(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 4 PARAMETERS
for ii=1:length(LPS10_Cilio0)
    LPS10_Cilio0_Area(ii,1)=LPS10_Cilio0(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS10_Cilio0_MajorAxisLength(ii,1)=LPS10_Cilio0(ii).MajorAxisLength*0.18;
    LPS10_Cilio0_MinorAxisLength(ii,1)=LPS10_Cilio0(ii).MinorAxisLength*0.18;
    LPS10_Cilio0_Eccentricity(ii,1)=LPS10_Cilio0(ii).Eccentricity;
    LPS10_Cilio0_Circularity(ii,1)=LPS10_Cilio0(ii).Circularity;
    LPS10_Cilio0_Extent(ii,1)=LPS10_Cilio0(ii).Extent;
    LPS10_Cilio0_Perimeter(ii,1)=LPS10_Cilio0(ii).Perimeter*0.18;
    LPS10_Cilio0_Elongation(ii,1)=(LPS10_Cilio0(ii).MajorAxisLength*0.18)/(LPS10_Cilio0(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 5 PARAMETERS
for ii=1:length(LPS10_Cilio20)
    LPS10_Cilio20_Area(ii,1)=LPS10_Cilio20(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS10_Cilio20_MajorAxisLength(ii,1)=LPS10_Cilio20(ii).MajorAxisLength*0.18;
    LPS10_Cilio20_MinorAxisLength(ii,1)=LPS10_Cilio20(ii).MinorAxisLength*0.18;
    LPS10_Cilio20_Eccentricity(ii,1)=LPS10_Cilio20(ii).Eccentricity;
    LPS10_Cilio20_Circularity(ii,1)=LPS10_Cilio20(ii).Circularity;
    LPS10_Cilio20_Extent(ii,1)=LPS10_Cilio20(ii).Extent;
    LPS10_Cilio20_Perimeter(ii,1)=LPS10_Cilio20(ii).Perimeter*0.18;
    LPS10_Cilio20_Elongation(ii,1)=(LPS10_Cilio20(ii).MajorAxisLength*0.18)/(LPS10_Cilio20(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%CONDITION 6 PARAMETERS
for ii=1:length(LPS10_Cilio50)
    LPS10_Cilio50_Area(ii,1)=LPS10_Cilio50(ii).Area*0.18; %multiply by pixel size to convert pixels to microns.
    LPS10_Cilio50_MajorAxisLength(ii,1)=LPS10_Cilio50(ii).MajorAxisLength*0.18;
    LPS10_Cilio50_MinorAxisLength(ii,1)=LPS10_Cilio50(ii).MinorAxisLength*0.18;
    LPS10_Cilio50_Eccentricity(ii,1)=LPS10_Cilio50(ii).Eccentricity;
    LPS10_Cilio50_Circularity(ii,1)=LPS10_Cilio50(ii).Circularity;
    LPS10_Cilio50_Extent(ii,1)=LPS10_Cilio50(ii).Extent;
    LPS10_Cilio50_Perimeter(ii,1)=LPS10_Cilio50(ii).Perimeter*0.18;
    LPS10_Cilio50_Elongation(ii,1)=(LPS10_Cilio50(ii).MajorAxisLength*0.18)/(LPS10_Cilio50(ii).MinorAxisLength*0.18); %calculation additional parameter.
end

%% REMOVE OUTLIERS

%CONDITION 1
for ii=1:length(LPS0_Cilio0)
    LPS0_Cilio0_Area=rmoutliers(LPS0_Cilio0_Area);
    LPS0_Cilio0_MajorAxisLength=rmoutliers(LPS0_Cilio0_MajorAxisLength);
    LPS0_Cilio0_MinorAxisLength=rmoutliers(LPS0_Cilio0_MinorAxisLength);
    LPS0_Cilio0_Eccentricity=rmoutliers(LPS0_Cilio0_Eccentricity);
    LPS0_Cilio0_Circularity=rmoutliers(LPS0_Cilio0_Circularity);
    LPS0_Cilio0_Extent=rmoutliers(LPS0_Cilio0_Extent);
    LPS0_Cilio0_Perimeter=rmoutliers(LPS0_Cilio0_Perimeter);
    LPS0_Cilio0_Elongation=rmoutliers(LPS0_Cilio0_Elongation);
end

%CONDITION 2
for ii=1:length(LPS0_Cilio20)
    LPS0_Cilio20_Area=rmoutliers(LPS0_Cilio20_Area);
    LPS0_Cilio20_MajorAxisLength=rmoutliers(LPS0_Cilio20_MajorAxisLength);
    LPS0_Cilio20_MinorAxisLength=rmoutliers(LPS0_Cilio20_MinorAxisLength);
    LPS0_Cilio20_Eccentricity=rmoutliers(LPS0_Cilio20_Eccentricity);
    LPS0_Cilio20_Circularity=rmoutliers(LPS0_Cilio20_Circularity);
    LPS0_Cilio20_Extent=rmoutliers(LPS0_Cilio20_Extent);
    LPS0_Cilio20_Perimeter=rmoutliers(LPS0_Cilio20_Perimeter);
    LPS0_Cilio20_Elongation=rmoutliers(LPS0_Cilio20_Elongation);
end

%CONDITION 3
for ii=1:length(LPS0_Cilio50)
    LPS0_Cilio50_Area=rmoutliers(LPS0_Cilio50_Area);
    LPS0_Cilio50_MajorAxisLength=rmoutliers(LPS0_Cilio50_MajorAxisLength);
    LPS0_Cilio50_MinorAxisLength=rmoutliers(LPS0_Cilio50_MinorAxisLength);
    LPS0_Cilio50_Eccentricity=rmoutliers(LPS0_Cilio50_Eccentricity);
    LPS0_Cilio50_Circularity=rmoutliers(LPS0_Cilio50_Circularity);
    LPS0_Cilio50_Extent=rmoutliers(LPS0_Cilio50_Extent);
    LPS0_Cilio50_Perimeter=rmoutliers(LPS0_Cilio50_Perimeter);
    LPS0_Cilio50_Elongation=rmoutliers(LPS0_Cilio50_Elongation);
end

%CONDITION 4
for ii=1:length(LPS10_Cilio0)
    LPS10_Cilio0_Area=rmoutliers(LPS10_Cilio0_Area);
    LPS10_Cilio0_MajorAxisLength=rmoutliers(LPS10_Cilio0_MajorAxisLength);
    LPS10_Cilio0_MinorAxisLength=rmoutliers(LPS10_Cilio0_MinorAxisLength);
    LPS10_Cilio0_Eccentricity=rmoutliers(LPS10_Cilio0_Eccentricity);
    LPS10_Cilio0_Circularity=rmoutliers(LPS10_Cilio0_Circularity);
    LPS10_Cilio0_Extent=rmoutliers(LPS10_Cilio0_Extent);
    LPS10_Cilio0_Perimeter=rmoutliers(LPS10_Cilio0_Perimeter);
    LPS10_Cilio0_Elongation=rmoutliers(LPS10_Cilio0_Elongation);
end

%CONDITION 5
for ii=1:length(LPS10_Cilio20)
    LPS10_Cilio20_Area=rmoutliers(LPS10_Cilio20_Area);
    LPS10_Cilio20_MajorAxisLength=rmoutliers(LPS10_Cilio20_MajorAxisLength);
    LPS10_Cilio20_MinorAxisLength=rmoutliers(LPS10_Cilio20_MinorAxisLength);
    LPS10_Cilio20_Eccentricity=rmoutliers(LPS10_Cilio20_Eccentricity);
    LPS10_Cilio20_Circularity=rmoutliers(LPS10_Cilio20_Circularity);
    LPS10_Cilio20_Extent=rmoutliers(LPS10_Cilio20_Extent);
    LPS10_Cilio20_Perimeter=rmoutliers(LPS10_Cilio20_Perimeter);
    LPS10_Cilio20_Elongation=rmoutliers(LPS10_Cilio20_Elongation);
end

%CONDITION 6
for ii=1:length(LPS10_Cilio50)
    LPS10_Cilio50_Area=rmoutliers(LPS10_Cilio50_Area);
    LPS10_Cilio50_MajorAxisLength=rmoutliers(LPS10_Cilio50_MajorAxisLength);
    LPS10_Cilio50_MinorAxisLength=rmoutliers(LPS10_Cilio50_MinorAxisLength);
    LPS10_Cilio50_Eccentricity=rmoutliers(LPS10_Cilio50_Eccentricity);
    LPS10_Cilio50_Circularity=rmoutliers(LPS10_Cilio50_Circularity);
    LPS10_Cilio50_Extent=rmoutliers(LPS10_Cilio50_Extent);
    LPS10_Cilio50_Perimeter=rmoutliers(LPS10_Cilio50_Perimeter);
    LPS10_Cilio50_Elongation=rmoutliers(LPS10_Cilio50_Elongation);
end

%% MERGE DATA FOR STATISTICAL ANALYSIS

%Merge data in order.
area=[LPS0_Cilio0_Area;LPS0_Cilio20_Area;LPS0_Cilio50_Area;LPS10_Cilio0_Area;LPS10_Cilio20_Area;LPS10_Cilio50_Area];
majorAxisLength=[LPS0_Cilio0_MajorAxisLength;LPS0_Cilio20_MajorAxisLength;LPS0_Cilio50_MajorAxisLength;LPS10_Cilio0_MajorAxisLength;LPS10_Cilio20_MajorAxisLength;LPS10_Cilio50_MajorAxisLength];
minorAxisLength=[LPS0_Cilio0_MinorAxisLength;LPS0_Cilio20_MinorAxisLength;LPS0_Cilio50_MinorAxisLength;LPS10_Cilio0_MinorAxisLength;LPS10_Cilio20_MinorAxisLength;LPS10_Cilio50_MinorAxisLength];
eccentricity=[LPS0_Cilio0_Eccentricity;LPS0_Cilio20_Eccentricity;LPS0_Cilio50_Eccentricity;LPS10_Cilio0_Eccentricity;LPS10_Cilio20_Eccentricity;LPS10_Cilio50_Eccentricity];
circularity=[LPS0_Cilio0_Circularity;LPS0_Cilio20_Circularity;LPS0_Cilio50_Circularity;LPS10_Cilio0_Circularity;LPS10_Cilio20_Circularity;LPS10_Cilio50_Circularity];
extent=[LPS0_Cilio0_Extent;LPS0_Cilio20_Extent;LPS0_Cilio50_Extent;LPS10_Cilio0_Extent;LPS10_Cilio20_Extent;LPS10_Cilio50_Extent];
perimeter=[LPS0_Cilio0_Perimeter;LPS0_Cilio20_Perimeter;LPS0_Cilio50_Perimeter;LPS10_Cilio0_Perimeter;LPS10_Cilio20_Perimeter;LPS10_Cilio50_Perimeter];
elongation=[LPS0_Cilio0_Elongation;LPS0_Cilio20_Elongation;LPS0_Cilio50_Elongation;LPS10_Cilio0_Elongation;LPS10_Cilio20_Elongation;LPS10_Cilio50_Elongation];

%Grouping variables.
lenAreaLPS0=length(LPS0_Cilio0_Area)+length(LPS0_Cilio20_Area)+length(LPS0_Cilio50_Area);
lenAreaLPS10=length(LPS10_Cilio0_Area)+length(LPS10_Cilio20_Area)+length(LPS10_Cilio50_Area);
LPS0_Area=ones(1,lenAreaLPS0);
LPS10_Area=2*ones(1,lenAreaLPS10);
LPS_Area=[LPS0_Area,LPS10_Area];

lenMajorAxisLengthLPS0=length(LPS0_Cilio0_MajorAxisLength)+length(LPS0_Cilio20_MajorAxisLength)+length(LPS0_Cilio50_MajorAxisLength);
lenMajorAxisLengthLPS10=length(LPS10_Cilio0_MajorAxisLength)+length(LPS10_Cilio20_MajorAxisLength)+length(LPS10_Cilio50_MajorAxisLength);
LPS0_MajorAxisLength=ones(1,lenMajorAxisLengthLPS0);
LPS10_MajorAxisLength=2*ones(1,lenMajorAxisLengthLPS10);
LPS_MajorAxisLength=[LPS0_MajorAxisLength,LPS10_MajorAxisLength];

lenMinorAxisLengthLPS0=length(LPS0_Cilio0_MinorAxisLength)+length(LPS0_Cilio20_MinorAxisLength)+length(LPS0_Cilio50_MinorAxisLength);
lenMinorAxisLengthLPS10=length(LPS10_Cilio0_MinorAxisLength)+length(LPS10_Cilio20_MinorAxisLength)+length(LPS10_Cilio50_MinorAxisLength);
LPS0_MinorAxisLength=ones(1,lenMinorAxisLengthLPS0);
LPS10_MinorAxisLength=2*ones(1,lenMinorAxisLengthLPS10);
LPS_MinorAxisLength=[LPS0_MinorAxisLength,LPS10_MinorAxisLength];

lenEccentricityLPS0=length(LPS0_Cilio0_Eccentricity)+length(LPS0_Cilio20_Eccentricity)+length(LPS0_Cilio50_Eccentricity);
lenEccentricityLPS10=length(LPS10_Cilio0_Eccentricity)+length(LPS10_Cilio20_Eccentricity)+length(LPS10_Cilio50_Eccentricity);
LPS0_Eccentricity=ones(1,lenEccentricityLPS0);
LPS10_Eccentricity=2*ones(1,lenEccentricityLPS10);
LPS_Eccentricity=[LPS0_Eccentricity,LPS10_Eccentricity];

lenCircularityLPS0=length(LPS0_Cilio0_Circularity)+length(LPS0_Cilio20_Circularity)+length(LPS0_Cilio50_Circularity);
lenCircularityLPS10=length(LPS10_Cilio0_Circularity)+length(LPS10_Cilio20_Circularity)+length(LPS10_Cilio50_Circularity);
LPS0_Circularity=ones(1,lenCircularityLPS0);
LPS10_Circularity=2*ones(1,lenCircularityLPS10);
LPS_Circularity=[LPS0_Circularity,LPS10_Circularity];

lenExtentLPS0=length(LPS0_Cilio0_Extent)+length(LPS0_Cilio20_Extent)+length(LPS0_Cilio50_Extent);
lenExtentLPS10=length(LPS10_Cilio0_Extent)+length(LPS10_Cilio20_Extent)+length(LPS10_Cilio50_Extent);
LPS0_Extent=ones(1,lenExtentLPS0);
LPS10_Extent=2*ones(1,lenExtentLPS10);
LPS_Extent=[LPS0_Extent,LPS10_Extent];

lenPerimeterLPS0=length(LPS0_Cilio0_Perimeter)+length(LPS0_Cilio20_Perimeter)+length(LPS0_Cilio50_Perimeter);
lenPerimeterLPS10=length(LPS10_Cilio0_Perimeter)+length(LPS10_Cilio20_Perimeter)+length(LPS10_Cilio50_Perimeter);
LPS0_Perimeter=ones(1,lenPerimeterLPS0);
LPS10_Perimeter=2*ones(1,lenPerimeterLPS10);
LPS_Perimeter=[LPS0_Perimeter,LPS10_Perimeter];

lenElongationLPS0=length(LPS0_Cilio0_Elongation)+length(LPS0_Cilio20_Elongation)+length(LPS0_Cilio50_Elongation);
lenElongationLPS10=length(LPS10_Cilio0_Elongation)+length(LPS10_Cilio20_Elongation)+length(LPS10_Cilio50_Elongation);
LPS0_Elongation=ones(1,lenElongationLPS0);
LPS10_Elongation=2*ones(1,lenElongationLPS10);
LPS_Elongation=[LPS0_Elongation,LPS10_Elongation];

lenAreaCilio0=length(LPS0_Cilio0_Area)+length(LPS10_Cilio0_Area);
lenAreaCilio20=length(LPS0_Cilio20_Area)+length(LPS10_Cilio20_Area);
lenAreaCilio50=length(LPS0_Cilio50_Area)+length(LPS10_Cilio50_Area);
Cilio0_Area=ones(1,lenAreaCilio0);
Cilio20_Area=2*ones(1,lenAreaCilio20);
Cilio50_Area=3*ones(1,lenAreaCilio50);
Cilio_Area=[Cilio0_Area,Cilio20_Area,Cilio50_Area];

lenMajorAxisLengthCilio0=length(LPS0_Cilio0_MajorAxisLength)+length(LPS10_Cilio0_MajorAxisLength);
lenMajorAxisLengthCilio20=length(LPS0_Cilio20_MajorAxisLength)+length(LPS10_Cilio20_MajorAxisLength);
lenMajorAxisLengthCilio50=length(LPS0_Cilio50_MajorAxisLength)+length(LPS10_Cilio50_MajorAxisLength);
Cilio0_MajorAxisLength=ones(1,lenMajorAxisLengthCilio0);
Cilio20_MajorAxisLength=2*ones(1,lenMajorAxisLengthCilio20);
Cilio50_MajorAxisLength=3*ones(1,lenMajorAxisLengthCilio50);
Cilio_MajorAxisLength=[Cilio0_MajorAxisLength,Cilio20_MajorAxisLength,Cilio50_MajorAxisLength];

lenMinorAxisLengthCilio0=length(LPS0_Cilio0_MinorAxisLength)+length(LPS10_Cilio0_MinorAxisLength);
lenMinorAxisLengthCilio20=length(LPS0_Cilio20_MinorAxisLength)+length(LPS10_Cilio20_MinorAxisLength);
lenMinorAxisLengthCilio50=length(LPS0_Cilio50_MinorAxisLength)+length(LPS10_Cilio50_MinorAxisLength);
Cilio0_MinorAxisLength=ones(1,lenMinorAxisLengthCilio0);
Cilio20_MinorAxisLength=2*ones(1,lenMinorAxisLengthCilio20);
Cilio50_MinorAxisLength=3*ones(1,lenMinorAxisLengthCilio50);
Cilio_MinorAxisLength=[Cilio0_MinorAxisLength,Cilio20_MinorAxisLength,Cilio50_MinorAxisLength];

lenEccentricityCilio0=length(LPS0_Cilio0_Eccentricity)+length(LPS10_Cilio0_Eccentricity);
lenEccentricityCilio20=length(LPS0_Cilio20_Eccentricity)+length(LPS10_Cilio20_Eccentricity);
lenEccentricityCilio50=length(LPS0_Cilio50_Eccentricity)+length(LPS10_Cilio50_Eccentricity);
Cilio0_Eccentricity=ones(1,lenEccentricityCilio0);
Cilio20_Eccentricity=2*ones(1,lenEccentricityCilio20);
Cilio50_Eccentricity=3*ones(1,lenEccentricityCilio50);
Cilio_Eccentricity=[Cilio0_Eccentricity,Cilio20_Eccentricity,Cilio50_Eccentricity];

lenCircularityCilio0=length(LPS0_Cilio0_Circularity)+length(LPS10_Cilio0_Circularity);
lenCircularityCilio20=length(LPS0_Cilio20_Circularity)+length(LPS10_Cilio20_Circularity);
lenCircularityCilio50=length(LPS0_Cilio50_Circularity)+length(LPS10_Cilio50_Circularity);
Cilio0_Circularity=ones(1,lenCircularityCilio0);
Cilio20_Circularity=2*ones(1,lenCircularityCilio20);
Cilio50_Circularity=3*ones(1,lenCircularityCilio50);
Cilio_Circularity=[Cilio0_Circularity,Cilio20_Circularity,Cilio50_Circularity];

lenExtentCilio0=length(LPS0_Cilio0_Extent)+length(LPS10_Cilio0_Extent);
lenExtentCilio20=length(LPS0_Cilio20_Extent)+length(LPS10_Cilio20_Extent);
lenExtentCilio50=length(LPS0_Cilio50_Extent)+length(LPS10_Cilio50_Extent);
Cilio0_Extent=ones(1,lenExtentCilio0);
Cilio20_Extent=2*ones(1,lenExtentCilio20);
Cilio50_Extent=3*ones(1,lenExtentCilio50);
Cilio_Extent=[Cilio0_Extent,Cilio20_Extent,Cilio50_Extent];

lenPerimeterCilio0=length(LPS0_Cilio0_Perimeter)+length(LPS10_Cilio0_Perimeter);
lenPerimeterCilio20=length(LPS0_Cilio20_Perimeter)+length(LPS10_Cilio20_Perimeter);
lenPerimeterCilio50=length(LPS0_Cilio50_Perimeter)+length(LPS10_Cilio50_Perimeter);
Cilio0_Perimeter=ones(1,lenPerimeterCilio0);
Cilio20_Perimeter=2*ones(1,lenPerimeterCilio20);
Cilio50_Perimeter=3*ones(1,lenPerimeterCilio50);
Cilio_Perimeter=[Cilio0_Perimeter,Cilio20_Perimeter,Cilio50_Perimeter];

lenElongationCilio0=length(LPS0_Cilio0_Elongation)+length(LPS10_Cilio0_Elongation);
lenElongationCilio20=length(LPS0_Cilio20_Elongation)+length(LPS10_Cilio20_Elongation);
lenElongationCilio50=length(LPS0_Cilio50_Elongation)+length(LPS10_Cilio50_Elongation);
Cilio0_Elongation=ones(1,lenElongationCilio0);
Cilio20_Elongation=2*ones(1,lenElongationCilio20);
Cilio50_Elongation=3*ones(1,lenElongationCilio50);
Cilio_Elongation=[Cilio0_Elongation,Cilio20_Elongation,Cilio50_Elongation];

%% PERFORM STATISTICS FOR EACH PARAMETER

%Perform two-way ANOVAs.
[p_area,~,stats_area]=anovan(area,{LPS_Area,Cilio_Area});
[p_majorAxisLength,~,stats_majorAxisLength]=anovan(majorAxisLength,{LPS_MajorAxisLength,Cilio_MajorAxisLength});
[p_minorAxisLength,~,stats_minorAxisLength]=anovan(minorAxisLength,{LPS_MinorAxisLength,Cilio_MinorAxisLength});
[p_eccentricity,~,stats_eccentricity]=anovan(eccentricity,{LPS_Eccentricity,Cilio_Eccentricity});
[p_circularity,~,stats_circularity]=anovan(circularity,{LPS_Circularity,Cilio_Circularity});
[p_extent,~,stats_extent]=anovan(extent,{LPS_Extent,Cilio_Extent});
[p_perimeter,~,stats_perimeter]=anovan(perimeter,{LPS_Perimeter,Cilio_Perimeter});
[p_elongation,~,stats_elongation]=anovan(elongation,{LPS_Elongation,Cilio_Elongation});
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

save('AnovanResults.mat','p_area','p_majorAxisLength','p_minorAxisLength','p_eccentricity','p_circularity','p_extent','p_perimeter','p_elongation');
save('multCompResultsAnovan.mat','c_area','c_majorAxisLength','c_minorAxisLength','c_eccentricity','c_circularity','c_extent','c_perimeter','c_elongation');

%% PLOT AREA FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Cilio0_Area));
x3=8*ones(1,length(LPS0_Cilio20_Area));
x5=14*ones(1,length(LPS0_Cilio50_Area));
x2=4*ones(1,length(LPS10_Cilio0_Area));
x4=10*ones(1,length(LPS10_Cilio20_Area));
x6=16*ones(1,length(LPS10_Cilio50_Area));

plot(x1,LPS0_Cilio0_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_Area),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_Area),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_Area,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_Area),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_Area,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_Area),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_Area,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_Area),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_Area,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_Area),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Cell area (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

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
            x1=5;
        elseif sig_area(ii,1)==3
            x1=8;
        elseif sig_area(ii,1)==4
            x1=11;
        elseif sig_area(ii,1)==5
            x1=14;
        else
            x1=17;
        end
        if sig_area(ii,2)==2
            x2=5;
        elseif sig_area(ii,2)==3
            x2=8;
        elseif sig_area(ii,2)==4
            x2=11;
        elseif sig_area(ii,2)==5
            x2=14;
        else
            x2=17;
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

x1=2*ones(1,length(LPS0_Cilio0_MajorAxisLength));
x3=8*ones(1,length(LPS0_Cilio20_MajorAxisLength));
x5=14*ones(1,length(LPS0_Cilio50_MajorAxisLength));
x2=4*ones(1,length(LPS10_Cilio0_MajorAxisLength));
x4=10*ones(1,length(LPS10_Cilio20_MajorAxisLength));
x6=16*ones(1,length(LPS10_Cilio50_MajorAxisLength));

plot(x1,LPS0_Cilio0_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_MajorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_MajorAxisLength,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_MajorAxisLength,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_MajorAxisLength,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_MajorAxisLength),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Cell MajorAxisLength (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

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
            x1=5;
        elseif sig_MajorAxisLength(ii,1)==3
            x1=8;
        elseif sig_MajorAxisLength(ii,1)==4
            x1=11;
        elseif sig_MajorAxisLength(ii,1)==5
            x1=14;
        else
            x1=17;
        end
        if sig_MajorAxisLength(ii,2)==2
            x2=5;
        elseif sig_MajorAxisLength(ii,2)==3
            x2=8;
        elseif sig_MajorAxisLength(ii,2)==4
            x2=11;
        elseif sig_MajorAxisLength(ii,2)==5
            x2=14;
        else
            x2=17;
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

x1=2*ones(1,length(LPS0_Cilio0_MinorAxisLength));
x3=8*ones(1,length(LPS0_Cilio20_MinorAxisLength));
x5=14*ones(1,length(LPS0_Cilio50_MinorAxisLength));
x2=4*ones(1,length(LPS10_Cilio0_MinorAxisLength));
x4=10*ones(1,length(LPS10_Cilio20_MinorAxisLength));
x6=16*ones(1,length(LPS10_Cilio50_MinorAxisLength));

plot(x1,LPS0_Cilio0_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_MinorAxisLength,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_MinorAxisLength,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_MinorAxisLength,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_MinorAxisLength,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_MinorAxisLength),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Cell MinorAxisLength (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

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
            x1=5;
        elseif sig_MinorAxisLength(ii,1)==3
            x1=8;
        elseif sig_MinorAxisLength(ii,1)==4
            x1=11;
        elseif sig_MinorAxisLength(ii,1)==5
            x1=14;
        else
            x1=17;
        end
        if sig_MinorAxisLength(ii,2)==2
            x2=5;
        elseif sig_MinorAxisLength(ii,2)==3
            x2=8;
        elseif sig_MinorAxisLength(ii,2)==4
            x2=11;
        elseif sig_MinorAxisLength(ii,2)==5
            x2=14;
        else
            x2=17;
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

x1=2*ones(1,length(LPS0_Cilio0_Eccentricity));
x3=8*ones(1,length(LPS0_Cilio20_Eccentricity));
x5=14*ones(1,length(LPS0_Cilio50_Eccentricity));
x2=4*ones(1,length(LPS10_Cilio0_Eccentricity));
x4=10*ones(1,length(LPS10_Cilio20_Eccentricity));
x6=16*ones(1,length(LPS10_Cilio50_Eccentricity));

plot(x1,LPS0_Cilio0_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_Eccentricity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_Eccentricity,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_Eccentricity,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_Eccentricity,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_Eccentricity),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Cell Eccentricity (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

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
            x1=5;
        elseif sig_eccentricity(ii,1)==3
            x1=8;
        elseif sig_eccentricity(ii,1)==4
            x1=11;
        elseif sig_eccentricity(ii,1)==5
            x1=14;
        else
            x1=17;
        end
        if sig_eccentricity(ii,2)==2
            x2=5;
        elseif sig_eccentricity(ii,2)==3
            x2=8;
        elseif sig_eccentricity(ii,2)==4
            x2=11;
        elseif sig_eccentricity(ii,2)==5
            x2=14;
        else
            x2=17;
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

x1=2*ones(1,length(LPS0_Cilio0_Circularity));
x3=8*ones(1,length(LPS0_Cilio20_Circularity));
x5=14*ones(1,length(LPS0_Cilio50_Circularity));
x2=4*ones(1,length(LPS10_Cilio0_Circularity));
x4=10*ones(1,length(LPS10_Cilio20_Circularity));
x6=16*ones(1,length(LPS10_Cilio50_Circularity));

plot(x1,LPS0_Cilio0_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_Circularity),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_Circularity),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_Circularity,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_Circularity),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_Circularity,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_Circularity),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_Circularity,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_Circularity),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_Circularity,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_Circularity),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Cell Circularity (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

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
            x1=5;
        elseif sig_circularity(ii,1)==3
            x1=8;
        elseif sig_circularity(ii,1)==4
            x1=11;
        elseif sig_circularity(ii,1)==5
            x1=14;
        else
            x1=17;
        end
        if sig_circularity(ii,2)==2
            x2=5;
        elseif sig_circularity(ii,2)==3
            x2=8;
        elseif sig_circularity(ii,2)==4
            x2=11;
        elseif sig_circularity(ii,2)==5
            x2=14;
        else
            x2=17;
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

x1=2*ones(1,length(LPS0_Cilio0_Extent));
x3=8*ones(1,length(LPS0_Cilio20_Extent));
x5=14*ones(1,length(LPS0_Cilio50_Extent));
x2=4*ones(1,length(LPS10_Cilio0_Extent));
x4=10*ones(1,length(LPS10_Cilio20_Extent));
x6=16*ones(1,length(LPS10_Cilio50_Extent));

plot(x1,LPS0_Cilio0_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_Extent),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_Extent),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_Extent,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_Extent),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_Extent,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_Extent),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_Extent,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_Extent),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_Extent,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_Extent),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Cell Extent','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

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
            x1=5;
        elseif sig_extent(ii,1)==3
            x1=8;
        elseif sig_extent(ii,1)==4
            x1=11;
        elseif sig_extent(ii,1)==5
            x1=14;
        else
            x1=17;
        end
        if sig_extent(ii,2)==2
            x2=5;
        elseif sig_extent(ii,2)==3
            x2=8;
        elseif sig_extent(ii,2)==4
            x2=11;
        elseif sig_extent(ii,2)==5
            x2=14;
        else
            x2=17;
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

x1=2*ones(1,length(LPS0_Cilio0_Perimeter));
x3=8*ones(1,length(LPS0_Cilio20_Perimeter));
x5=14*ones(1,length(LPS0_Cilio50_Perimeter));
x2=4*ones(1,length(LPS10_Cilio0_Perimeter));
x4=10*ones(1,length(LPS10_Cilio20_Perimeter));
x6=16*ones(1,length(LPS10_Cilio50_Perimeter));

plot(x1,LPS0_Cilio0_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_Perimeter),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_Perimeter),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_Perimeter,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_Perimeter),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_Perimeter,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_Perimeter),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_Perimeter,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_Perimeter),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_Perimeter,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_Perimeter),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Cell Perimeter (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

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
            x1=5;
        elseif sig_perimeter(ii,1)==3
            x1=8;
        elseif sig_perimeter(ii,1)==4
            x1=11;
        elseif sig_perimeter(ii,1)==5
            x1=14;
        else
            x1=17;
        end
        if sig_perimeter(ii,2)==2
            x2=5;
        elseif sig_perimeter(ii,2)==3
            x2=8;
        elseif sig_perimeter(ii,2)==4
            x2=11;
        elseif sig_perimeter(ii,2)==5
            x2=14;
        else
            x2=17;
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

x1=2*ones(1,length(LPS0_Cilio0_Elongation));
x3=8*ones(1,length(LPS0_Cilio20_Elongation));
x5=14*ones(1,length(LPS0_Cilio50_Elongation));
x2=4*ones(1,length(LPS10_Cilio0_Elongation));
x4=10*ones(1,length(LPS10_Cilio20_Elongation));
x6=16*ones(1,length(LPS10_Cilio50_Elongation));

plot(x1,LPS0_Cilio0_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_Elongation),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_Elongation),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_Elongation,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_Elongation),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_Elongation,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_Elongation),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_Elongation,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_Elongation),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_Elongation,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_Elongation),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Cell Elongation (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

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
            x1=5;
        elseif sig_elongation(ii,1)==3
            x1=8;
        elseif sig_elongation(ii,1)==4
            x1=11;
        elseif sig_elongation(ii,1)==5
            x1=14;
        else
            x1=17;
        end
        if sig_elongation(ii,2)==2
            x2=5;
        elseif sig_elongation(ii,2)==3
            x2=8;
        elseif sig_elongation(ii,2)==4
            x2=11;
        elseif sig_elongation(ii,2)==5
            x2=14;
        else
            x2=17;
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