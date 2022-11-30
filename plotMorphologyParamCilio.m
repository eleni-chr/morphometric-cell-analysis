function plotMorphologyParamCilio
%% PREPARE DATA
load('dataAnalysis.mat','data');

%EXCLUDE CELLS FROM ANALYSIS
data=[data(1:70,:);data(72:91,:);data(93:106,:);data(108:159,:);data(161:170,:);data(172:end,:)]; %from repeat 1.

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

LPS0_Cilio0_CHA=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_CHA=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_CHA=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_CHA=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_CHA=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_CHA=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_CHP=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_CHP=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_CHP=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_CHP=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_CHP=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_CHP=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_Density=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_Density=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_Density=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_Density=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_Density=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_Density=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_Roughness=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_Roughness=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_Roughness=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_Roughness=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_Roughness=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_Roughness=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_CHspan=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_CHspan=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_CHspan=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_CHspan=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_CHspan=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_CHspan=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_CHC=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_CHC=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_CHC=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_CHC=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_CHC=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_CHC=zeros(length(LPS10_Cilio50),1);

LPS0_Cilio0_FD=zeros(length(LPS0_Cilio0),1);
LPS0_Cilio20_FD=zeros(length(LPS0_Cilio20),1);
LPS0_Cilio50_FD=zeros(length(LPS0_Cilio50),1);
LPS10_Cilio0_FD=zeros(length(LPS10_Cilio0),1);
LPS10_Cilio20_FD=zeros(length(LPS10_Cilio20),1);
LPS10_Cilio50_FD=zeros(length(LPS10_Cilio50),1);

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
    LPS0_Cilio0_CHA(ii,1)=LPS0_Cilio0(ii).CHA*0.18;
    LPS0_Cilio0_CHP(ii,1)=LPS0_Cilio0(ii).CHP*0.18;
    LPS0_Cilio0_Density(ii,1)=LPS0_Cilio0(ii).Density*0.18;
    LPS0_Cilio0_Roughness(ii,1)=LPS0_Cilio0(ii).Roughness*0.18;
    LPS0_Cilio0_CHspan(ii,1)=LPS0_Cilio0(ii).CHspan;
    LPS0_Cilio0_CHC(ii,1)=LPS0_Cilio0(ii).CHC;
    LPS0_Cilio0_FD(ii,1)=LPS0_Cilio0(ii).FD;
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
    LPS0_Cilio20_CHA(ii,1)=LPS0_Cilio20(ii).CHA*0.18;
    LPS0_Cilio20_CHP(ii,1)=LPS0_Cilio20(ii).CHP*0.18;
    LPS0_Cilio20_Density(ii,1)=LPS0_Cilio20(ii).Density*0.18;
    LPS0_Cilio20_Roughness(ii,1)=LPS0_Cilio20(ii).Roughness*0.18;
    LPS0_Cilio20_CHspan(ii,1)=LPS0_Cilio20(ii).CHspan;
    LPS0_Cilio20_CHC(ii,1)=LPS0_Cilio20(ii).CHC;
    LPS0_Cilio20_FD(ii,1)=LPS0_Cilio20(ii).FD;
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
    LPS0_Cilio50_CHA(ii,1)=LPS0_Cilio50(ii).CHA*0.18;
    LPS0_Cilio50_CHP(ii,1)=LPS0_Cilio50(ii).CHP*0.18;
    LPS0_Cilio50_Density(ii,1)=LPS0_Cilio50(ii).Density*0.18;
    LPS0_Cilio50_Roughness(ii,1)=LPS0_Cilio50(ii).Roughness*0.18;
    LPS0_Cilio50_CHspan(ii,1)=LPS0_Cilio50(ii).CHspan;
    LPS0_Cilio50_CHC(ii,1)=LPS0_Cilio50(ii).CHC;
    LPS0_Cilio50_FD(ii,1)=LPS0_Cilio50(ii).FD;
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
    LPS10_Cilio0_CHA(ii,1)=LPS10_Cilio0(ii).CHA*0.18;
    LPS10_Cilio0_CHP(ii,1)=LPS10_Cilio0(ii).CHP*0.18;
    LPS10_Cilio0_Density(ii,1)=LPS10_Cilio0(ii).Density*0.18;
    LPS10_Cilio0_Roughness(ii,1)=LPS10_Cilio0(ii).Roughness*0.18;
    LPS10_Cilio0_CHspan(ii,1)=LPS10_Cilio0(ii).CHspan;
    LPS10_Cilio0_CHC(ii,1)=LPS10_Cilio0(ii).CHC;
    LPS10_Cilio0_FD(ii,1)=LPS10_Cilio0(ii).FD;
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
    LPS10_Cilio20_CHA(ii,1)=LPS10_Cilio20(ii).CHA*0.18;
    LPS10_Cilio20_CHP(ii,1)=LPS10_Cilio20(ii).CHP*0.18;
    LPS10_Cilio20_Density(ii,1)=LPS10_Cilio20(ii).Density*0.18;
    LPS10_Cilio20_Roughness(ii,1)=LPS10_Cilio20(ii).Roughness*0.18;
    LPS10_Cilio20_CHspan(ii,1)=LPS10_Cilio20(ii).CHspan;
    LPS10_Cilio20_CHC(ii,1)=LPS10_Cilio20(ii).CHC;
    LPS10_Cilio20_FD(ii,1)=LPS10_Cilio20(ii).FD;
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
    LPS10_Cilio50_CHA(ii,1)=LPS10_Cilio50(ii).CHA*0.18;
    LPS10_Cilio50_CHP(ii,1)=LPS10_Cilio50(ii).CHP*0.18;
    LPS10_Cilio50_Density(ii,1)=LPS10_Cilio50(ii).Density*0.18;
    LPS10_Cilio50_Roughness(ii,1)=LPS10_Cilio50(ii).Roughness*0.18;
    LPS10_Cilio50_CHspan(ii,1)=LPS10_Cilio50(ii).CHspan;
    LPS10_Cilio50_CHC(ii,1)=LPS10_Cilio50(ii).CHC;
    LPS10_Cilio50_FD(ii,1)=LPS10_Cilio50(ii).FD;
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
    LPS0_Cilio0_CHA=rmoutliers(LPS0_Cilio0_CHA);
    LPS0_Cilio0_CHP=rmoutliers(LPS0_Cilio0_CHP);
    LPS0_Cilio0_Density=rmoutliers(LPS0_Cilio0_Density);
    LPS0_Cilio0_Roughness=rmoutliers(LPS0_Cilio0_Roughness);
    LPS0_Cilio0_CHspan=rmoutliers(LPS0_Cilio0_CHspan);
    LPS0_Cilio0_CHC=rmoutliers(LPS0_Cilio0_CHC);
    LPS0_Cilio0_FD=rmoutliers(LPS0_Cilio0_FD);
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
    LPS0_Cilio20_CHA=rmoutliers(LPS0_Cilio20_CHA);
    LPS0_Cilio20_CHP=rmoutliers(LPS0_Cilio20_CHP);
    LPS0_Cilio20_Density=rmoutliers(LPS0_Cilio20_Density);
    LPS0_Cilio20_Roughness=rmoutliers(LPS0_Cilio20_Roughness);
    LPS0_Cilio20_CHspan=rmoutliers(LPS0_Cilio20_CHspan);
    LPS0_Cilio20_CHC=rmoutliers(LPS0_Cilio20_CHC);
    LPS0_Cilio20_FD=rmoutliers(LPS0_Cilio20_FD);
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
    LPS0_Cilio50_CHA=rmoutliers(LPS0_Cilio50_CHA);
    LPS0_Cilio50_CHP=rmoutliers(LPS0_Cilio50_CHP);
    LPS0_Cilio50_Density=rmoutliers(LPS0_Cilio50_Density);
    LPS0_Cilio50_Roughness=rmoutliers(LPS0_Cilio50_Roughness);
    LPS0_Cilio50_CHspan=rmoutliers(LPS0_Cilio50_CHspan);
    LPS0_Cilio50_CHC=rmoutliers(LPS0_Cilio50_CHC);
    LPS0_Cilio50_FD=rmoutliers(LPS0_Cilio50_FD);
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
    LPS10_Cilio0_CHA=rmoutliers(LPS10_Cilio0_CHA);
    LPS10_Cilio0_CHP=rmoutliers(LPS10_Cilio0_CHP);
    LPS10_Cilio0_Density=rmoutliers(LPS10_Cilio0_Density);
    LPS10_Cilio0_Roughness=rmoutliers(LPS10_Cilio0_Roughness);
    LPS10_Cilio0_CHspan=rmoutliers(LPS10_Cilio0_CHspan);
    LPS10_Cilio0_CHC=rmoutliers(LPS10_Cilio0_CHC);
    LPS10_Cilio0_FD=rmoutliers(LPS10_Cilio0_FD);
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
    LPS10_Cilio20_CHA=rmoutliers(LPS10_Cilio20_CHA);
    LPS10_Cilio20_CHP=rmoutliers(LPS10_Cilio20_CHP);
    LPS10_Cilio20_Density=rmoutliers(LPS10_Cilio20_Density);
    LPS10_Cilio20_Roughness=rmoutliers(LPS10_Cilio20_Roughness);
    LPS10_Cilio20_CHspan=rmoutliers(LPS10_Cilio20_CHspan);
    LPS10_Cilio20_CHC=rmoutliers(LPS10_Cilio20_CHC);
    LPS10_Cilio20_FD=rmoutliers(LPS10_Cilio20_FD);
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
    LPS10_Cilio50_CHA=rmoutliers(LPS10_Cilio50_CHA);
    LPS10_Cilio50_CHP=rmoutliers(LPS10_Cilio50_CHP);
    LPS10_Cilio50_Density=rmoutliers(LPS10_Cilio50_Density);
    LPS10_Cilio50_Roughness=rmoutliers(LPS10_Cilio50_Roughness);
    LPS10_Cilio50_CHspan=rmoutliers(LPS10_Cilio50_CHspan);
    LPS10_Cilio50_CHC=rmoutliers(LPS10_Cilio50_CHC);
    LPS10_Cilio50_FD=rmoutliers(LPS10_Cilio50_FD);
end

%% PERFORM STATISTICS FOR EACH PARAMETER

%Prepare data for Kruskall Wallis tests.
g1=ones(length(LPS0_Cilio0_Area),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_Area),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_Area),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_Area),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_Area),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_Area),1); %condition 6.
group_Area=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_MajorAxisLength),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_MajorAxisLength),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_MajorAxisLength),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_MajorAxisLength),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_MajorAxisLength),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_MajorAxisLength),1); %condition 6.
group_MajorAxisLength=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_MinorAxisLength),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_MinorAxisLength),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_MinorAxisLength),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_MinorAxisLength),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_MinorAxisLength),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_MinorAxisLength),1); %condition 6.
group_MinorAxisLength=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_Eccentricity),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_Eccentricity),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_Eccentricity),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_Eccentricity),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_Eccentricity),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_Eccentricity),1); %condition 6.
group_Eccentricity=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_Circularity),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_Circularity),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_Circularity),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_Circularity),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_Circularity),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_Circularity),1); %condition 6.
group_Circularity=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_Extent),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_Extent),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_Extent),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_Extent),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_Extent),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_Extent),1); %condition 6.
group_Extent=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_Perimeter),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_Perimeter),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_Perimeter),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_Perimeter),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_Perimeter),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_Perimeter),1); %condition 6.
group_Perimeter=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_Elongation),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_Elongation),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_Elongation),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_Elongation),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_Elongation),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_Elongation),1); %condition 6.
group_Elongation=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_CHA),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_CHA),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_CHA),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_CHA),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_CHA),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_CHA),1); %condition 6.
group_CHA=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_CHP),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_CHP),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_CHP),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_CHP),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_CHP),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_CHP),1); %condition 6.
group_CHP=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_Density),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_Density),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_Density),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_Density),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_Density),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_Density),1); %condition 6.
group_Density=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_Roughness),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_Roughness),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_Roughness),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_Roughness),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_Roughness),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_Roughness),1); %condition 6.
group_Roughness=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_CHspan),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_CHspan),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_CHspan),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_CHspan),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_CHspan),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_CHspan),1); %condition 6.
group_CHspan=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_CHC),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_CHC),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_CHC),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_CHC),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_CHC),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_CHC),1); %condition 6.
group_CHC=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

g1=ones(length(LPS0_Cilio0_FD),1); %condition 1.
g2=2*ones(length(LPS10_Cilio0_FD),1); %condition 2.
g3=3*ones(length(LPS0_Cilio20_FD),1); %condition 3.
g4=4*ones(length(LPS10_Cilio20_FD),1); %condition 4.
g5=5*ones(length(LPS0_Cilio50_FD),1); %condition 5.
g6=6*ones(length(LPS10_Cilio50_FD),1); %condition 6.
group_FD=[g1;g2;g3;g4;g5;g6]; %merge conditions in order.

%Merge data in order.
Area=[LPS0_Cilio0_Area;LPS10_Cilio0_Area;LPS0_Cilio20_Area;LPS10_Cilio20_Area;LPS0_Cilio50_Area;LPS10_Cilio50_Area];
MajorAxisLength=[LPS0_Cilio0_MajorAxisLength;LPS10_Cilio0_MajorAxisLength;LPS0_Cilio20_MajorAxisLength;LPS10_Cilio20_MajorAxisLength;LPS0_Cilio50_MajorAxisLength;LPS10_Cilio50_MajorAxisLength];
MinorAxisLength=[LPS0_Cilio0_MinorAxisLength;LPS10_Cilio0_MinorAxisLength;LPS0_Cilio20_MinorAxisLength;LPS10_Cilio20_MinorAxisLength;LPS0_Cilio50_MinorAxisLength;LPS10_Cilio50_MinorAxisLength];
Eccentricity=[LPS0_Cilio0_Eccentricity;LPS10_Cilio0_Eccentricity;LPS0_Cilio20_Eccentricity;LPS10_Cilio20_Eccentricity;LPS0_Cilio50_Eccentricity;LPS10_Cilio50_Eccentricity];
Circularity=[LPS0_Cilio0_Circularity;LPS10_Cilio0_Circularity;LPS0_Cilio20_Circularity;LPS10_Cilio20_Circularity;LPS0_Cilio50_Circularity;LPS10_Cilio50_Circularity];
Extent=[LPS0_Cilio0_Extent;LPS10_Cilio0_Extent;LPS0_Cilio20_Extent;LPS10_Cilio20_Extent;LPS0_Cilio50_Extent;LPS10_Cilio50_Extent];
Perimeter=[LPS0_Cilio0_Perimeter;LPS10_Cilio0_Perimeter;LPS0_Cilio20_Perimeter;LPS10_Cilio20_Perimeter;LPS0_Cilio50_Perimeter;LPS10_Cilio50_Perimeter];
Elongation=[LPS0_Cilio0_Elongation;LPS10_Cilio0_Elongation;LPS0_Cilio20_Elongation;LPS10_Cilio20_Elongation;LPS0_Cilio50_Elongation;LPS10_Cilio50_Elongation];
CHA=[LPS0_Cilio0_CHA;LPS10_Cilio0_CHA;LPS0_Cilio20_CHA;LPS10_Cilio20_CHA;LPS0_Cilio50_CHA;LPS10_Cilio50_CHA];
CHP=[LPS0_Cilio0_CHP;LPS10_Cilio0_CHP;LPS0_Cilio20_CHP;LPS10_Cilio20_CHP;LPS0_Cilio50_CHP;LPS10_Cilio50_CHP];
Density=[LPS0_Cilio0_Density;LPS10_Cilio0_Density;LPS0_Cilio20_Density;LPS10_Cilio20_Density;LPS0_Cilio50_Density;LPS10_Cilio50_Density];
Roughness=[LPS0_Cilio0_Roughness;LPS10_Cilio0_Roughness;LPS0_Cilio20_Roughness;LPS10_Cilio20_Roughness;LPS0_Cilio50_Roughness;LPS10_Cilio50_Roughness];
CHspan=[LPS0_Cilio0_CHspan;LPS10_Cilio0_CHspan;LPS0_Cilio20_CHspan;LPS10_Cilio20_CHspan;LPS0_Cilio50_CHspan;LPS10_Cilio50_CHspan];
CHC=[LPS0_Cilio0_CHC;LPS10_Cilio0_CHC;LPS0_Cilio20_CHC;LPS10_Cilio20_CHC;LPS0_Cilio50_CHC;LPS10_Cilio50_CHC];
FD=[LPS0_Cilio0_FD;LPS10_Cilio0_FD;LPS0_Cilio20_FD;LPS10_Cilio20_FD;LPS0_Cilio50_FD;LPS10_Cilio50_FD];

%Perform Kruskal-Wallis test.
[p_Area,~,stats_Area]=kruskalwallis(Area,group_Area);
[p_MajorAxisLength,~,stats_MajorAxisLength]=kruskalwallis(MajorAxisLength,group_MajorAxisLength);
[p_MinorAxisLength,~,stats_MinorAxisLength]=kruskalwallis(MinorAxisLength,group_MinorAxisLength);
[p_Eccentricity,~,stats_Eccentricity]=kruskalwallis(Eccentricity,group_Eccentricity);
[p_Circularity,~,stats_Circularity]=kruskalwallis(Circularity,group_Circularity);
[p_Extent,~,stats_Extent]=kruskalwallis(Extent,group_Extent);
[p_Perimeter,~,stats_Perimeter]=kruskalwallis(Perimeter,group_Perimeter);
[p_Elongation,~,stats_Elongation]=kruskalwallis(Elongation,group_Elongation);
[p_CHA,~,stats_CHA]=kruskalwallis(CHA,group_CHA);
[p_CHP,~,stats_CHP]=kruskalwallis(CHP,group_CHP);
[p_Density,~,stats_Density]=kruskalwallis(Density,group_Density);
[p_Roughness,~,stats_Roughness]=kruskalwallis(Roughness,group_Roughness);
[p_CHspan,~,stats_CHspan]=kruskalwallis(CHspan,group_CHspan);
[p_CHC,~,stats_CHC]=kruskalwallis(CHC,group_CHC);
[p_FD,~,stats_FD]=kruskalwallis(FD,group_FD);
close all

%Perform multiple comparison tests for all parameters.
c_Area=multcompare(stats_Area,'CType','hsd','Display','off'); %first two columns are the groups being compared. Last column is the p-value.
c_MajorAxisLength=multcompare(stats_MajorAxisLength,'CType','hsd','Display','off');
c_MinorAxisLength=multcompare(stats_MinorAxisLength,'CType','hsd','Display','off');
c_Eccentricity=multcompare(stats_Eccentricity,'CType','hsd','Display','off');
c_Circularity=multcompare(stats_Circularity,'CType','hsd','Display','off');
c_Extent=multcompare(stats_Extent,'CType','hsd','Display','off');
c_Perimeter=multcompare(stats_Perimeter,'CType','hsd','Display','off');
c_Elongation=multcompare(stats_Elongation,'CType','hsd','Display','off');
c_CHA=multcompare(stats_CHA,'CType','hsd','Display','off');
c_CHP=multcompare(stats_CHP,'CType','hsd','Display','off');
c_Density=multcompare(stats_Density,'CType','hsd','Display','off');
c_Roughness=multcompare(stats_Roughness,'CType','hsd','Display','off');
c_CHspan=multcompare(stats_CHspan,'CType','hsd','Display','off');
c_CHC=multcompare(stats_CHC,'CType','hsd','Display','off');
c_FD=multcompare(stats_FD,'CType','hsd','Display','off');

save('KruskalWallisResults.mat','p_Area','p_MajorAxisLength','p_MinorAxisLength','p_Eccentricity','p_Circularity','p_Extent','p_Perimeter','p_Elongation','p_CHA','p_CHP','p_Density','p_Roughness','p_CHspan','p_CHC','p_FD');
save('multCompResults.mat','c_Area','c_MajorAxisLength','c_MinorAxisLength','c_Eccentricity','c_Circularity','c_Extent','c_Perimeter','c_Elongation','c_CHA','c_CHP','c_Density','c_Roughness','c_CHspan','c_CHC','c_FD');

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
ylabel('Cell Area (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

%Indicate significance on plot, if any.
sig=c_Area(:,6)<0.05; %extract only significant results.
sig_Area=c_Area(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_Area,1)
    if sig_Area(ii,6)<0.05
        if sig_Area(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_Area(ii,1)==2
            x1=4;
        elseif sig_Area(ii,1)==3
            x1=8;
        elseif sig_Area(ii,1)==4
            x1=10;
        elseif sig_Area(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_Area(ii,2)==2
            x2=4;
        elseif sig_Area(ii,2)==3
            x2=8;
        elseif sig_Area(ii,2)==4
            x2=10;
        elseif sig_Area(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_Area(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Area(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Area(ii,6)<0.01
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
sig=c_MajorAxisLength(:,6)<0.05; %extract only significant results.
sig_MajorAxisLength=c_MajorAxisLength(sig,:);
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
            x1=8;
        elseif sig_MajorAxisLength(ii,1)==4
            x1=10;
        elseif sig_MajorAxisLength(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_MajorAxisLength(ii,2)==2
            x2=4;
        elseif sig_MajorAxisLength(ii,2)==3
            x2=8;
        elseif sig_MajorAxisLength(ii,2)==4
            x2=10;
        elseif sig_MajorAxisLength(ii,2)==5
            x2=14;
        else
            x2=16;
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
sig=c_MinorAxisLength(:,6)<0.05; %extract only significant results.
sig_MinorAxisLength=c_MinorAxisLength(sig,:);
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
            x1=8;
        elseif sig_MinorAxisLength(ii,1)==4
            x1=10;
        elseif sig_MinorAxisLength(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_MinorAxisLength(ii,2)==2
            x2=4;
        elseif sig_MinorAxisLength(ii,2)==3
            x2=8;
        elseif sig_MinorAxisLength(ii,2)==4
            x2=10;
        elseif sig_MinorAxisLength(ii,2)==5
            x2=14;
        else
            x2=16;
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
sig=c_Eccentricity(:,6)<0.05; %extract only significant results.
sig_Eccentricity=c_Eccentricity(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_Eccentricity,1)
    if sig_Eccentricity(ii,6)<0.05
        if sig_Eccentricity(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_Eccentricity(ii,1)==2
            x1=4;
        elseif sig_Eccentricity(ii,1)==3
            x1=8;
        elseif sig_Eccentricity(ii,1)==4
            x1=10;
        elseif sig_Eccentricity(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_Eccentricity(ii,2)==2
            x2=4;
        elseif sig_Eccentricity(ii,2)==3
            x2=8;
        elseif sig_Eccentricity(ii,2)==4
            x2=10;
        elseif sig_Eccentricity(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_Eccentricity(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Eccentricity(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Eccentricity(ii,6)<0.01
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
sig=c_Circularity(:,6)<0.05; %extract only significant results.
sig_Circularity=c_Circularity(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_Circularity,1)
    if sig_Circularity(ii,6)<0.05
        if sig_Circularity(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_Circularity(ii,1)==2
            x1=4;
        elseif sig_Circularity(ii,1)==3
            x1=8;
        elseif sig_Circularity(ii,1)==4
            x1=10;
        elseif sig_Circularity(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_Circularity(ii,2)==2
            x2=4;
        elseif sig_Circularity(ii,2)==3
            x2=8;
        elseif sig_Circularity(ii,2)==4
            x2=10;
        elseif sig_Circularity(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_Circularity(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Circularity(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Circularity(ii,6)<0.01
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
sig=c_Extent(:,6)<0.05; %extract only significant results.
sig_Extent=c_Extent(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_Extent,1)
    if sig_Extent(ii,6)<0.05
        if sig_Extent(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_Extent(ii,1)==2
            x1=4;
        elseif sig_Extent(ii,1)==3
            x1=8;
        elseif sig_Extent(ii,1)==4
            x1=10;
        elseif sig_Extent(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_Extent(ii,2)==2
            x2=4;
        elseif sig_Extent(ii,2)==3
            x2=8;
        elseif sig_Extent(ii,2)==4
            x2=10;
        elseif sig_Extent(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_Extent(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Extent(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Extent(ii,6)<0.01
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
sig=c_Perimeter(:,6)<0.05; %extract only significant results.
sig_Perimeter=c_Perimeter(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_Perimeter,1)
    if sig_Perimeter(ii,6)<0.05
        if sig_Perimeter(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_Perimeter(ii,1)==2
            x1=4;
        elseif sig_Perimeter(ii,1)==3
            x1=8;
        elseif sig_Perimeter(ii,1)==4
            x1=10;
        elseif sig_Perimeter(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_Perimeter(ii,2)==2
            x2=4;
        elseif sig_Perimeter(ii,2)==3
            x2=8;
        elseif sig_Perimeter(ii,2)==4
            x2=10;
        elseif sig_Perimeter(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_Perimeter(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Perimeter(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Perimeter(ii,6)<0.01
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
sig=c_Elongation(:,6)<0.05; %extract only significant results.
sig_Elongation=c_Elongation(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_Elongation,1)
    if sig_Elongation(ii,6)<0.05
        if sig_Elongation(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_Elongation(ii,1)==2
            x1=4;
        elseif sig_Elongation(ii,1)==3
            x1=8;
        elseif sig_Elongation(ii,1)==4
            x1=10;
        elseif sig_Elongation(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_Elongation(ii,2)==2
            x2=4;
        elseif sig_Elongation(ii,2)==3
            x2=8;
        elseif sig_Elongation(ii,2)==4
            x2=10;
        elseif sig_Elongation(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_Elongation(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Elongation(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Elongation(ii,6)<0.01
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

%% PLOT CHA FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Cilio0_CHA));
x3=8*ones(1,length(LPS0_Cilio20_CHA));
x5=14*ones(1,length(LPS0_Cilio50_CHA));
x2=4*ones(1,length(LPS10_Cilio0_CHA));
x4=10*ones(1,length(LPS10_Cilio20_CHA));
x6=16*ones(1,length(LPS10_Cilio50_CHA));

plot(x1,LPS0_Cilio0_CHA,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_CHA),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_CHA,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_CHA),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_CHA,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_CHA),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_CHA,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_CHA),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_CHA,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_CHA),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_CHA,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_CHA),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Convex hull area (μm^2)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

%Indicate significance on plot, if any.
sig=c_CHA(:,6)<0.05; %extract only significant results.
sig_CHA=c_CHA(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_CHA,1)
    if sig_CHA(ii,6)<0.05
        if sig_CHA(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_CHA(ii,1)==2
            x1=4;
        elseif sig_CHA(ii,1)==3
            x1=8;
        elseif sig_CHA(ii,1)==4
            x1=10;
        elseif sig_CHA(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_CHA(ii,2)==2
            x2=4;
        elseif sig_CHA(ii,2)==3
            x2=8;
        elseif sig_CHA(ii,2)==4
            x2=10;
        elseif sig_CHA(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_CHA(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_CHA(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_CHA(ii,6)<0.01
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
savefig(fig,'CellCHA'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'CellCHA.svg'); %save figure as an SVG file.
saveas(fig,'CellCHA.jpg'); %save figure as an JPEG file.
close

%% PLOT CHP FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Cilio0_CHP));
x3=8*ones(1,length(LPS0_Cilio20_CHP));
x5=14*ones(1,length(LPS0_Cilio50_CHP));
x2=4*ones(1,length(LPS10_Cilio0_CHP));
x4=10*ones(1,length(LPS10_Cilio20_CHP));
x6=16*ones(1,length(LPS10_Cilio50_CHP));

plot(x1,LPS0_Cilio0_CHP,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_CHP),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_CHP,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_CHP),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_CHP,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_CHP),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_CHP,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_CHP),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_CHP,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_CHP),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_CHP,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_CHP),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Convex hull perimeter (μm)','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

%Indicate significance on plot, if any.
sig=c_CHP(:,6)<0.05; %extract only significant results.
sig_CHP=c_CHP(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_CHP,1)
    if sig_CHP(ii,6)<0.05
        if sig_CHP(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_CHP(ii,1)==2
            x1=4;
        elseif sig_CHP(ii,1)==3
            x1=8;
        elseif sig_CHP(ii,1)==4
            x1=10;
        elseif sig_CHP(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_CHP(ii,2)==2
            x2=4;
        elseif sig_CHP(ii,2)==3
            x2=8;
        elseif sig_CHP(ii,2)==4
            x2=10;
        elseif sig_CHP(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_CHP(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_CHP(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_CHP(ii,6)<0.01
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
savefig(fig,'CellCHP'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'CellCHP.svg'); %save figure as an SVG file.
saveas(fig,'CellCHP.jpg'); %save figure as an JPEG file.
close

%% PLOT Density FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Cilio0_Density));
x3=8*ones(1,length(LPS0_Cilio20_Density));
x5=14*ones(1,length(LPS0_Cilio50_Density));
x2=4*ones(1,length(LPS10_Cilio0_Density));
x4=10*ones(1,length(LPS10_Cilio20_Density));
x6=16*ones(1,length(LPS10_Cilio50_Density));

plot(x1,LPS0_Cilio0_Density,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_Density),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_Density,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_Density),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_Density,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_Density),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_Density,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_Density),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_Density,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_Density),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_Density,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_Density),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Cell density','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

%Indicate significance on plot, if any.
sig=c_Density(:,6)<0.05; %extract only significant results.
sig_Density=c_Density(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_Density,1)
    if sig_Density(ii,6)<0.05
        if sig_Density(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_Density(ii,1)==2
            x1=4;
        elseif sig_Density(ii,1)==3
            x1=8;
        elseif sig_Density(ii,1)==4
            x1=10;
        elseif sig_Density(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_Density(ii,2)==2
            x2=4;
        elseif sig_Density(ii,2)==3
            x2=8;
        elseif sig_Density(ii,2)==4
            x2=10;
        elseif sig_Density(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_Density(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Density(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Density(ii,6)<0.01
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
savefig(fig,'CellDensity'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'CellDensity.svg'); %save figure as an SVG file.
saveas(fig,'CellDensity.jpg'); %save figure as an JPEG file.
close

%% PLOT Roughness FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Cilio0_Roughness));
x3=8*ones(1,length(LPS0_Cilio20_Roughness));
x5=14*ones(1,length(LPS0_Cilio50_Roughness));
x2=4*ones(1,length(LPS10_Cilio0_Roughness));
x4=10*ones(1,length(LPS10_Cilio20_Roughness));
x6=16*ones(1,length(LPS10_Cilio50_Roughness));

plot(x1,LPS0_Cilio0_Roughness,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_Roughness),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_Roughness,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_Roughness),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_Roughness,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_Roughness),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_Roughness,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_Roughness),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_Roughness,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_Roughness),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_Roughness,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_Roughness),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Cell roughness','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

%Indicate significance on plot, if any.
sig=c_Roughness(:,6)<0.05; %extract only significant results.
sig_Roughness=c_Roughness(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_Roughness,1)
    if sig_Roughness(ii,6)<0.05
        if sig_Roughness(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_Roughness(ii,1)==2
            x1=4;
        elseif sig_Roughness(ii,1)==3
            x1=8;
        elseif sig_Roughness(ii,1)==4
            x1=10;
        elseif sig_Roughness(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_Roughness(ii,2)==2
            x2=4;
        elseif sig_Roughness(ii,2)==3
            x2=8;
        elseif sig_Roughness(ii,2)==4
            x2=10;
        elseif sig_Roughness(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_Roughness(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Roughness(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_Roughness(ii,6)<0.01
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
savefig(fig,'CellRoughness'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'CellRoughness.svg'); %save figure as an SVG file.
saveas(fig,'CellRoughness.jpg'); %save figure as an JPEG file.
close

%% PLOT CHspan FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Cilio0_CHspan));
x3=8*ones(1,length(LPS0_Cilio20_CHspan));
x5=14*ones(1,length(LPS0_Cilio50_CHspan));
x2=4*ones(1,length(LPS10_Cilio0_CHspan));
x4=10*ones(1,length(LPS10_Cilio20_CHspan));
x6=16*ones(1,length(LPS10_Cilio50_CHspan));

plot(x1,LPS0_Cilio0_CHspan,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_CHspan),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_CHspan,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_CHspan),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_CHspan,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_CHspan),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_CHspan,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_CHspan),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_CHspan,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_CHspan),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_CHspan,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_CHspan),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Convex hull span','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

%Indicate significance on plot, if any.
sig=c_CHspan(:,6)<0.05; %extract only significant results.
sig_CHspan=c_CHspan(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_CHspan,1)
    if sig_CHspan(ii,6)<0.05
        if sig_CHspan(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_CHspan(ii,1)==2
            x1=4;
        elseif sig_CHspan(ii,1)==3
            x1=8;
        elseif sig_CHspan(ii,1)==4
            x1=10;
        elseif sig_CHspan(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_CHspan(ii,2)==2
            x2=4;
        elseif sig_CHspan(ii,2)==3
            x2=8;
        elseif sig_CHspan(ii,2)==4
            x2=10;
        elseif sig_CHspan(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_CHspan(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_CHspan(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_CHspan(ii,6)<0.01
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
savefig(fig,'CellCHspan'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'CellCHspan.svg'); %save figure as an SVG file.
saveas(fig,'CellCHspan.jpg'); %save figure as an JPEG file.
close

%% PLOT CHC FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Cilio0_CHC));
x3=8*ones(1,length(LPS0_Cilio20_CHC));
x5=14*ones(1,length(LPS0_Cilio50_CHC));
x2=4*ones(1,length(LPS10_Cilio0_CHC));
x4=10*ones(1,length(LPS10_Cilio20_CHC));
x6=16*ones(1,length(LPS10_Cilio50_CHC));

plot(x1,LPS0_Cilio0_CHC,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_CHC),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_CHC,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_CHC),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_CHC,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_CHC),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_CHC,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_CHC),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_CHC,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_CHC),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_CHC,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_CHC),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Convex hull circularity','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

%Indicate significance on plot, if any.
sig=c_CHC(:,6)<0.05; %extract only significant results.
sig_CHC=c_CHC(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_CHC,1)
    if sig_CHC(ii,6)<0.05
        if sig_CHC(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_CHC(ii,1)==2
            x1=4;
        elseif sig_CHC(ii,1)==3
            x1=8;
        elseif sig_CHC(ii,1)==4
            x1=10;
        elseif sig_CHC(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_CHC(ii,2)==2
            x2=4;
        elseif sig_CHC(ii,2)==3
            x2=8;
        elseif sig_CHC(ii,2)==4
            x2=10;
        elseif sig_CHC(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_CHC(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_CHC(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_CHC(ii,6)<0.01
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
savefig(fig,'CellCHC'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'CellCHC.svg'); %save figure as an SVG file.
saveas(fig,'CellCHC.jpg'); %save figure as an JPEG file.
close

%% PLOT FD FOR ALL CONDITIONS

fig=figure;

x1=2*ones(1,length(LPS0_Cilio0_FD));
x3=8*ones(1,length(LPS0_Cilio20_FD));
x5=14*ones(1,length(LPS0_Cilio50_FD));
x2=4*ones(1,length(LPS10_Cilio0_FD));
x4=10*ones(1,length(LPS10_Cilio20_FD));
x6=16*ones(1,length(LPS10_Cilio50_FD));

plot(x1,LPS0_Cilio0_FD,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 1.
set(gcf,'WindowStyle','docked'); %dock figure.
hold on;
plot(x1,mean(LPS0_Cilio0_FD),'kx','MarkerSize',18); %plot mean for condition 1.
plot(x3,LPS0_Cilio20_FD,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 2.
plot(x3,mean(LPS0_Cilio20_FD),'kx','MarkerSize',18); %plot mean for condition 2.
plot(x5,LPS0_Cilio50_FD,'ko','MarkerEdgeColor',[.75 0 .75],'MarkerSize',5); %plot condition 3.
plot(x5,mean(LPS0_Cilio50_FD),'kx','MarkerSize',18); %plot mean for condition 3.
plot(x2,LPS10_Cilio0_FD,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 4.
plot(x2,mean(LPS10_Cilio0_FD),'kx','MarkerSize',18); %plot mean for condition 4.
plot(x4,LPS10_Cilio20_FD,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 5.
plot(x4,mean(LPS10_Cilio20_FD),'kx','MarkerSize',18); %plot mean for condition 5.
plot(x6,LPS10_Cilio50_FD,'ko','MarkerEdgeColor',[0 .25 .25],'MarkerSize',5); %plot condition 6.
plot(x6,mean(LPS10_Cilio50_FD),'kx','MarkerSize',18); %plot mean for condition 6.

%set(gca,'box','off'); %remove top x-axis and right y-axis.
xlim([0 17]);
xticks(0:17);
ylabel('Fractal dimension','FontSize',12); %add y-axis label.
gca.FontSize=12; %set axes fontsize.
xticklabels({'','','','0 μg/ml cilioD','','','','','','20 μg/ml cilioD','','','','','','50 μg/ml cilioD',''});

%Indicate significance on plot, if any.
sig=c_FD(:,6)<0.05; %extract only significant results.
sig_FD=c_FD(sig,:);
yt=yticks; %extract values of yticks.
ydiff=yt(end)-yt(end-1); %size of yticks.
a=ydiff;
for ii=1:size(sig_FD,1)
    if sig_FD(ii,6)<0.05
        if sig_FD(ii,1)==1 %find x-coordinates of significance line.
            x1=2;
        elseif sig_FD(ii,1)==2
            x1=4;
        elseif sig_FD(ii,1)==3
            x1=8;
        elseif sig_FD(ii,1)==4
            x1=10;
        elseif sig_FD(ii,1)==5
            x1=14;
        else
            x1=16;
        end
        if sig_FD(ii,2)==2
            x2=4;
        elseif sig_FD(ii,2)==3
            x2=8;
        elseif sig_FD(ii,2)==4
            x2=10;
        elseif sig_FD(ii,2)==5
            x2=14;
        else
            x2=16;
        end
        
        %Plot significance line.
        ylim([0,yt(end)+ydiff]); %increase limit of y-axis.
        y=[max(ylim),max(ylim)]; %y-cooridnates of line.
        x=[x1,x2]; %x-cooridnates of line.
        plot(x,y,'k'); %plot significance line.
        ydiff=ydiff+a; %increase y-coordinates for next line.
        
        %Plot significance asterisks.
        if sig_FD(ii,6)<0.0001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25,x(1)+1.75]; %x-coordinates of asterisks.
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_FD(ii,6)<0.001
            x=[x(1)+0.25,x(1)+0.75,x(1)+1.25]; %x-coordinates of asterisks.
            
            
            %Find y-coordinates.
            yrange=a+a;
            y8=yrange/8;
            ypos=max(ylim)+y8;
            y=[ypos,ypos,ypos]; %y-coordinates of asterisks.
            
            plot(x,y,'k*'); %plot asterisks.
        elseif sig_FD(ii,6)<0.01
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
savefig(fig,'CellFD'); %save figure as a FIG file in the working directory.
fig.Renderer='painters'; %force MATLAB to render the image as a vector.
saveas(fig,'CellFD.svg'); %save figure as an SVG file.
saveas(fig,'CellFD.jpg'); %save figure as an JPEG file.
close
end