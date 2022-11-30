function prepareDataDyna
%% PREPARE DATA
load('dataAnalysis.mat','data');

%EXCLUDE CELLS FROM ANALYSIS
data=[data(1:63,:);data(65,:);data(67,:);data(69:70,:);data(72:80,:);data(82,:);data(84:100,:);data(102:111,:);data(155:end,:)];

%SEPARATE DATA ACCORDING TO CONDITION
%Initialise variables.
LPS0_Dyn0=[];
LPS0_Dyn25=[];
LPS10_Dyn0=[];
LPS10_Dyn25=[];

for ii=1:length(data)
    if ~isempty(strfind(data{ii,1},'_10ugLPS+DMSO'))
        LPS10_Dyn0=[LPS10_Dyn0;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'_10ugLPS+Dynarrestin'))
        LPS10_Dyn25=[LPS10_Dyn25;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'_0ugLPS+DMSO'))
        LPS0_Dyn0=[LPS0_Dyn0;data{ii,3}];
    elseif ~isempty(strfind(data{ii,1},'_0ugLPS+Dyanarrestin'))
        LPS0_Dyn25=[LPS0_Dyn25;data{ii,3}];
    end
end

sheet1=struct2table(LPS0_Dyn0);
sheet2=struct2table(LPS0_Dyn25);
sheet3=struct2table(LPS10_Dyn0);
sheet4=struct2table(LPS10_Dyn25);

%Multiply by pixel size to convert pixels to microns.
sheet1.Area=sheet1.Area*0.18;
sheet1.MajorAxisLength=sheet1.MajorAxisLength*0.18;
sheet1.MinorAxisLength=sheet1.MinorAxisLength*0.18;
sheet1.Perimeter=sheet1.Perimeter*0.18;
sheet1.CHA=sheet1.CHA*0.18;
sheet1.CHP=sheet1.CHP*0.18;
sheet1.Density=sheet1.Density*0.18;
sheet1.Roughness=sheet1.Roughness*0.18;
sheet1.Elongation=sheet1.MajorAxisLength./sheet1.MinorAxisLength; %add additional parameter.

sheet2.Area=sheet2.Area*0.18;
sheet2.MajorAxisLength=sheet2.MajorAxisLength*0.18;
sheet2.MinorAxisLength=sheet2.MinorAxisLength*0.18;
sheet2.Perimeter=sheet2.Perimeter*0.18;
sheet2.CHA=sheet2.CHA*0.18;
sheet2.CHP=sheet2.CHP*0.18;
sheet2.Density=sheet2.Density*0.18;
sheet2.Roughness=sheet2.Roughness*0.18;
sheet2.Elongation=sheet2.MajorAxisLength./sheet2.MinorAxisLength; %add additional parameter.

sheet3.Area=sheet3.Area*0.18;
sheet3.MajorAxisLength=sheet3.MajorAxisLength*0.18;
sheet3.MinorAxisLength=sheet3.MinorAxisLength*0.18;
sheet3.Perimeter=sheet3.Perimeter*0.18;
sheet3.CHA=sheet3.CHA*0.18;
sheet3.CHP=sheet3.CHP*0.18;
sheet3.Density=sheet3.Density*0.18;
sheet3.Roughness=sheet3.Roughness*0.18;
sheet3.Elongation=sheet3.MajorAxisLength./sheet3.MinorAxisLength; %add additional parameter.

sheet4.Area=sheet4.Area*0.18;
sheet4.MajorAxisLength=sheet4.MajorAxisLength*0.18;
sheet4.MinorAxisLength=sheet4.MinorAxisLength*0.18;
sheet4.Perimeter=sheet4.Perimeter*0.18;
sheet4.CHA=sheet4.CHA*0.18;
sheet4.CHP=sheet4.CHP*0.18;
sheet4.Density=sheet4.Density*0.18;
sheet4.Roughness=sheet4.Roughness*0.18;
sheet4.Elongation=sheet4.MajorAxisLength./sheet4.MinorAxisLength; %add additional parameter.

writetable(sheet1,'morphologyDataDyna.xlsx','WriteMode','overwritesheet','sheet','LPS0_Dyna0');
writetable(sheet2,'morphologyDataDyna.xlsx','WriteMode','overwritesheet','sheet','LPS0_Dyna25');
writetable(sheet3,'morphologyDataDyna.xlsx','WriteMode','overwritesheet','sheet','LPS10_Dyna0');
writetable(sheet4,'morphologyDataDyna.xlsx','WriteMode','overwritesheet','sheet','LPS10_Dyna25');
end