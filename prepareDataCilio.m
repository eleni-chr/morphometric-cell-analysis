function prepareDataCilio
%% PREPARE DATA
load('dataAnalysis.mat','data');

%EXCLUDE CELLS FROM ANALYSIS
data=[data(1:70,:);data(72:91,:);data(93:106,:);data(108:159,:);data(161:170,:);data(172:end,:)];

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

sheet1=struct2table(LPS0_Cilio0);
sheet2=struct2table(LPS0_Cilio20);
sheet3=struct2table(LPS0_Cilio50);
sheet4=struct2table(LPS10_Cilio0);
sheet5=struct2table(LPS10_Cilio20);
sheet6=struct2table(LPS10_Cilio50);

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

sheet5.Area=sheet5.Area*0.18;
sheet5.MajorAxisLength=sheet5.MajorAxisLength*0.18;
sheet5.MinorAxisLength=sheet5.MinorAxisLength*0.18;
sheet5.Perimeter=sheet5.Perimeter*0.18;
sheet5.CHA=sheet5.CHA*0.18;
sheet5.CHP=sheet5.CHP*0.18;
sheet5.Density=sheet5.Density*0.18;
sheet5.Roughness=sheet5.Roughness*0.18;
sheet5.Elongation=sheet5.MajorAxisLength./sheet5.MinorAxisLength; %add additional parameter.

sheet6.Area=sheet6.Area*0.18;
sheet6.MajorAxisLength=sheet6.MajorAxisLength*0.18;
sheet6.MinorAxisLength=sheet6.MinorAxisLength*0.18;
sheet6.Perimeter=sheet6.Perimeter*0.18;
sheet6.CHA=sheet6.CHA*0.18;
sheet6.CHP=sheet6.CHP*0.18;
sheet6.Density=sheet6.Density*0.18;
sheet6.Roughness=sheet6.Roughness*0.18;
sheet6.Elongation=sheet6.MajorAxisLength./sheet6.MinorAxisLength; %add additional parameter.

writetable(sheet1,'morphologyDataCilio.xlsx','WriteMode','overwritesheet','sheet','LPS0_Cilio0');
writetable(sheet2,'morphologyDataCilio.xlsx','WriteMode','overwritesheet','sheet','LPS0_Cilio20');
writetable(sheet3,'morphologyDataCilio.xlsx','WriteMode','overwritesheet','sheet','LPS0_Cilio50');
writetable(sheet4,'morphologyDataCilio.xlsx','WriteMode','overwritesheet','sheet','LPS10_Cilio0');
writetable(sheet5,'morphologyDataCilio.xlsx','WriteMode','overwritesheet','sheet','LPS10_Cilio20');
writetable(sheet6,'morphologyDataCilio.xlsx','WriteMode','overwritesheet','sheet','LPS10_Cilio50');
end