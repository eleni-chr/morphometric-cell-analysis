function visInspect
%% Function written by Eleni Christoforidou in MATLAB R2020b.

%This function is used to visually inspect the image processing performed
%by the custom function "morphology". The function "morphology" must be run
%before running this function.

%Run this function from inside the folder containing the subfolders with
%the FIG files produced by the custom function "morphology".

%% FIND ALL IMAGES TO ANALYSE
list=dir('*/cell*.fig');
wd=cd;

for ii=1:length(list)
    cd(list(ii).folder);
    fn=list(ii).name;
    openfig(fn)
    fprintf('Image number %d\n',ii);
    pause; %press any key to move to the next figure in the loop.
    close
end
cd(wd);
end