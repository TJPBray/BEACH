%Import mDixon Quant FF map and convert to usable FF and R2* volumes
%Colourmap as normal

clear all

%Select folder
foldername=uigetdir;
folderinfo=dir(foldername);
n=numel(folderinfo)-2;

%Acquire images from DICOM files
for k=1:n
    filename=folderinfo(k+2).name;
    filename=fullfile(foldername,filename);
    image=dicomread(filename);
    stack(:,:,k)=image;
    
    %try
    %info(:,:,k)=dicominfo(filename);
    %catch
    %    ;
    %end
    
end

newanal2(stack)