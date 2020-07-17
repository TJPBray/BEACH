%Joint analysis
%sz vol gives image dimensions
%Slice numbers in 'slice'
%Masks in BW
%Slice names are labstrc

%Get ROI data for analysis (assume image vol already loaded)

%Create masked images
for j=1:numel(slice)

sl=slice(1,j)

Mask=BW{1,j};  %{} used for cell array

Im=volp(:,:,sl);
Im=double(Im); %Converts image to double to allow for multiplication

Maskedim=Im.*Mask;

Maskedim(~Mask) = NaN;

All(:,:,j)=Maskedim;
end

%Calculate threshold parameters for each ROI
for k=1:numel(slice)

   pixelvals=All(:,:,k);
   mean(k,1)=nanmean(pixelvals(:));
   [abovethresh,fract]=inflamquant(pixelvals(:))
   Fraction(k,1)=fract
   Abovethreshold(k,1)=abovethresh
end

% Bar chart
%Get labels
for l=1:numel(slice)
lab=labstrc{1,l} %use 'curly brackets' to extract text from cell array
labels(1,l)=lab
end

   bar(Fraction);
set(gca,'XTickLabel',labels)

