# BEACH
MATLAB tool for analysis of qMRI data to quantify oedema and fat

This code implements the tool for ADC measurement described in:
Histographic analysis of oedema and fat in inflamed bone marrow based on qMRI. Bray et al. European Radiology April 2020

The code enables 
(1) semi-automated ROI placement once the observer has defined the joint
(2) histographic analysis of voxel values in the ROI
(3) BEACH_FF_R2star analyses FF and R2* maps simultaneously - the joint is defined on PDFF maps and the ROIs are propagated onto both PDFF maps (visible) and R2* maps (not visible) to enable 2D histogram / density analysis.

To analyse a stack of ADC values, use BEACH_ADC(ADC)

To analyse two stacks of FF and R2* values acquired together (and generated 2D-histograms / density maps), use BEACH_FF_R2(FF,R2star)

BEACH_FF_R2 generates 2D histograms / heatmaps and could be used for cluster analysis or similar although this is not fully implemented at present. 

Examples of all three data types (ADC, FF, R2star) are provided in the testdata folder.

When generating the ROIs, the outer width can be set to 7 and the inner width should be set to 3 (this could also be modified depending on image dimensions if needed). Labels such as R1, L1, R2, L2... can be used to specify the left and right joints on the relevant slices. 

Histographic analysis is implemented by the 'Histogram' pushbutton and can be modified in the corresponding functions in the BEACH_ADC.m and BEACH_FF_R2.m files. 


