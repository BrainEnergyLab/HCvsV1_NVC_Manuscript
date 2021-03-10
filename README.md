# HCvsV1_NVC_Manuscript

All scripts used to analyse data for the Shaw et al. manuscript: "Hippocampus has lower oxygenation and weaker control of brain blood flow than cortex, due to microvascular differences."
https://www.biorxiv.org/content/10.1101/835728v1 


# Figure 2

A-C: extractProbeData.mat, extractRestAndLocoData_probe.mat, avgRestAndLocoData_probe 
Run these code in the order above to convert the raw oxy-CBF probe data into a mat file, and then find the average CMRO2, CBF and SO2 values during rest periods. 

*NB/ The code to extract baseline haemodynamic data during rest for figure 2, is also used in supplementary figure 10.*

F-J: extractLineScanData.mat, findDiam_lineScan.mat, linescanDiamAnalysis.mat, linescanVelocityAnalysis.mat, findAvg_RBCV_diam.mat, compareGrps_avgRBCV_diam.mat, extractHaematocrit.mat, extractRBCflux.mat
Run the top function 'extractLineScanData.mat' to call on the subfunctions which extract the diameter and RBCV from linescan 2P data. Call 'extractHaematocrit.mat' to extract the haematocrit over time, and 'extractRBCflux.mat' to extract the RBCs per second. 
Once linescan data is extracted, call 'getAvgLinescanAll.mat' to get the average values for each measure during rest periods (i.e. removes locomotion, and finds mean, std, and coefficient of variation). 

*NB/ The code to extract linescan diameter and velocity traces for figure 2, is also used in supplementary figures 1 & 7.*


# Figure 3

xyFWHM.mat, cutByCalcEvents_xy_longer.mat, findSignalSpikesInTrace.mat, avgCalcPeaks_Fig3.mat 

The 2P vessel diameter data can be extracted using 'xyFWHM.mat' - this loads the vessel channel into matlab, and asks the user to specify the vessel skeleton and vessel branches. The code will then run a perpendicular line along the vessel skeleton and extract and average the full width half maximum for points along this skeleton - to get the vessel diameter (converted to um from the pixel size info). 

The calcium data was extracted using the CellSort package in matlab. Mukamel, E. A., Nimmerjahn, A. & Schnitzer, M. J. Automated analysis of cellular signals from large-scale calcium imaging data. Neuron 63, 747–60 (2009). 


Once vessel and calcium time series were extracted, the calcium peaks were detected using the 'findSignalSpikesInTrace.mat' function, which is called from the top script 'cutByCalcEvents_xy'. This script also then cuts around all calcium peaks for neuronal calcium traces and corresponding vessel diameter, locomotion and stimulus traces. 


The 'avgCalcPeaks' function takes all the detected calcium peaks across experimental directories, and sorts them to find responsive vessel dilations, categorise by locomotion/stimulation condition, neuron type or brain layer, etc. 

*NB/ The code to extract xy vessel diameter and neuronal calcium traces for figure 3, is also used in supplementary figures 2-4 & supplementary figures 6-7.*


# Figure 5

entireNetTrace_PPM_szPks.mat, findPeaksInTrace.mat, cutByCMRO2Events_probe.mat, avgCMRO2Peaks.mat 


A-D: Run the top function 'entireNetTrace_PPM_szPks.mat', and this will find the calcium peaks in the net calcium trace for each wide-field recording inside the top directory. It will call the subfunction 'findPeaksInTrace.mat' in order to detect the peaks in the net trace used for the analysis. Once the peaks have been detected, the size of the net trace peaks will be averaged across all detected peaks. And the correlation between individual ROI traces will be correlated. 


E-J: First call the function 'cutByCMRO2Events_probe.mat' to find all the peaks in the CMRO2 traces for each recording within the top directory. This will cut around these CMRO2 events, and also the corresponding haemodynamic parameters, i.e. Hbt. Once you have extracted these CMRO2 peaks for both groups you wish to compare (i.e. HC and V1) - call the function 'avgCMRO2Peaks.mat' to compare the average size of the CMRO2 traces and corresponding Hbt traces between groups. 


# Figure 7 

capSpacing_fromDistMap.mat, runpdesinglecapDwall.mat, pdesinglecapDwall.mat

Run the 'capSpacing_fromDistMap.mat' to find all the distance maps in the experimental directory and find the capillary spacing for multiple percentiles. This will then get average capillary spacings per recording, and then an overall average (from these averages) per region. This cap spacing info can then be fed into the oxygen diffusion model. 

Run the top function 'runpdesinglecapDwall.mat' to run the oxygen diffusion model for both regions, based on the capillary spacings inputted. This will call the subfunction 'pdesinglecapDwall.mat' which calculates the oxygen diffusion across the tissue based on a single capillary with a wall in an arena with specified spacings to the nearest capillary. 


# Supplementary Figure 1

extractLScalcium.mat, cutByLScalc.mat

Run 'extractLScalcium.mat' to search for the calcium.tif file and create a trace of light intensity changes over time (i.e. when calcium activity increases vs baseline). You must also have run the mat file to extract RBCV or diameter over time within the same directory (NB/ see info for figure 2). 

Use the cutByLScalc.mat file to detect peaks in the calcium, and cut around corresponding RBCV or diameter traces (need to edit prefs / matfile search to look for velocity or diameter file). 





