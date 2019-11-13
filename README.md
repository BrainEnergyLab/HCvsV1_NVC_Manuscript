# HCvsV1_NVC_Manuscript

All scripts used to analyse data for the Shaw et al. manuscript: "Hippocampus has lower oxygenation and weaker control of brain blood flow than cortex, due to microvascular differences."
https://www.biorxiv.org/content/10.1101/835728v1 


# Figure 2

A-C: extractProbeData.mat, extractRestAndLocoData_probe.mat, avgRestAndLocoData_probe 
Run these code in the order above to convert the raw oxy-CBF probe data into a mat file, and then find the average CMRO2, CBF and SO2 values during rest periods. 

F-H: extractLineScanData.mat, findDiam_lineScan.mat, linescanDiamAnalysis.mat, linescanVelocityAnalysis.mat, findAvg_RBCV_diam.mat, compareGrps_avgRBCV_diam.mat
Run the top function 'extractLineScanData.mat' to call on the subfunctions which extract the diameter and RBCV from linescan 2P data. 
Once RBCV data is extracted, call 'findAvg_RBCV_diam.mat' to average the RBCV during rest periods for individual vessels, and then use 'compareGrps_avgRBCV_diam.mat' to compare all average RBCV readings between the HC and V1. 


# Figure 3

xyFWHM.mat, cutByCalcEvents_xy_longer.mat, findSignalSpikesInTrace.mat, avgCalcPeaks_Fig3.mat 

The 2P vessel diameter data can be extracted using 'xyFWHM.mat' - this loads the vessel channel into matlab, and asks the user to specify the vessel skeleton and vessel branches. The code will then run a perpendicular line along the vessel skeleton and extract and average the full width half maximum for points along this skeleton - to get the vessel diameter (converted to um from the pixel size info). 
The calcium data was extracted using the CellSort package in matlab. Mukamel, E. A., Nimmerjahn, A. & Schnitzer, M. J. Automated analysis of cellular signals from large-scale calcium imaging data. Neuron 63, 747–60 (2009). 


Once vessel and calcium time series were extracted, the calcium peaks were detected using the 'findSignalSpikesInTrace.mat' function, which is called from the top script 'cutByCalcEvents_xy'. This script also then cuts around all calcium peaks for neuronal calcium traces and corresponding vessel diameter, locomotion and stimulus traces. 


The 'avgCalcPeaks' function takes all the detected calcium peaks across experimental directories, and sorts them to find responsive vessel dilations, categorise by locomotion/stimulation condition, neuron type or brain layer, etc. 




