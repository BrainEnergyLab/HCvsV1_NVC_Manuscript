# HCvsV1_NVC_Manuscript

All scripts used to analyse data for the Shaw et al. manuscript: "Hippocampus has lower oxygenation and weaker control of brain blood flow than cortex, due to microvascular differences."

# Figure 2

A-C: extractProbeData.mat, extractRestAndLocoData_probe.mat, avgRestAndLocoData_probe 
Run these code in the order above to convert the raw oxy-CBF probe data into a mat file, and then find the average CMRO2, CBF and SO2 values during rest periods. 

F-H: extractLineScanData.mat, findDiam_lineScan.mat, linescanDiamAnalysis.mat, linescanVelocityAnalysis.mat, findAvg_RBCV_diam.mat, compareGrps_avgRBCV_diam.mat
Run the top function 'extractLineScanData.mat' to call on the subfunctions which extract the diameter and RBCV from linescan 2P data. 
Once RBCV data is extracted, call 'findAvg_RBCV_diam.mat' to average the RBCV during rest periods for individual vessels, and then use 'compareGrps_avgRBCV_diam.mat' to compare all average RBCV readings between the HC and V1. 

# Figure 3


