# Personal Directory: Goto

root -b -l -q build_event_goto.C+'("Data/run_20044_merge.root")'  
root -b -l -q build_event_all.C+'("Data/run_20044_merge.root")'  
root -b -l -q event_display_goto.C+'("Data/run_20044_merge_BUILD_EVENT.root")'  
root -b -l -q event_display_all.C+'("Data/run_20044_merge_BUILD_EVENT.root")'  
root -b -l -q merge_analysis.C+'("Data/run_20044_merge.root")'  
root -b -l -q merge_root_fev.C+'("Data/run_20044_")'  
root -b -l -q merge_root_all.C+'("Data/run_20044_")'  


### ShowerSampleData/
run_42003 : Shower sample data with LAL  
run_XXXXX_merge_BUILD_EVENT.root : Shower event by build_event_all.C  
run_XXXXX_merge_BUILD_EVENT_EvDisplay.root : TH3F shower event display by event_display.C  

### AngleSampleData/
run_20044 : Change angle sample data

### Result/
run_XXXXX_merge_analysis.pdf : Analysis for cross talk by merge_analysis.C  

### Pedestal_Map/
Pedestal tree for Only FEV13 in Kyushu  

### map_chip/
map chip information
