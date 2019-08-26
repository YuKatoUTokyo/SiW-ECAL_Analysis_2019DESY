# Personal Directory: Goto

root -b -l -q build_event_goto.C+'(".../run_XXXXX_merge.root")'  
root -b -l -q build_event_all.C+'(".../run_XXXXX_merge.root")'  
root -b -l -q event_display_goto.C+'(".../run_XXXXX_merge_BUILD_EVENT.root")'  
root -b -l -q event_display_all.C+'(".../run_XXXXX_merge_BUILD_EVENT.root")'  
root -b -l -q merge_analysis.C+'(".../run_XXXXX_merge.root")'  
root -b -l -q merge_root_fev.C+'(".../run_XXXXX_")'  
root -b -l -q merge_root_all.C+'(".../run_XXXXX_")'  


### ShowerSampleData/
run_42003 : Shower sample data with LAL  
run_XXXXX_merge_BUILD_EVENT.root : Shower event (by build_event_all.C)  
run_XXXXX_merge_BUILD_EVENT_EvDisplay.root : TH3F shower event display (by event_display.C)  

### AngleSampleData/
run_20044 : Change angle sample data

### ChangeAngleResult/
run_XXXXX_merge_analysis.pdf : Analysis for cross talk and sensitive thickness (by merge_analysis.C)  

### pedetal/
Pedestal tree 

### map_chip/
map chip information
