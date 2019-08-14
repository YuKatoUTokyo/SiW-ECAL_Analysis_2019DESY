# Personal Directory: Goto

root -b -l -q build_event_goto.C+'("Data/run_20044_merge.root")'  
root -b -l -q build_event_all.C+'("Data/run_20044_merge.root")'  
root -b -l -q event_display_goto.C+'("Data/run_20044_merge_BUILD_EVENT.root")'  
root -b -l -q event_display_all.C+'("Data/run_20044_merge_BUILD_EVENT.root")'  
root -b -l -q merge_analysis.C+'("Data/run_20044_merge.root")'  
root -b -l -q merge_root_fev.C+'("Data/run_20044_")'  
root -b -l -q merge_root_all.C+'("Data/run_20044_")'  

