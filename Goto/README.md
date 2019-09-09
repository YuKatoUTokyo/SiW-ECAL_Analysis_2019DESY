# Personal Directory: Goto

root -b -l -q build_event_goto.C+'(".../run_XXXXX_merge.root")'  
root -b -l -q build_event_all.C+'(".../run_XXXXX_merge.root")'  
root -b -l -q event_display_goto.C+'(".../run_XXXXX_merge_BUILD_EVENT.root")'  
root -b -l -q event_display_all.C+'(".../run_XXXXX_merge_BUILD_EVENT.root")'  
root -b -l -q merge_root_fev.C+'(".../run_XXXXX_")'  
root -b -l -q merge_root_all.C+'(".../run_XXXXX_")'  


### Analysis_Angle/
run_20044 : Change angle sample data   
root -b -l -q cross_talk_analysis.C+'(".../run_XXXXX_merge_BUILD_EVENT.root")'  
#### Analysis_Angle/Result/Cross_Talk/
analysis result for cross talk (ongoing)   

### Analysis_MIP/
root -b -l -q mip_chip.C+'(".../run_XXXXX_merge.root")'  
root -b -l -q mip_map.C+  
root -b -l -q mip_merge.C+'(".../run_XXXXX_merge.root")'  
root -b -l -q mip_merge_range.C+'(".../run_XXXXX_merge.root")'  
root -b -l -q mip_sn.C+  
#### Analysis_MIP/MIP_Map/ (by mip_map.C)
mip_mean_chip.txt : MIP mean information (slab chip channel mip_mean mip_chi mip_ndf)   
MIP_Map.root : MIP information in Tree (mip_mean, mip_chi2ndf)
#### Analysis_MIP/Result/
each analysis result for MIP mean  
##### Analysis_MIP/Result/Slab_Result/ 
##### Analysis_MIP/Result/Chip_Result/
##### Analysis_MIP/Result/Channel_Result/ 
#### Analysis_MIP/S_over_N_ratio/ (by mip_sn.C)
S/N ratio with fllowing cutting conditions,   
mip chi2/ndf < 3   
mip mean > 100 (650 µm)   
mip mean > 50 (320 µm)   

### Analysis_Pedestal/
root -b -l -q ped_memory.C+'(".../run_XXXXX_merge.root")'  
root -b -l -q ped_time.C+  
#### Analysis_Pedestal/Memory_Result/
#### Analysis_Pedestal/Pedestal_Shift/
pedestal difference between ADC mode and TDC mode (fast check)   
the difference depend on memory cells   
#### Analysis_Pedestal/Time_Dependances/
time variation of pedestal and temperature  
##### Analysis_Pedestal/Time_Dependances/Time_Dependance_goto
using the data at 24th June   
##### Analysis_Pedestal/Time_Dependances/Time_Dependance_kato
using the data at 28th June   
##### Analysis_Pedestal/Time_Dependances/Time_Dependance_lowgain
using pedestal calculated by charge_lowGain of ADC mode  
##### Analysis_Pedestal/Time_Dependances/Time_Dependance_rela and slide
pedestal and temperature   

### Analysis_Shower/
run_42003 : Shower sample data with LAL  
run_XXXXX_merge_BUILD_EVENT.root : Shower event (by build_event_all.C)  
run_XXXXX_merge_BUILD_EVENT_EvDisplay.root : TH3F shower event display (by event_display.C)  

### pedetal/
pedestal information in Tree 

### map_chip/
map chip information

### run_data/
run information   
