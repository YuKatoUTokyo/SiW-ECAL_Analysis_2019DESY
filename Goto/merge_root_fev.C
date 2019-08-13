#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <bitset>
#include <math.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCut.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TH1.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"



using namespace::std;

void merge_root_fev(std::string str, bool filterEvents = false){

  TString fullpath = str.data();

  // Set Loop Max
  Int_t MaxSlab = 5;
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxMemory = 15;
  Int_t Entries[MaxSlab];
  Int_t Entry[MaxSlab];
  Int_t MaxEntries = 0;

  TFile *f[MaxSlab];
  TTree *tree[MaxSlab];

  // New TTree Components
  Int_t acqNumber;
  Int_t chipid[MaxSlab][MaxChip];
  Int_t numCol[MaxSlab][MaxChip];
  Int_t nhits[MaxSlab][MaxChip][MaxMemory];
  Int_t bcid[MaxSlab][MaxChip][MaxMemory];
  Int_t corrected_bcid[MaxSlab][MaxChip][MaxMemory];
  Int_t badbcid[MaxSlab][MaxChip][MaxMemory];
  Int_t charge_lowGain[MaxSlab][MaxChip][MaxMemory][MaxChannel];
  Int_t charge_hiGain[MaxSlab][MaxChip][MaxMemory][MaxChannel];
  Int_t gain_hit_low[MaxSlab][MaxChip][MaxMemory][MaxChannel];
  Int_t gain_hit_high[MaxSlab][MaxChip][MaxMemory][MaxChannel];

  Int_t acqNumber_in[MaxSlab];
  Int_t chipid_in[MaxSlab][MaxChip];
  Int_t numCol_in[MaxSlab][MaxChip];
  Int_t nhits_in[MaxSlab][MaxChip][MaxMemory];
  Int_t bcid_in[MaxSlab][MaxChip][MaxMemory];
  Int_t corrected_bcid_in[MaxSlab][MaxChip][MaxMemory];
  Int_t badbcid_in[MaxSlab][MaxChip][MaxMemory];
  Int_t charge_lowGain_in[MaxSlab][MaxChip][MaxMemory][MaxChannel];
  Int_t charge_hiGain_in[MaxSlab][MaxChip][MaxMemory][MaxChannel];
  Int_t gain_hit_low_in[MaxSlab][MaxChip][MaxMemory][MaxChannel];
  Int_t gain_hit_high_in[MaxSlab][MaxChip][MaxMemory][MaxChannel];

  // ================================================================================================================================================ //
  for(Int_t islab=0; islab<MaxSlab; islab++){
    Entries[islab] = 0;
    Entry[islab] = 0;

    TString dif="";
    if(islab == 0) dif = "dif_1_1_1";
    if(islab == 1) dif = "dif_1_1_2";
    if(islab == 2) dif = "dif_1_1_3";
    if(islab == 3) dif = "dif_1_1_4";
    if(islab == 4) dif = "dif_1_1_5";

    TString filename = fullpath;
    filename += dif;
    filename += ".raw.root";

    cout << "merging " << filename << endl;
    f[islab] = new TFile(filename,"read");
    if(f[islab]->IsOpen()){

      tree[islab] = (TTree*)f[islab]->Get("fev10");

      tree[islab]->SetBranchAddress("acqNumber", &acqNumber_in[islab]);
      tree[islab]->SetBranchAddress("chipid", chipid_in[islab]);
      tree[islab]->SetBranchAddress("nColumns", numCol_in[islab]);
      tree[islab]->SetBranchAddress("nhits", nhits_in[islab]);
      tree[islab]->SetBranchAddress("bcid", bcid_in[islab]);
      tree[islab]->SetBranchAddress("corrected_bcid", corrected_bcid_in[islab]);
      tree[islab]->SetBranchAddress("badbcid", badbcid_in[islab]);
      tree[islab]->SetBranchAddress("charge_lowGain", charge_lowGain_in[islab]);
      tree[islab]->SetBranchAddress("charge_hiGain", charge_hiGain_in[islab]);
      tree[islab]->SetBranchAddress("gain_hit_low", gain_hit_low_in[islab]);
      tree[islab]->SetBranchAddress("gain_hit_high", gain_hit_high_in[islab]);

      Entries[islab] = tree[islab]->GetEntries();
      if(Entries[islab]>MaxEntries) MaxEntries = Entries[islab];
      }
    }

  // ================================================================================================================================================ //
  // Create New TTree
  TFile *fout = new TFile(fullpath + "merge.root", "recreate");
  TTree *treeout = new TTree("fev10","fev10");
  TString name;

  treeout->Branch("acqNumber", &acqNumber, "acqNumber/I");
  treeout->Branch("chipid", chipid, TString::Format("chipid[%i][%i]/I", MaxSlab, MaxChip));
  treeout->Branch("nColumns", numCol, TString::Format("nColumns[%i][%i]/I", MaxSlab, MaxChip));
  treeout->Branch("bcid", bcid, TString::Format("bcid[%i][%i][%i]/I", MaxSlab, MaxChip, MaxMemory));
  treeout->Branch("corrected_bcid", corrected_bcid, TString::Format("corrected_bcid[%i][%i][%i]/I", MaxSlab, MaxChip, MaxMemory));
  treeout->Branch("badbcid", badbcid, TString::Format("badbcid[%i][%i][%i]/I", MaxSlab, MaxChip, MaxMemory));
  treeout->Branch("nhits", nhits, TString::Format("nhits[%i][%i][%i]/I", MaxSlab, MaxChip, MaxMemory));
  treeout->Branch("charge_lowGain", charge_lowGain, TString::Format("lowGain[%i][%i][%i][%i]/I", MaxSlab, MaxChip, MaxMemory, MaxChannel));
  treeout->Branch("charge_hiGain", charge_hiGain, TString::Format("highGain[%i][%i][%i][%i]/I", MaxSlab, MaxChip, MaxMemory, MaxChannel));
  treeout->Branch("gain_hit_low", gain_hit_low, TString::Format("gain_hit_low[%i][%i][%i][%i]/I", MaxSlab, MaxChip, MaxMemory, MaxChannel));
  treeout->Branch("gain_hit_high", gain_hit_high, TString::Format("gain_hit_high[%i][%i][%i][%i]/I", MaxSlab, MaxChip, MaxMemory, MaxChannel));

  // ================================================================================================================================================ //
  // Fill TTree
  Int_t currentSpill=0;
  bool continueScan = true;
  Int_t MaxEvents = 2000000000;
  Int_t evt = 0;
  bool isSlabEnd[MaxSlab];
  for(Int_t islab=0; islab<MaxSlab; islab++) isSlabEnd[islab] = false;

  while(continueScan && evt<MaxEvents){

    currentSpill = 2000000000;
    for(Int_t islab=0; islab<MaxSlab; islab++){
      if(f[islab]->IsOpen()){
        tree[islab]->GetEntry(Entry[islab]);
        if(currentSpill>acqNumber_in[islab]&&!isSlabEnd[islab]) currentSpill = acqNumber_in[islab];
      }
    }

    for(Int_t islab=0; islab<MaxSlab; islab++){
      if(currentSpill==acqNumber_in[islab]&&Entries[islab]>Entry[islab]){
        Entry[islab]++;
        for(Int_t ichip=0; ichip<MaxChip; ichip++){
          Int_t imemory = 0;
          Int_t loopBCID = 0;
          for(Int_t iraw=0; iraw<MaxMemory; iraw++){
            if(filterEvents){
              if(iraw==0) imemory = iraw;
              else if(bcid_in[islab][ichip][iraw] - bcid_in[islab][ichip][iraw-1]>5) imemory = iraw;
              else if(bcid_in[islab][ichip][iraw] - bcid_in[islab][ichip][iraw-1]<0 && bcid_in[islab][ichip][iraw]+4096 - bcid_in[islab][ichip][iraw-1]>5){
                imemory = iraw;
                loopBCID++;
              }
              else continue;
	    }
            else imemory = iraw;
            bcid[islab][ichip][imemory] = bcid_in[islab][ichip][imemory];
            badbcid[islab][ichip][imemory] = badbcid_in[islab][ichip][imemory];
            corrected_bcid[islab][ichip][imemory] = corrected_bcid_in[islab][ichip][imemory];
            nhits[islab][ichip][imemory] = nhits_in[islab][ichip][imemory];
            for(Int_t ichn=0; ichn<MaxChannel; ichn++){
              charge_lowGain[islab][ichip][imemory][ichn] = charge_lowGain_in[islab][ichip][imemory][ichn];
              charge_hiGain[islab][ichip][imemory][ichn] = charge_hiGain_in[islab][ichip][imemory][ichn];
              gain_hit_low[islab][ichip][imemory][ichn] = gain_hit_low_in[islab][ichip][imemory][ichn];
              gain_hit_high[islab][ichip][imemory][ichn] = gain_hit_high_in[islab][ichip][imemory][ichn];
	    }
	  }
          chipid[islab][ichip] = chipid_in[islab][ichip];
          numCol[islab][ichip] = numCol_in[islab][ichip];
	}
        acqNumber = currentSpill;
      }
      else{
	//treeInitslab(islab);
        for(Int_t ichip=0; ichip<MaxChip; ichip++){
          chipid[islab][ichip] = -999;
          numCol[islab][ichip] = -999;
          for(Int_t imemory=0; imemory<MaxMemory; imemory++){
            bcid[islab][ichip][imemory] = -999;
            badbcid[islab][ichip][imemory] = 0;
            corrected_bcid[islab][ichip][imemory] = -999;
            nhits[islab][ichip][imemory] = -999;
            for(Int_t ichn=0; ichn<MaxChannel; ichn++){
              charge_lowGain[islab][ichip][imemory][ichn] = -999;
              charge_hiGain[islab][ichip][imemory][ichn] = -999;
              gain_hit_low[islab][ichip][imemory][ichn] = -999;
              gain_hit_high[islab][ichip][imemory][ichn] = -999;
	    }
	  }
	}
      }
    }
    treeout->Fill();
    evt++;
    continueScan = false;
    for(Int_t islab=0; islab<MaxSlab; islab++){
      isSlabEnd[islab] = true;
      if(Entries[islab]>Entry[islab]){
        continueScan = true;
        isSlabEnd[islab]=false;
      }
    }
  }
  fout->cd();
  fout->Write(0);
  fout->Close();
}

