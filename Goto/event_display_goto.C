#include <stdlib.h>
#include <fstream>
#include <math.h>
#include "TMath.h"
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
#include "TH3F.h"
#include "TSpectrum.h"

using namespace::std;

void event_display_goto(std::string str){
  // Cut File name
  int path_i = str.find_last_of("/")+1;
  int ext_i = str.find_last_of(".");
  std::string pathname = str.substr(0, path_i);
  std::string filename_ = str.substr(path_i, ext_i-path_i);
  //std::string extroot = filename_.substr(ext_i, str.size()-ext_i);

  TString fullpath = str.data();
  TString filename = pathname + filename_;

  // ================================================================================================================================================ //
  // TTree open  
  TFile *file = TFile::Open(fullpath); 
  TTree *Event_Tree = (TTree*)file->Get("Event_Tree");

  // ================================================================================================================================================ //
  // Set Loop Max
  Int_t MaxEvent = Event_Tree->GetEntries();
  Int_t NOver5Hits = Event_Tree->GetEntries("nhit_channel>5");
  Int_t MaxSlab = 5;
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxMemory = 15;
  Int_t NMAX = MaxSlab*MaxChip*MaxChannel*MaxSlab;

  // ================================================================================================================================================ //
  // Reading ADC data from Tree
  Int_t acqNumber;
  Event_Tree->SetBranchAddress("acqNumber", &acqNumber);

  Int_t nhit_channel;
  Event_Tree->SetBranchAddress("nhit_channel", &nhit_channel);

  Int_t bcid_hit;
  Event_Tree->SetBranchAddress("bcid_hit", &bcid_hit);

  Int_t slab_number[NMAX];
  Event_Tree->SetBranchAddress("slab_number", slab_number);

  Int_t chip_number[NMAX];
  Event_Tree->SetBranchAddress("chip_number", chip_number);

  Int_t sca_number[NMAX];
  Event_Tree->SetBranchAddress("sca_number", sca_number);

  Int_t channel_number[NMAX];
  Event_Tree->SetBranchAddress("channel_number", channel_number);

  Float_t hit_x[NMAX];
  Event_Tree->SetBranchAddress("hit_x", hit_x);

  Float_t hit_y[NMAX];
  Event_Tree->SetBranchAddress("hit_y", hit_y);

  Float_t hit_z[NMAX];
  Event_Tree->SetBranchAddress("hit_z", hit_z);

  Int_t charge_hiGain_tdc[NMAX];
  Event_Tree->SetBranchAddress("charge_hiGain_tdc", charge_hiGain_tdc);

  Int_t charge_lowGain_adc[NMAX];
  Event_Tree->SetBranchAddress("charge_lowGain_adc", charge_lowGain_adc);

  // ================================================================================================================================================ //
  // Make Event Display
  TH3F* mip_evdisp[NOver5Hits];
  for(Int_t i=0; i<NOver5Hits; i++) mip_evdisp[i]= new TH3F(TString::Format("mip_evdisp_%i",i),TString::Format("mip_evdisp_%i",i),32,-90,90,20,-0.5,19.5,32,-90,90);
  // ================================================================================================================================================ //
  // Fill Event Display
  Int_t cuthits = 0;
  for (Int_t ievent=0; ievent<MaxEvent;ievent++){
    Event_Tree->GetEntry(ievent);
    if(float(ievent%2000)==0) cout << "Event" << ievent << "End" << endl;

    if(nhit_channel>5){

      // =============================================================================================== //
      // Cut condition loop
      Int_t nslabs[5] = {0};
      Int_t nchips[16] = {0};
      Int_t nchannels[64] = {0};
      Int_t nadc_charge = 0;
      bool bslab = true, bchip = false, bchannel = false;
      for(Int_t ihit=0; ihit<nhit_channel; ihit++){
	//if(slab_number[ihit]>MaxSlab) continue;
	//if(chip_number[ihit]>MaxChip) continue;
	nslabs[slab_number[ihit]]++;
	nchips[chip_number[ihit]]++;
	nchannels[channel_number[ihit]]++;
	if(charge_lowGain_adc[ihit]<200&&charge_lowGain_adc[ihit]>100) nadc_charge++;
      }
      for(Int_t islab=0; islab<MaxSlab; islab++){
        if(nslabs[islab]==0) bslab = false;
      }
      for(Int_t ichip=0; ichip<MaxChip; ichip++){
        if(nchips[ichip]>4){
	  bchip = true;
	  //cout << "break chip cut loop!" << endl;
	  break;
	}
      }
      for(Int_t ichn=0; ichn<MaxChannel; ichn++){
        if(nchannels[ichn]>4){
	  bchannel = true;
	  //cout << "break channel cut loop!" << endl;
	  break;
	}
      }
      // Cut slab condition
      if(bslab==false) continue;
      // Cut chip condition
      if(bchip==false) continue;
      // Cut channel condition
      if(bchannel==false) continue;
      // Cut charge condition
      if(nadc_charge<5) continue;

      // =============================================================================================== //
      for(Int_t ihit=0; ihit<nhit_channel; ihit++){
	mip_evdisp[cuthits]->Fill(hit_x[ihit], hit_z[ihit], hit_y[ihit], charge_lowGain_adc[ihit]);
      }
      cuthits++;
    }

  }
  // ================================================================================================================================================ //
  // Write Event Display
  TFile *fout = new TFile(filename+"_EvDisplay.root" , "RECREATE");
  fout->cd();
  for(int i=0; i<cuthits; i++) {
    mip_evdisp[i]->GetXaxis()->SetTitle("X");
    mip_evdisp[i]->GetYaxis()->SetTitle("Z");
    mip_evdisp[i]->GetZaxis()->SetTitle("Y");
    mip_evdisp[i]->Write();
  }
  cout << "Event" << cuthits << "Created" << endl;
  fout->Close();
}
