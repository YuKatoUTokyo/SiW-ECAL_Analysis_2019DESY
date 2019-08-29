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
#include "TSpectrum.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"

using namespace::std;

void ped_time_chip(){

  Int_t MaxSlab = 5;
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxMemory = 15;
  Int_t chip_mapper[64] = {16,  7, 19, 25, 31, 35, 38, 47,
                           10, 22, 13, 28, 32, 41, 44, 53,
  		            9,  3, 15, 27, 33, 45, 50, 60,
  		           21, 12, 18, 30, 36, 48, 54, 39,
  		            6,  5,  4, 20, 40, 42, 58, 56,
  		            2,  8, 11, 23, 37, 49, 63, 57,
		            1, 24, 17, 29, 43, 52, 55, 51,
                            0, 14, 26, 34, 46, 59, 62, 61};

  TGraph *ped_var[5][16][15][64];

  for(Int_t islab=0; islab<MaxSlab; islab++){
    for(Int_t ichip=0; ichip<MaxChip; ichip++){
      for(Int_t imemory=0; imemory<MaxMemory; imemory++){
        for(Int_t ichn=0; ichn<MaxChannel; ichn++){
	  ped_var[islab][ichip][imemory][ichn] = new TGraph(0);
	}
      }
    }
  }

  // Read Time
  TString time_path = "./run_second.txt";
  std::ifstream reading_time(time_path, std::ios::in);
  if(!reading_time){
    cout << "run_second is not found" << endl;
  }
  TString tmp_run;
  Float_t tmp_first = 0, tmp_final = 0, tmp_acqtime = 0;
  Int_t MaxRun[5][15][64] = {0};
  std::vector<std::vector< std::vector<std::vector<Float_t> > > > time_label;
  std::vector<std::vector< std::vector<std::vector<Int_t> > > > xdata;
  std::vector<std::vector< std::vector<std::vector<TString> > > > run_label;
  for(Int_t islab=0; islab<MaxSlab; islab++){
    std::vector<std::vector<std::vector<Float_t> > > temp_time;
    std::vector<std::vector<std::vector<Int_t> > > temp_xdata;
    std::vector<std::vector<std::vector<TString> > > temp_run;
    for(Int_t imemory=0; imemory<MaxMemory; imemory++){
      std::vector<std::vector<Float_t> > temp_time2;
      std::vector<std::vector<Int_t> > temp_xdata2;
      std::vector<std::vector<TString> > temp_run2;
      for(Int_t ichn=0; ichn<MaxChannel; ichn++){
	std::vector<Float_t> time_label2;
	std::vector<Int_t> xdata2;
	std::vector<TString> run_label2;
	temp_time2.push_back(time_label2);
	temp_xdata2.push_back(xdata2);
	temp_run2.push_back(run_label2);
      }
      temp_time.push_back(temp_time2);
      temp_xdata.push_back(temp_xdata2);
      temp_run.push_back(temp_run2);
    }
    time_label.push_back(temp_time);
    xdata.push_back(temp_xdata);
    run_label.push_back(temp_run);
  }

  while(reading_time){
    reading_time >> tmp_run >> tmp_first >> tmp_final >> tmp_acqtime;

    // TTree open
    TString filename = "./Chip_Result/" + tmp_run + "_merge_PedestalMap.root";
    TFile *file = TFile::Open(filename); 
    if(file){
    
      TTree *Pedestal_Tree = (TTree*)file->Get("Pedestal_Tree");

      Int_t MaxEvent = Pedestal_Tree->GetEntries();
      Float_t runtime = tmp_first + tmp_acqtime/2;

      // Read data from Pedestal_Tree
      Int_t most_hit_chips[5];
      Pedestal_Tree->SetBranchAddress("most_hit_chips", most_hit_chips);
      Double_t pedestal_mean[5][15][64];
      Pedestal_Tree->SetBranchAddress("pedestal_mean", pedestal_mean);

      for(Int_t ievent=0; ievent<MaxEvent; ievent++){
        Pedestal_Tree->GetEntry(ievent);
        for(Int_t islab=0; islab<MaxSlab; islab++){
	  Int_t ichip = most_hit_chips[islab];
	  cout <<  "most hit chip " << most_hit_chips[0] << endl;
          for(Int_t imemory=0; imemory<MaxMemory; imemory++){
	    for(Int_t ichn=0; ichn<MaxChannel; ichn++){
              MaxRun[islab][imemory][ichn]++;
              xdata.at(islab).at(imemory).at(ichn).push_back(MaxRun[islab][imemory][ichn]);
              time_label.at(islab).at(imemory).at(ichn).push_back(runtime);
              run_label.at(islab).at(imemory).at(ichn).push_back(tmp_run);

	      ped_var[islab][ichip][imemory][ichn]->SetPoint(ped_var[islab][ichip][imemory][ichn]->GetN(), runtime, pedestal_mean[islab][imemory][ichn]);
	    }
	  }
        }
      }
    }
  }
 
  Int_t imemory = 0;
  Int_t Color[5] = {1, 2, 3, 4, 6};
  TCanvas *canvas_pedestal_map[16];
  TLegend *legend[16][64];
  TText *t[5][15][64];
  for(Int_t ichip=0; ichip<MaxChip; ichip++){
    canvas_pedestal_map[ichip]= new TCanvas(TString::Format("pedestal_map_chip%i",ichip),
		                            TString::Format("pedestal_map_chip%i",ichip), 1200, 800);
    canvas_pedestal_map[ichip]->Print(TString::Format("./Time_Dependance/pedestal_vs_time_chip%i.pdf[", ichip));
    for(Int_t ichn=0; ichn<MaxChannel; ichn++){
      canvas_pedestal_map[ichip]->DrawFrame(0, 200, 835000, 400);
      legend[ichip][ichn] = new TLegend( 0.6, 0.7, 0.99, 0.90);
      for(Int_t islab=0; islab<MaxSlab; islab++){
        ped_var[islab][ichip][imemory][ichn]->GetXaxis()->SetTitle("Time[second]");
        ped_var[islab][ichip][imemory][ichn]->GetYaxis()->SetTitle("Pedestal Mean");
        //ped_var[islab][ichip][imemory][ichn]->GetXaxis()->SetRangeUser(0, 835000);
        //ped_var[islab][ichip][imemory][ichn]->GetYaxis()->SetRangeUser(200, 400);
        //ped_var[islab][ichip][imemory][ichn]->SetLineColor(Color[islab]);
        ped_var[islab][ichip][imemory][ichn]->SetMarkerColor(Color[islab]);
	if(islab==0) ped_var[islab][ichip][imemory][ichn]->Draw("Psame");
	else ped_var[islab][ichip][imemory][ichn]->Draw("Psame");
        t[islab][imemory][ichn] = new TText();
        t[islab][imemory][ichn]->SetTextAngle(0);
        t[islab][imemory][ichn]->SetTextSize(0.04);
        t[islab][imemory][ichn]->SetTextAlign(-10);
        for(Int_t irun=0; irun<MaxRun[islab][imemory][ichn]; irun++){
          t[islab][imemory][ichn]->DrawText(-0.42, time_label[islab][imemory][ichn][irun], run_label[islab][imemory][ichn][irun].Data());
	}
        legend[ichip][ichn]->AddEntry(ped_var[islab][ichip][imemory][ichn], TString::Format("ped_slab%i_chip%i_chn%i_sca%i", islab, ichip, ichn, imemory), "p");
	gPad->Modified();
        gPad->Update();
      }
      legend[ichip][ichn]->Draw("same");
      canvas_pedestal_map[ichip]->RedrawAxis();
      canvas_pedestal_map[ichip]->Print(TString::Format("./Time_Dependance/pedestal_vs_time_chip%i.pdf", ichip));
    }
    canvas_pedestal_map[ichip]->Print(TString::Format("./Time_Dependance/pedestal_vs_time_chip%i.pdf]", ichip));
  }

}


