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

using namespace::std;

void ped_analysis_tdc(std::string str){

  // Cut File name
  int path_i = str.find_last_of("/")+1;
  int ext_i = str.find_last_of(".");
  std::string pathname = str.substr(0, path_i);
  std::string extroot = str.substr(ext_i, str.size()-ext_i);
  std::string filename_i = str.substr(path_i, ext_i-path_i);
  int ext_j = filename_i.find_last_of(".");
  std::string extraw = filename_i.substr(ext_j, filename_i.size()-ext_j);
  std::string true_filename = filename_i.substr(0, ext_j);

  TString fullpath = str.data();
  TString filename = pathname + true_filename;

  // TTree open  
  TFile *file = TFile::Open(fullpath); 
  TTree *fev10 = (TTree*)file->Get("fev10");

  // Read Chip Map
  TString mapchip_path = "/home/goto/ILC/SiWECAL_2019/Analysis/Pedestal_Stability/Macro/map_chip.dat";
  std::ifstream reading_file(mapchip_path, std::ios::in);
  if(!reading_file){
    cout << "map_chip.dat is not found" << endl;
  }
  
  Float_t map_pointX[16][64], map_pointY[16][64];
  Int_t tmp_chip = 0, tmp_channel = 0;
  Float_t tmp_x0 = 0, tmp_y0 = 0, tmp_x = 0, tmp_y = 0;
  while(reading_file){
    reading_file >> tmp_chip >> tmp_x0 >> tmp_y0 >> tmp_channel >> tmp_x >> tmp_y;
    map_pointX[tmp_chip][tmp_channel] = -tmp_x;
    map_pointY[tmp_chip][tmp_channel] = tmp_y;
  }
  
  Int_t MaxEvent = fev10->GetEntries();
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxSca = 1;

  TH2F* pedestal_map[15];
  TH2F* pedestal_width_map[15];
  TH2F* pedestal_map_badbcid_0[15];
  TH2F* pedestal_width_map_badbcid_0[15];
  TH2F* pedestal_map_badbcid_not_0[15];
  TH2F* pedestal_width_map_badbcid_not_0[15];

  for(int isca=0; isca<MaxSca; isca++){
    pedestal_map[isca]= new TH2F(TString::Format("pedestal_map_sca%i",isca),
		                 TString::Format("pedestal_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_map_badbcid_0[isca]= new TH2F(TString::Format("pedestal_map_badbcid_0_sca%i",isca),
		    			   TString::Format("pedestal_map_badbcid_0_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_map_badbcid_not_0[isca]= new TH2F(TString::Format("pedestal_map_badbcid_not_0_sca%i",isca),
		    			       TString::Format("pedestal_map_badbcid_not_0_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_width_map[isca]= new TH2F(TString::Format("pedestal_width_map_sca%i",isca),
		                       TString::Format("pedestal_width_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_width_map_badbcid_0[isca]= new TH2F(TString::Format("pedestal_width_map_badbcid_0_sca%i",isca),
		                                 TString::Format("pedestal_width_map_badbcid_0_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
    pedestal_width_map_badbcid_not_0[isca]= new TH2F(TString::Format("pedestal_width_map_badbcid_not_0_sca%i",isca),
		                                     TString::Format("pedestal_width_map_badbcid_not_0_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
  }

  // Read data from fev10
  Int_t charge_lowGain[16][15][64];
  fev10->SetBranchAddress("charge_lowGain", charge_lowGain);

  Int_t gain_hit_low[16][15][64];
  fev10->SetBranchAddress("gain_hit_low", gain_hit_low);

  Int_t badbcid[16][15];
  fev10->SetBranchAddress("badbcid", badbcid);

  // Fill to Tree Data
  Double_t pedestal_mean[16][15][64] = {-1};
  Double_t pedestal_width[16][15][64] = {-1};
  Double_t pedestal_error[16][15][64] = {-1};
  Double_t pedestal_chi2ndf[16][15][64] = {-1};
  Double_t pedestal_mean_badbcid_0[16][15][64] = {-1};
  Double_t pedestal_width_badbcid_0[16][15][64] = {-1};
  Double_t pedestal_mean_badbcid_not_0[16][15][64] = {-1};
  Double_t pedestal_width_badbcid_not_0[16][15][64] = {-1};

  // Create New TTree
  TFile *fout = new TFile("/home/goto/ILC/SiWECAL_2019/Analysis/Pedestal_Stability/Pedestal_Map/" + filename + "_PedestalMap.root" , "RECREATE");
  TTree *Pedestal_Tree = new TTree("Pedestal_Tree", "Pedestal_Tree");
  Pedestal_Tree->Branch("pedestal_mean", pedestal_mean,
		        "pedestal_mean[16][15][64]/D");
  Pedestal_Tree->Branch("pedestal_width", pedestal_width,
		        "pedestal_width[16][15][64]/D");
  Pedestal_Tree->Branch("pedestal_error", pedestal_error,
		        "pedestal_error[16][15][64]/D");
  Pedestal_Tree->Branch("pedestal_chi2ndf", pedestal_chi2ndf,
		        "pedestal_chi2ndf[16][15][64]/D");
  Pedestal_Tree->Branch("pedestal_mean_badbcid_0", pedestal_mean_badbcid_0,
		        "pedestal_mean_badbcid_0[16][15][64]/D");
  Pedestal_Tree->Branch("pedestal_width_badbcid_0", pedestal_width_badbcid_0,
		        "pedestal_width_badbcid_0[16][15][64]/D");
  Pedestal_Tree->Branch("pedestal_mean_badbcid_not_0", pedestal_mean_badbcid_not_0,
		        "pedestal_mean_badbcid_not_0[16][15][64]/D");
  Pedestal_Tree->Branch("pedestal_width_badbcid_not_0", pedestal_width_badbcid_not_0,
		        "pedestal_width_badbcid_not_0[16][15][64]/D");

  TDirectory *hist_sca[16];
  Double_t hit_number[16][64] = {0};

 
  // --------------------
  // all sca
  std::vector<std::vector<std::vector<TH1F*> > > ped_sca;
  std::vector<std::vector<std::vector<TH1F*> > > ped_sca_badbcid_0;
  std::vector<std::vector<std::vector<TH1F*> > > ped_sca_badbcid_not_0;

  std::vector<TH1F*> pedestal_chip ;
  std::vector<TH1F*> pedestal_diff_chip ;

  for(int ichip=0; ichip<MaxChip; ichip++){
    TH1F *ped_chip = new TH1F(TString::Format("ped_chip%i", ichip), 
		    	      TString::Format("ped_chip%i", ichip), 1000, 0.5, 1000.5);
    pedestal_chip.push_back(ped_chip);

    TH1F *ped_diff_chip = new TH1F(TString::Format("ped_diff_chip%i", ichip),
		                   TString::Format("ped_diff_chip%i", ichip), 1002, -500, 500);
    pedestal_diff_chip.push_back(ped_diff_chip);
    
    std::vector<std::vector<TH1F*> >pedtemp_sca;
    std::vector<std::vector<TH1F*> >pedtemp_sca_badbcid_0;
    std::vector<std::vector<TH1F*> >pedtemp_sca_badbcid_not_0;

    for(int ichn=0; ichn<MaxChannel; ichn++){
      std::vector<TH1F*> pedtemp_sca2;
      std::vector<TH1F*> pedtemp_sca2_badbcid_0;
      std::vector<TH1F*> pedtemp_sca2_badbcid_not_0;

      for(int isca=0; isca<MaxSca; isca++){
	TH1F *ped_sca2 = new TH1F(TString::Format("ped_chip%i_chn%i_sca%i", ichip, ichn, isca),
			          TString::Format("ped_chip%i_chn%i_sca%i", ichip, ichn, isca), 1000, 0.5, 1000.5);
	TH1F *ped_sca2_badbcid_0 = new TH1F(TString::Format("ped_badbcid_0_chip%i_chn%i_sca%i", ichip, ichn, isca),
			                    TString::Format("ped_badbcid_0_chip%i_chn%i_sca%i", ichip, ichn, isca), 1000, 0.5, 1000.5);
	TH1F *ped_sca2_badbcid_not_0 = new TH1F(TString::Format("ped_badbcid_not_0_chip%i_chn%i_sca%i", ichip, ichn, isca),
			                        TString::Format("ped_badbcid_not_0_chip%i_chn%i_sca%i", ichip, ichn, isca), 1000, 0.5, 1000.5);
	pedtemp_sca2.push_back(ped_sca2);
	pedtemp_sca2_badbcid_0.push_back(ped_sca2_badbcid_0);
	pedtemp_sca2_badbcid_not_0.push_back(ped_sca2_badbcid_not_0);
      }
      pedtemp_sca.push_back(pedtemp_sca2);
      pedtemp_sca_badbcid_0.push_back(pedtemp_sca2_badbcid_0);
      pedtemp_sca_badbcid_not_0.push_back(pedtemp_sca2_badbcid_not_0);
    }
    ped_sca.push_back(pedtemp_sca);
    ped_sca_badbcid_0.push_back(pedtemp_sca_badbcid_0);
    ped_sca_badbcid_not_0.push_back(pedtemp_sca_badbcid_not_0);
  }

  // -----------------------------------------------------------------------------------------------------   
  // SCA analysis
  for (int ievent=0; ievent<MaxEvent; ievent++){
    fev10->GetEntry(ievent);

    for(int ichip=0; ichip<MaxChip; ichip++){

      for(int isca=0; isca<MaxSca; isca++){
        if(badbcid[ichip][isca]!=0){
	  for(int ichn=0; ichn<MaxChannel; ichn++){
	    if(charge_lowGain[ichip][isca][ichn]>10 && gain_hit_low[ichip][isca][ichn]==0){
	      ped_sca.at(ichip).at(ichn).at(isca)->Fill(charge_lowGain[ichip][isca][ichn]);
	      ped_sca_badbcid_not_0.at(ichip).at(ichn).at(isca)->Fill(charge_lowGain[ichip][isca][ichn]);
	    }
	    else if(gain_hit_low[ichip][isca][ichn]==1) hit_number[ichip][ichn] += gain_hit_low[ichip][isca][ichn];
	  }
	}
	else if(badbcid[ichip][isca]==0){
	  for(int ichn=0; ichn<MaxChannel; ichn++){
	    if(charge_lowGain[ichip][isca][ichn]>10 && gain_hit_low[ichip][isca][ichn]==0){
	      ped_sca.at(ichip).at(ichn).at(isca)->Fill(charge_lowGain[ichip][isca][ichn]);
	      ped_sca_badbcid_0.at(ichip).at(ichn).at(isca)->Fill(charge_lowGain[ichip][isca][ichn]);
	    }
	    else if(gain_hit_low[ichip][isca][ichn]==1) hit_number[ichip][ichn] += gain_hit_low[ichip][isca][ichn];
	  }
	}

      }//isca

    }//ichip 
   
  }  // end first loop analysis to fill pedestal historgrams

  std::vector<std::vector<std::vector<Double_t> > > ped_mean;
  std::vector<std::vector<std::vector<Double_t> > > ped_error;
  std::vector<std::vector<std::vector<Double_t> > > ped_width;
  std::vector<std::vector<std::vector<Double_t> > > ped_mean_badbcid_0;
  std::vector<std::vector<std::vector<Double_t> > > ped_width_badbcid_0;
  std::vector<std::vector<std::vector<Double_t> > > ped_mean_badbcid_not_0;
  std::vector<std::vector<std::vector<Double_t> > > ped_width_badbcid_not_0;
  Double_t ped_mean_integral_badbcid_0[16][15][64] = {0};
  Double_t ped_mean_integral_badbcid_not_0[16][15][64] = {0};

  //initialize pedestal vectors
  for(int i=0; i<16; i++){
    std::vector<std::vector<Double_t> > chip_ped_mean;
    std::vector<std::vector<Double_t> > chip_ped_error;
    std::vector<std::vector<Double_t> > chip_ped_width;

    for(int j=0; j<64; j++){
      std::vector<Double_t> chn_ped_mean;
      std::vector<Double_t> chn_ped_error;
      std::vector<Double_t> chn_ped_width;
      chip_ped_mean.push_back(chn_ped_mean);
      chip_ped_error.push_back(chn_ped_error);
      chip_ped_width.push_back(chn_ped_width);
    }  
    ped_mean.push_back(chip_ped_mean);
    ped_error.push_back(chip_ped_error);
    ped_width.push_back(chip_ped_width);
    ped_mean_badbcid_0.push_back(chip_ped_mean);
    ped_width_badbcid_0.push_back(chip_ped_width);
    ped_mean_badbcid_not_0.push_back(chip_ped_mean);
    ped_width_badbcid_not_0.push_back(chip_ped_width);
  }

  // do pedestal (chip/channel/sca based) analysis
  for(int ichip=0; ichip<MaxChip; ichip++){
    cout << "CHIP " << ichip << " START" << endl;
    TString TDirectory_Name = Form("ALL_PEDESTAL_%02dCHIP", ichip);
    hist_sca[ichip] = fout->mkdir(TDirectory_Name);
    hist_sca[ichip]->cd();
    for(int ichn=0; ichn<MaxChannel; ichn++){
      for(int isca=0; isca<MaxSca; isca++){

	ped_sca.at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca));
	ped_sca.at(ichip).at(ichn).at(isca)->SetName(TString::Format("ped_chip%i_chn%i_sca%i",ichip,ichn,isca));

	ped_sca_badbcid_0.at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("ped_badbcid_0_chip%i_chn%i_sca%i",ichip,ichn,isca));
	ped_sca_badbcid_0.at(ichip).at(ichn).at(isca)->SetName(TString::Format("ped_badbcid_0_chip%i_chn%i_sca%i",ichip,ichn,isca));

	ped_sca_badbcid_not_0.at(ichip).at(ichn).at(isca)->SetTitle(TString::Format("ped_badbcid_not_0_chip%i_chn%i_sca%i",ichip,ichn,isca));
	ped_sca_badbcid_not_0.at(ichip).at(ichn).at(isca)->SetName(TString::Format("ped_badbcid_not_0_chip%i_chn%i_sca%i",ichip,ichn,isca));

	ped_mean.at(ichip).at(ichn).push_back(0.);
	ped_width.at(ichip).at(ichn).push_back(0.);

	ped_mean_badbcid_0.at(ichip).at(ichn).push_back(0.);
	ped_width_badbcid_0.at(ichip).at(ichn).push_back(0.);

	ped_mean_badbcid_not_0.at(ichip).at(ichn).push_back(0.);
	ped_width_badbcid_not_0.at(ichip).at(ichn).push_back(0.);

	//====================================================================================================
	if(ped_sca.at(ichip).at(ichn).at(isca)->GetEntries()>10){ //max_entries/2 ){
	  TSpectrum *s = new TSpectrum();
	  int npeaks = s->Search(ped_sca.at(ichip).at(ichn).at(isca), 2, "", 0.2); 

	  if(npeaks > 0){

            Double_t *mean_peak=s->GetPositionX();
            Double_t *mean_high=s->GetPositionY();
            double mean_peak_higher=0;
            double mean_high_higher=0;
	    int npeak_max=0;

            for(int ipeak=0; ipeak<npeaks; ipeak++){
              if(mean_high[ipeak]>mean_high_higher && mean_high[ipeak]>50){
                mean_high_higher=mean_high[ipeak];
                mean_peak_higher=mean_peak[ipeak];
		npeak_max=ipeak;
              }
            }

            for(int ipeak=0; ipeak<npeaks; ipeak++){
	      if(ipeak != npeak_max) pedestal_diff_chip.at(ichip)->Fill(mean_peak[npeak_max] - mean_peak[ipeak]);
	    }

	  }


	  TH1 *pedestal_histogram = ped_sca.at(ichip).at(ichn).at(isca);
	  Double_t mean_peak = pedestal_histogram->GetBinCenter(pedestal_histogram->GetMaximumBin());
	    
	  TF1 *f0 = new TF1("f0", "gaus", mean_peak-2*ped_sca.at(ichip).at(ichn).at(isca)->GetRMS(), 
			      mean_peak+2*ped_sca.at(ichip).at(ichn).at(isca)->GetRMS());
	  ped_sca.at(ichip).at(ichn).at(isca)->Fit("f0", "RQNOC");

	  TF1 *f1 = new TF1("f1", "gaus", f0->GetParameter(1)-2.*f0->GetParameter(2), 
			      f0->GetParameter(1)+2.*f0->GetParameter(2));
	  ped_sca.at(ichip).at(ichn).at(isca)->Fit("f1","RQME");

  	  gStyle->SetOptFit(1111);
	  ped_sca.at(ichip).at(ichn).at(isca)->Write();
	  
	  ped_mean.at(ichip).at(ichn).at(isca)=f1->GetParameter(1);
	  ped_width.at(ichip).at(ichn).at(isca)=f1->GetParameter(2);

	  pedestal_map[isca]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(1));
	  pedestal_width_map[isca]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f1->GetParameter(2));
	      
	  pedestal_mean[ichip][isca][ichn] = f1->GetParameter(1);
	  pedestal_width[ichip][isca][ichn] = f1->GetParameter(2);

	}

	//====================================================================================================
	if(ped_sca_badbcid_0.at(ichip).at(ichn).at(isca)->GetEntries()>10){ //max_entries/2 ){
	  TSpectrum *s2 = new TSpectrum();
	  int npeaks = s2->Search(ped_sca_badbcid_0.at(ichip).at(ichn).at(isca), 2, "", 0.2); 

	  if(npeaks > 0){

            Double_t *mean_peak2=s2->GetPositionX();
            Double_t *mean_high2=s2->GetPositionY();
            double mean_peak_higher2=0;
            double mean_high_higher2=0;
	    int npeak_max=0;

	  }

	  TH1 *pedestal_histogram_badbcid_0 = ped_sca_badbcid_0.at(ichip).at(ichn).at(isca);
	  Double_t mean_peak2 = pedestal_histogram_badbcid_0->GetBinCenter(pedestal_histogram_badbcid_0->GetMaximumBin());
	    
	  TF1 *f12 = new TF1("f12", "gaus",
	  		     pedestal_mean[ichip][isca][ichn]-2.*pedestal_width[ichip][isca][ichn],
	  		     pedestal_mean[ichip][isca][ichn]+2.*pedestal_width[ichip][isca][ichn]);
	  ped_sca_badbcid_0.at(ichip).at(ichn).at(isca)->Fit("f12","RQME");

  	  gStyle->SetOptFit(1111);
	  ped_sca_badbcid_0.at(ichip).at(ichn).at(isca)->Write();

	  ped_mean_badbcid_0.at(ichip).at(ichn).at(isca)=f12->GetParameter(1);
	  ped_width_badbcid_0.at(ichip).at(ichn).at(isca)=f12->GetParameter(2);

	  pedestal_map_badbcid_0[isca]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f12->GetParameter(1));
	  pedestal_width_map_badbcid_0[isca]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f12->GetParameter(2));
	      
	  pedestal_mean_badbcid_0[ichip][isca][ichn] = f12->GetParameter(1);
	  pedestal_width_badbcid_0[ichip][isca][ichn] = f12->GetParameter(2);

          ped_mean_integral_badbcid_0[ichip][isca][ichn] = ped_sca_badbcid_0.at(ichip).at(ichn).at(isca)->Integral(
	  		       int(pedestal_mean_badbcid_0[ichip][isca][ichn]-5.*pedestal_width_badbcid_0[ichip][isca][ichn]),
	  		       int(pedestal_mean_badbcid_0[ichip][isca][ichn]+5.*pedestal_width_badbcid_0[ichip][isca][ichn]));
	
	}

	//====================================================================================================
	if(ped_sca_badbcid_not_0.at(ichip).at(ichn).at(isca)->GetEntries()>10){ //max_entries/2 ){
	  TSpectrum *s3 = new TSpectrum();
	  int npeaks = s3->Search(ped_sca_badbcid_not_0.at(ichip).at(ichn).at(isca), 2, "", 0.2); 

	  if(npeaks > 0){

            Double_t *mean_peak3=s3->GetPositionX();
            Double_t *mean_high3=s3->GetPositionY();
            double mean_peak_higher3=0;
            double mean_high_higher3=0;
	    int npeak_max=0;

	  }

	  TH1 *pedestal_histogram_badbcid_not_0 = ped_sca_badbcid_not_0.at(ichip).at(ichn).at(isca);
	  Double_t mean_peak3 = pedestal_histogram_badbcid_not_0->GetBinCenter(pedestal_histogram_badbcid_not_0->GetMaximumBin());
	    
	  TF1 *f13 = new TF1("f13", "gaus",
	  		     pedestal_mean[ichip][isca][ichn]-2.*pedestal_width[ichip][isca][ichn],
	  		     pedestal_mean[ichip][isca][ichn]+2.*pedestal_width[ichip][isca][ichn]);
	  ped_sca_badbcid_not_0.at(ichip).at(ichn).at(isca)->Fit("f13","RQME");

  	  gStyle->SetOptFit(1111);
	  ped_sca_badbcid_not_0.at(ichip).at(ichn).at(isca)->Write();

	  ped_mean_badbcid_not_0.at(ichip).at(ichn).at(isca)=f13->GetParameter(1);
	  ped_width_badbcid_not_0.at(ichip).at(ichn).at(isca)=f13->GetParameter(2);

	  pedestal_map_badbcid_not_0[isca]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f13->GetParameter(1));
	  pedestal_width_map_badbcid_not_0[isca]->Fill(map_pointX[ichip][ichn] , map_pointY[ichip][ichn] , f13->GetParameter(2));
	      
	  pedestal_mean_badbcid_not_0[ichip][isca][ichn] = f13->GetParameter(1);
	  pedestal_width_badbcid_not_0[ichip][isca][ichn] = f13->GetParameter(2);

          ped_mean_integral_badbcid_not_0[ichip][isca][ichn] = ped_sca_badbcid_not_0.at(ichip).at(ichn).at(isca)->Integral(
	  		       int(pedestal_mean_badbcid_not_0[ichip][isca][ichn]-5.*pedestal_width_badbcid_not_0[ichip][isca][ichn]),
	  		       int(pedestal_mean_badbcid_not_0[ichip][isca][ichn]+5.*pedestal_width_badbcid_not_0[ichip][isca][ichn]));
		  
	}

      }

    }
    cout << "CHIP " << ichip << " END" << endl;
  }

  // Make Hit count histgram
  TH1F* hit_count_h = new TH1F("hit_count", "hit_count", 4999, 1, 5000);
  for(Int_t chip=0; chip<16; chip++){
    for(Int_t channel=0; channel<64; channel++){
      hit_count_h->Fill(hit_number[chip][channel]);
    }
  }


  fout->cd();
  Pedestal_Tree->Fill();
  Pedestal_Tree->Write();
 
  // good pedestal events (not tagged events)
  TCanvas *canvas_pedestal_map = new TCanvas("pedestal_map", "pedestal_map", 1200, 1200);
  Int_t isca = 0;
  canvas_pedestal_map->Divide(3, 3);
  canvas_pedestal_map->cd(1);
  pedestal_map[isca]->SetStats(kFALSE);
  pedestal_map[isca]->SetTitle("pedestal_map_all");
  pedestal_map[isca]->GetXaxis()->SetTitle("x");
  pedestal_map[isca]->GetYaxis()->SetTitle("y");
  pedestal_map[isca]->GetZaxis()->SetRangeUser(200, 400);
  pedestal_map[isca]->Draw("colz");
  pedestal_map[isca]->Write();

  canvas_pedestal_map->cd(2);
  pedestal_map_badbcid_0[isca]->SetStats(kFALSE);
  pedestal_map_badbcid_0[isca]->SetTitle("pedestal_map_badbcid_0");
  pedestal_map_badbcid_0[isca]->GetXaxis()->SetTitle("x");
  pedestal_map_badbcid_0[isca]->GetYaxis()->SetTitle("y");
  pedestal_map_badbcid_0[isca]->GetZaxis()->SetRangeUser(200, 400);
  pedestal_map_badbcid_0[isca]->Draw("colz");
  pedestal_map_badbcid_0[isca]->Write();
  
  canvas_pedestal_map->cd(3);
  pedestal_map_badbcid_not_0[isca]->SetStats(kFALSE);
  pedestal_map_badbcid_not_0[isca]->SetTitle("pedestal_map_badbcid_not_0");
  pedestal_map_badbcid_not_0[isca]->GetXaxis()->SetTitle("x");
  pedestal_map_badbcid_not_0[isca]->GetYaxis()->SetTitle("y");
  pedestal_map_badbcid_not_0[isca]->GetZaxis()->SetRangeUser(200, 400);
  pedestal_map_badbcid_not_0[isca]->Draw("colz");
  pedestal_map_badbcid_not_0[isca]->Write();
  
  canvas_pedestal_map->cd(4);
  pedestal_width_map[isca]->SetStats(kFALSE);
  pedestal_width_map[isca]->SetTitle("pedestal_width_map_all");
  pedestal_map[isca]->GetXaxis()->SetTitle("x");
  pedestal_width_map[isca]->GetXaxis()->SetTitle("x");
  pedestal_width_map[isca]->GetYaxis()->SetTitle("y");
  pedestal_width_map[isca]->GetZaxis()->SetRangeUser(0,7);
  pedestal_width_map[isca]->Draw("colz");
  pedestal_width_map[isca]->Write();

  canvas_pedestal_map->cd(5);
  pedestal_width_map_badbcid_0[isca]->SetStats(kFALSE);
  pedestal_width_map_badbcid_0[isca]->SetTitle("pedestal_width_map_badbcid_0");
  pedestal_map_badbcid_0[isca]->GetXaxis()->SetTitle("x");
  pedestal_width_map_badbcid_0[isca]->GetXaxis()->SetTitle("x");
  pedestal_width_map_badbcid_0[isca]->GetYaxis()->SetTitle("y");
  pedestal_width_map_badbcid_0[isca]->GetZaxis()->SetRangeUser(0,7);
  pedestal_width_map_badbcid_0[isca]->Draw("colz");
  pedestal_width_map_badbcid_0[isca]->Write();

  canvas_pedestal_map->cd(6);
  pedestal_width_map_badbcid_not_0[isca]->SetStats(kFALSE);
  pedestal_width_map_badbcid_not_0[isca]->SetTitle("pedestal_width_map_badbcid_not_0");
  pedestal_map_badbcid_not_0[isca]->GetXaxis()->SetTitle("x");
  pedestal_width_map_badbcid_not_0[isca]->GetXaxis()->SetTitle("x");
  pedestal_width_map_badbcid_not_0[isca]->GetYaxis()->SetTitle("y");
  pedestal_width_map_badbcid_not_0[isca]->GetZaxis()->SetRangeUser(0,7);
  pedestal_width_map_badbcid_not_0[isca]->Draw("colz");
  pedestal_width_map_badbcid_not_0[isca]->Write();

  canvas_pedestal_map->cd(7);
  // Make Hit Map
  TH2F* mappp = new TH2F("hit_map", "Hit Map;X[mm];Y[mm]", 32, -90, 90, 32, -90, 90);
  TH2F* chip_map = new TH2F("chip_map", "Chip Map;X[mm];Y[mm]", 32, -90, 90, 32, -90, 90);
  TH2F* asic_map = new TH2F("asic_map", "Asic Map;X[mm];Y[mm]", 4, -90, 90, 4, -90, 90);
  for(Int_t chip=0; chip<16; chip++){
    asic_map->Fill(map_pointX[chip][0], map_pointY[chip][0], chip);
    for(Int_t channel=0; channel<64; channel++){
      mappp->Fill(map_pointX[chip][channel], map_pointY[chip][channel], hit_number[chip][channel]);
      chip_map->Fill(map_pointX[chip][channel], map_pointY[chip][channel], channel);
    }
  }
  gStyle->SetOptStat(0);
  mappp->Draw("colz");
  chip_map->Draw("text same");
  asic_map->Draw("text same");
   TLine *border[2][3];
   for(Int_t i=0; i<2; i++){ 
     for(Int_t j=0; j<3; j++){ 
       if (i==0)
 	 border[i][j] = new TLine(-45+45*j,-90,-45+45*j,90);
       if (i==1)
	 border[i][j] = new TLine(-90,-45+45*j,90,-45+45*j);
       border[i][j]->SetLineColor(kGray);
       border[i][j]->SetLineWidth(1);
       border[i][j]->SetLineStyle(2);
       border[i][j]->Draw("same");
     }
   }
   if(1){
     TLatex latex;
     latex.SetTextAlign(22);  //centered
     latex.SetTextColor(17); 
     for(Int_t i=0; i<16; i++){
       if(i<8)
	 latex.DrawLatex(67.5-45*((Int_t)i/2),-22.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
       else
	 latex.DrawLatex(67.5-45*((Int_t)(i-8)/2),67.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
     }
   }
  mappp->Write();

  canvas_pedestal_map->cd(8);
  TH2F* pedestal_mean_difference = new TH2F("pedestal_mean_difference", "Pedestal_mean_difference;X[mm];Y[mm]", 32, -90, 90, 32, -90, 90);
  for(Int_t chip=0; chip<16; chip++){
    for(Int_t channel=0; channel<64; channel++){
      pedestal_mean_difference->Fill(map_pointX[chip][channel], map_pointY[chip][channel], ped_mean_badbcid_not_0[chip][channel][isca]-ped_mean_badbcid_0[chip][channel][isca]);
    }
  }
  pedestal_mean_difference->GetZaxis()->SetRangeUser(-20, 20);
  pedestal_mean_difference->Draw("colz");
  pedestal_mean_difference->Write();

  canvas_pedestal_map->cd(9);
  TH2F* pedestal_mean_ratio = new TH2F("pedestal_mean_ratio", "Pedestal_mean_ratio;X[mm];Y[mm]", 32, -90, 90, 32, -90, 90);
  for(Int_t chip=0; chip<16; chip++){
    for(Int_t channel=0; channel<64; channel++){
      pedestal_mean_ratio->Fill(map_pointX[chip][channel], map_pointY[chip][channel], 
		               ped_mean_integral_badbcid_not_0[chip][isca][channel]/ped_mean_integral_badbcid_0[chip][isca][channel]);
    }
  }
  pedestal_mean_ratio->GetZaxis()->SetRangeUser(0, 1.5);
  pedestal_mean_ratio->Draw("colz");
  pedestal_mean_ratio->Write();

  canvas_pedestal_map->Print("/home/goto/ILC/SiWECAL_2019/Analysis/Pedestal_Stability/Pdf/" + filename + "_Map.pdf");

  fout->Close();

}


