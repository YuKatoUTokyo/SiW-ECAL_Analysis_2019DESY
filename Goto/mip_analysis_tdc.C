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

using namespace::std;


// ================================================================================================================================================ //
// ================================================================================================================================================ //
 
Double_t langaufun(Double_t *x, Double_t *par) {
     Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
     Double_t mpshift  = -0.22278298;       // Landau maximum location
     Double_t np = 100.0;      // number of convolution steps
     Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
     Double_t xx;
     Double_t mpc;
     Double_t fland;
     Double_t sum = 0.0;
     Double_t xlow,xupp;
     Double_t step;
     Double_t i;
     mpc = par[1] - mpshift * par[0];
     xlow = x[0] - sc * par[3];
     xupp = x[0] + sc * par[3];
     step = (xupp-xlow) / np;
     for(i=1.0; i<=np/2; i++) {
        xx = xlow + (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
        xx = xupp - (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
     }
     return (par[2] * step * sum * invsq2pi / par[3]);
}


TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
  Int_t i;
  Char_t FunName[100];
  sprintf(FunName,"Fitfcn_%s",his->GetName());
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;
  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MP","Area","GSigma");
  for (i=0; i<4; i++) {
     ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }
  his->Fit(FunName,"RQOCB0");   // fit within specified range, use ParLimits, do not plot
  ffit->GetParameters(fitparams);    // obtain fit parameters
  for (i=0; i<4; i++) {
     fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
  }
  ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
  NDF[0] = ffit->GetNDF();           // obtain ndf
  return (ffit);              // return fit function
}


Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {
  Double_t p, x = 0.0, fy, fxr, fxl;
  Double_t step;
  Double_t l, lold;
  Int_t i = 0;
  Int_t MAXCALLS = 10000;
  p = params[1] - 0.1 * params[0];
  step = 0.05 * params[0];
  lold = -2.0;
  l    = -1.0;
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = langaufun(&x,params);
     if (l < lold)
        step = -step/10;
     p += step;
  }
  if (i == MAXCALLS)
    return (-1);
  maxx = x;
  fy = l/2;
  p = maxx + params[0];
  step = params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
      step = -step/10;
      p += step;
  }
  if (i == MAXCALLS)
    return (-2);
  fxr = x;
  p = maxx - 0.5 * params[0];
  step = -params[0];
  lold = -2.0;
  l    = -1e300;
  i    = 0;
  while ( (l != lold) && (i < MAXCALLS) ) {
    i++;
    lold = l;
    x = p + step;
    l = TMath::Abs(langaufun(&x,params) - fy);
    if (l > lold)
      step = -step/10;
    p += step;
  }
  if (i == MAXCALLS)
    return (-3);
  fxl = x;
  FWHM = fxr - fxl;
  return (0);
}

// ================================================================================================================================================ //
// ================================================================================================================================================ //

void mip_analysis_tdc(std::string str, std::string slab){
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

  // ================================================================================================================================================ //
 
  // Read Chip Map
  TString mapchip_path = "map_chip.dat";
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
    map_pointY[tmp_chip][tmp_channel] = -tmp_y;
  }
  
  // ================================================================================================================================================ //

  Int_t MaxEvent = fev10->GetEntries();
  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxSca = 15;

  TH2F* mip_map[15];

  for(int isca=0; isca<MaxSca; isca++) {
    mip_map[isca]= new TH2F(TString::Format("mip_map_sca%i",isca),TString::Format("mip_map_sca%i;X[mm];Y[mm]",isca),32,-90,90,32,-90,90);
  }

  // Read data from fev10
  Int_t charge_lowGain[16][15][64];
  fev10->SetBranchAddress("charge_lowGain", charge_lowGain);

  Int_t gain_hit_low[16][15][64];
  fev10->SetBranchAddress("gain_hit_low", gain_hit_low);

  Int_t badbcid[16][15];
  fev10->SetBranchAddress("badbcid", badbcid);

  // ================================================================================================================================================ //
  
  // Pedestal Tree open  
  TString pedestal_filename = "/home/goto/ILC/SiWECAL_2019/Analysis/Pedestal_Map/" + slab + "_Pedestal.root";
  TFile *file_1 = TFile::Open(pedestal_filename); 
  TTree *Pedestal_Tree = (TTree*)file_1->Get("Pedestal_Tree");

  // Reading Calibration data from Tree
  Double_t pedestal_mean[16][15][64];
  Pedestal_Tree->SetBranchAddress("pedestal_mean", pedestal_mean);
  Double_t pedestal_mean_copy[16][15][64];

  for(Int_t event=0; event<1; event++){
    Pedestal_Tree->GetEntry(event);
    for(Int_t chip=0; chip<MaxChip; chip++){
      for(Int_t memory=0; memory<MaxSca; memory++){
	for(Int_t channel=0; channel<MaxChannel; channel++){
	  pedestal_mean_copy[chip][memory][channel]
	  = TMath::Nint(pedestal_mean[chip][memory][channel]);
	}
      }
    }
  }

  TString foutpath = "/home/goto/ILC/SiWECAL_2019/Analysis/MIP_Analysis/MIP_Map/" + true_filename + "_MIPMap.root";
  //TFile *fout = new TFile("MIP_Map/" + true_filename.data() + "_MIPMap.root" , "RECREATE");
  TFile *fout = new TFile(foutpath, "RECREATE");

  // ================================================================================================================================================ //

  // Fill to Tree Data
  /*Int_t charge_lowGain_mip[16][15][64] = {-1};
  Double_t charge_lowGain_ped_calib[16][15][64] = {-1};

  // Create New TTree
  TFile *fout = new TFile(filename + "_MIPMap.root" , "RECREATE");
  TTree *MIP_Tree = new TTree("MIP_Tree", "MIP_Tree");
  MIP_Tree->Branch("gain_hit_high", gain_hit_high,
		        "gain_hit_high[16][15][64]/I");
  MIP_Tree->Branch("badbcid", badbcid,
		        "badbcid[16][15]/I");
  MIP_Tree->Branch("charge_lowGain_ped_calib", charge_lowGain_ped_calib,
             	        "charge_lowGain_ped_calib[16][15][64]/D");
  */

  // ================================================================================================================================================ //

  TDirectory *chip_dir[16];
 
  // --------------------
  // all sca
  std::vector<std::vector<TH1F*> > mip_sca;

  for(int ichip=0; ichip<MaxChip; ichip++) {

    std::vector<TH1F*> miptemp_sca;

    for(int ichn=0; ichn<MaxChannel; ichn++) {
      TH1F *miptemp_sca2 = new TH1F(TString::Format("mip_chip%i_chn%i", ichip, ichn),
			            TString::Format("mip_chip%i_chn%i", ichip, ichn), 500, 50, 550);
      miptemp_sca.push_back(miptemp_sca2);
    }
    mip_sca.push_back(miptemp_sca);
  }

  // ================================================================================================================================================ //

  // SCA analysis
  for (int ievent=0; ievent<MaxEvent; ievent++) {
    fev10->GetEntry(ievent);

    for(int ichip=0; ichip<MaxChip; ichip++) {

      for(int isca=0; isca<MaxSca; isca++) {

	for(int ichn=0; ichn<MaxChannel; ichn++) {

	  if(charge_lowGain[ichip][isca][ichn]>10 && gain_hit_low[ichip][isca][ichn]==1){
	    mip_sca.at(ichip).at(ichn)->Fill(charge_lowGain[ichip][isca][ichn]-pedestal_mean_copy[ichip][isca][ichn]);
	  }

	} 

      }//isca

    }//ichip 
   
  }  // end first loop analysis to fill mip historgrams

  // do mip (chip/channel based) analysis
  TCanvas *canvas_mip_sca[16][64];
  for(int ichip=0; ichip<MaxChip; ichip++) {
    cout << "CHIP " << ichip << " START" << endl;
    TString TDirectory_Name = Form("ALL_MIP_%02dCHIP", ichip);
    chip_dir[ichip] = fout->mkdir(TDirectory_Name);
    chip_dir[ichip]->cd();
    for(int ichn=0; ichn<MaxChannel; ichn++) {
	canvas_mip_sca[ichip][ichn] = new TCanvas(TString::Format("mip_map_chip%i_channel%i",ichip, ichn),
		                                  TString::Format("mip_map_chip%i_channel%i",ichip, ichn), 1200, 1200);
	canvas_mip_sca[ichip][ichn]->cd();
        Double_t fr[2];
        Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
        fr[0]=0.6*mip_sca.at(ichip).at(ichn)->GetMean();
        fr[1]=2.0*mip_sca.at(ichip).at(ichn)->GetMean();
        pllo[0]=0.0; 
        pllo[1]=mip_sca.at(ichip).at(ichn)->GetMean() - 40;
        pllo[2]=0.0; 
        pllo[3]=0.0;;
        plhi[0]=0.0; 
        plhi[1]=mip_sca.at(ichip).at(ichn)->GetMean() + 10;
        plhi[2]=0.0; 
        plhi[3]=20;
        sv[0]=5.0;
        sv[1]=mip_sca.at(ichip).at(ichn)->GetMean();
        sv[2]=10000000;
        sv[3]=10;
        Double_t chisqr;
        Int_t    ndf;
        Double_t SNRPeak, SNRFWHM;
        langaupro(fp,SNRPeak,SNRFWHM);

	mip_sca.at(ichip).at(ichn)->SetTitle(TString::Format("mip_chip%i_chn%i",ichip,ichn));
	mip_sca.at(ichip).at(ichn)->SetName(TString::Format("mip_chip%i_chn%i",ichip,ichn));

	if(mip_sca.at(ichip).at(ichn)->GetEntries()>10){ //max_entries/2 ) {

            TF1 *fitsnr = langaufit(mip_sca.at(ichip).at(ichn), fr, sv, pllo, plhi, fp, fpe, &chisqr, &ndf);

            gStyle->SetOptStat(1111);
            gStyle->SetOptFit(111);
	    mip_sca.at(ichip).at(ichn)->Draw();
	    fitsnr->Draw("lsame");
	    canvas_mip_sca[ichip][ichn]->Modified();
	    canvas_mip_sca[ichip][ichn]->Update();
	    canvas_mip_sca[ichip][ichn]->Write();

	}


    }
    cout << "CHIP " << ichip << " END" << endl;
  }

  fout->cd();
  fout->Close();

}


