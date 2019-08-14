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

void mip_pdf(std::string str){
  // Cut File name
  int path_i = str.find_last_of("/")+1;
  int ext_i = str.find_last_of(".");
  std::string pathname = str.substr(0, path_i);
  std::string filename_ = str.substr(path_i, ext_i-path_i);

  TString fullpath = str.data();
  TString filename = pathname + filename_;

  TFile *file = TFile::Open(fullpath); 

  // ================================================================================================================================================ //

  Int_t MaxChip = 16;
  Int_t MaxChannel = 64;
  Int_t MaxSca = 15;

  // ================================================================================================================================================ //

  TCanvas *canvas_mip_sca[16][64];
  TDirectoryFile *chip_dir[16];
  for(int ichip=0; ichip<MaxChip; ichip++){
    cout << "CHIP " << ichip << " START" << endl;
    TString TDirectory_Name = Form("ALL_MIP_%02dCHIP", ichip);
    chip_dir[ichip] = new TDirectoryFile();
    file->GetObject(TDirectory_Name + ";1", chip_dir[ichip]);
    if(chip_dir[ichip]==nullptr) continue;
    chip_dir[ichip]->cd();
    TCanvas *c1 = new TCanvas();
    c1->Print(filename + "_" + TDirectory_Name + ".pdf" + "[", "pdf");
    for(int ichn=0; ichn<MaxChannel; ichn++){
      TString TH1F_Name = Form("mip_chip%i_chn%i", ichip, ichn);
      canvas_mip_sca[ichip][ichn] = new TCanvas();
      gDirectory->GetObject(TString::Format("mip_map_chip%i_channel%i;1",ichip, ichn), canvas_mip_sca[ichip][ichn]);
      if(canvas_mip_sca[ichip][ichn]==nullptr) continue;
      cout << ((TH1F*)canvas_mip_sca[ichip][ichn]->GetPad(0)->FindObject(TH1F_Name))->GetEntries() << endl;
      if(((TH1F*)canvas_mip_sca[ichip][ichn]->GetPad(0)->FindObject(TH1F_Name))->GetEntries()<1000.) continue;
      canvas_mip_sca[ichip][ichn]->Draw();
      canvas_mip_sca[ichip][ichn]->Print(filename + "_" + TDirectory_Name + ".pdf", "pdf");
    }
    c1->Print(filename + "_" + TDirectory_Name + ".pdf" + "]", "pdf");
    c1->Close();

    cout << "CHIP " << ichip << " END" << endl;
  }

}


