#include "InfoChip.cc"

void hitmap_kato(std::string str, Bool_t iSave=false){
  
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

  // TTree_open  
  TFile *file = TFile::Open(fullpath); 
  TTree *fev10 = (TTree*)file->Get("fev10");

  // Define
  Int_t corrected_bcid[16][15];
  fev10->SetBranchAddress("corrected_bcid", corrected_bcid);

  Int_t charge_hiGain[16][15][64];
  fev10->SetBranchAddress("charge_hiGain", charge_hiGain); 
  
  Int_t charge_lowGain[16][15][64];
  fev10->SetBranchAddress("charge_lowGain", charge_lowGain);   
  
  Int_t gain_hit_high[16][15][64];
  fev10->SetBranchAddress("gain_hit_high", gain_hit_high);
  
  Int_t gain_hit_low[16][15][64];
  fev10->SetBranchAddress("gain_hit_low", gain_hit_low);
  
  Int_t badbcid[16][15];
  fev10->SetBranchAddress("badbcid", badbcid);

  Int_t event_number = fev10->GetEntries();
  cout << "event_number" << event_number << endl;

  /*TString mapchip_path = "map_chip.dat";
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
  */

  InfoChip *info = new InfoChip();
  
  // Calculating hit number
  Double_t hit_number[16][64]={0};
  for(Int_t event = 0; event<event_number; event++){
    fev10->GetEntry(event);
    for(Int_t chip=0; chip<16; chip++){
      for(Int_t sca=0; sca<15; sca++){
	if(badbcid[chip][sca]==1)continue;
	for(Int_t channel=0; channel<64; channel++){
	  if(gain_hit_high[chip][sca][channel]<0)continue;
	  hit_number[chip][channel] += gain_hit_high[chip][sca][channel];
	}
      }
    }
  }

  // Make Hit count histgram
  TH1F* hit_count_h = new TH1F("hit_count", "hit_count", 4999, 1, 5000);
  for(Int_t chip=0; chip<16; chip++){
    for(Int_t channel=0; channel<64; channel++){
      //cout << "hit_number[" << chip << "][" << channel << "] = ," << hit_number[chip][channel] << endl;
      // cout<<hit_number[chip][channel]<<endl;
      hit_count_h->Fill(hit_number[chip][channel]);
    }
  }

  // Make Hit Map
  TH2F* mappp = new TH2F("hit_map", Form("Hit Map (%s);X[mm];Y[mm]",true_filename.c_str()), 32, -90, 90, 32, -90, 90);
  TH2F* chip_map = new TH2F("chip_map", "Chip Map;X[mm];Y[mm]", 32, -90, 90, 32, -90, 90);
  TH2F* asic_map = new TH2F("asic_map", "Asic Map;X[mm];Y[mm]", 4, -90, 90, 4, -90, 90);
  for(Int_t chip=0; chip<16; chip++){
    //asic_map->Fill(map_pointX[chip][0], map_pointY[chip][0], chip);
    asic_map->Fill(info->GetX(chip,0), info->GetY(chip,0), chip);
    for(Int_t channel=0; channel<64; channel++){
      
      // Cut Condition	   
      //if(hit_number[chip][channel]>=1500){
      //  continue;
      //  }

      mappp->Fill(info->GetX(chip,channel), info->GetY(chip,channel), hit_number[chip][channel]);
      chip_map->Fill(info->GetX(chip,channel), info->GetY(chip,channel), channel);
    }
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptLogz(0);
  mappp->Draw("colz");
  chip_map->Draw("text same");
  //asic_map->Draw("text same");
  
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
	//latex.DrawLatex(67.5-45*((Int_t)i/2),67.5-45*(i%2),Form("#bf{#scale[2]{%d}}",i));
	latex.DrawLatex(67.5-45*((Int_t)i/2),-22.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
      else
	//latex.DrawLatex(67.5-45*((Int_t)(i-8)/2),-22.5-45*(i%2),Form("#bf{#scale[2]{%d}}",i));
	latex.DrawLatex(67.5-45*((Int_t)(i-8)/2),67.5-45*((i+1)%2),Form("#bf{#scale[2]{%d}}",i));
    }
  }

  if(iSave){
    TFile *fout = new TFile(Form("%s/%s_hitmap.root",pathname.c_str(),true_filename.c_str()),"RECREATE");
    mappp->Write();
    chip_map->Write();
    asic_map->Write();
    hit_count_h->Write();
    fout->Close();
  }
}
