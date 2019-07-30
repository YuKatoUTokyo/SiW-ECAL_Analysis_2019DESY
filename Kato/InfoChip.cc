#ifndef INFOCHIP_CC
#define INFOCHIP_CC

#include <iostream>
#include <fstream>

using namespace std;


/*   GEOMETRY 
               by Y.Kato 
 ^+y ________________________ ______  
 |  |15 |13 |11 | 9 |~~~~|   /      |
 |  |___|___|___|___|~~~~|   / DIF  |
 |  |14 |12 |10 | 8 |~~~~|   /      |
_|_ |___|___|___|___|~~~~|   /______|
 |  | 7 | 5 | 3 | 1 |    |      |     
 |  |___|___|___|___|    | SMB  |    
 |  | 6 | 4 | 2 | 0 |    |      |
 |  |___|___|___|___|____|______|
 |／      FEV13
(X)---------+--------- > +x

board (0,0) at center of FEV13, beam direction is from front to back.

*/


class InfoChip {

public:
  int GetASUChipNumberFromChipID(int chipID){
    // Assuming the slow control numbers the chips with the PCB numeration
    return chipID;
  }

  double GetX(int chip, int channel){
    return posX[chip][channel];
  }

  double GetY(int chip, int channel){
    return posY[chip][channel];
  }

  void FindChipChannel(double x, double y, int& chip, int& channel){
    chip = -1;
    channel = -1;

    for (int i=0; i<NCHIPS; i++) {
      for (int j=0; j<NCHANNELS; j++) {
	//cout << i << "  " << j << "  " << posX[i][j] << "  " << posY[i][j] << endl;
	if (x==posX[i][j] && y==posY[i][j]) {
	  //cout << "--> OK" << endl;
	  chip = i;
	  channel = j;
	}
      }
    }
  }

  InfoChip () {
    posX = new double* [NCHIPS];
    posY = new double* [NCHIPS];
    for(int i = 0;i<NCHIPS;i++) {
      posX[i] = new double[NCHANNELS];
      posY[i] = new double[NCHANNELS];
    }

    // locations for quadrant center (1 chip=1 quadrant)
	// locations for pads in a Left or Right handed quadrant  (chip rotated by
	// 180o)
    ifstream fileChipPos("PositionChips_FEV13.txt");
    ifstream filePadPos("PositionPads_FEV13.txt");


    double val;

    // Ignore column headers (4 fields)
    char tmp[100];
    for (int i=0; i<4 && !fileChipPos.eof() && !filePadPos.eof(); i++) {
      fileChipPos >> tmp;
      filePadPos >> tmp;

    }
    // Read Chips positions
    double ChipPosX[NCHIPS];
    double ChipPosY[NCHIPS];
    for (int i=0; i<NCHIPS; i++) {
      if(!fileChipPos.eof()) { fileChipPos >> val; }
      if(!fileChipPos.eof()) { fileChipPos >> val; }
      if(!fileChipPos.eof()) { fileChipPos >> ChipPosX[i]; }
      if(!fileChipPos.eof()) { fileChipPos >> ChipPosY[i]; }
      if(fileChipPos.eof()) { cout << "Error in Chip's position file" << endl; }
    }
    // Read Pads positions  LEFT
    double PadPosX[NCHANNELS];
    double PadPosY[NCHANNELS];
	for (int j=0; j<NCHANNELS; j++) {
      if(!filePadPos.eof()) { filePadPos >> val; }
      if(!filePadPos.eof()) { filePadPos >> val; }
      if(!filePadPos.eof()) { filePadPos >> PadPosX[ (int) (val)];}
      if(!filePadPos.eof()) { filePadPos >> PadPosY[ (int) (val)];}
      if(filePadPos.eof()) { cout << "Error in Pad's position file" << endl; }
    }

    // Add them up
    for (int i=0; i<NCHIPS; i++) {
      for (int j=0; j<NCHANNELS; j++) {
	    //if (i/4==0 || i/4==2) { //for old numbering convention
	    if (i%2==0) { //calicoes June'15
		posX[i][j] = ChipPosX[i] + PadPosX[j];
		posY[i][j] = ChipPosY[i] + PadPosY[j];
		}
		else {
		posX[i][j] = ChipPosX[i] - PadPosX[j];
		posY[i][j] = ChipPosY[i] - PadPosY[j];
		}
      }
    }
    // DEBUG
    //for (int i=0; i<NCHIPS; i++) {
    //  for (int j=0; j<NCHANNELS; j++) {cout<<posX[i][j]<<" "<<posY[i][j]<<endl;}}
  }
  ~InfoChip(){
    delete posX;
    delete posY;
  }

  double **posX;
  double **posY;

private:
  enum {
    MEMDEPTH=15,
    NCHANNELS=64,
    NCHIPS=16
  };

};

#endif
