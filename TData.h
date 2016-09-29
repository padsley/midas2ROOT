// this is the data struct for UK MIDAS data - it includes some minor name modifications to fit with some of my (P Adsley's) previous sort codes
// arrays are set to the current maximums expected in the setup
#ifndef TDATA_H
#define TDATA_H

struct dataHeader {
  Char_t id[8];
  Int_t sequence;
  Short_t stream;
  Short_t tape;
  Short_t myEndian;
  Short_t dataEndian;
  Int_t dataLen;
};
typedef struct dataHeader Header_t;

class TData {
public:
  TData();
  ~TData();
Int_t SetMap(const Char_t * = NULL,Option_t * = NULL);
Int_t SetOffset(const Char_t * = NULL);
Int_t SetGain(const Char_t * = NULL);
void GetOffset(Double_t *);
void GetGain(Double_t *);
void GetAdcMap(Short_t *);
void GetTdcMap(Short_t *);
void rootTreeW(Char_t * = NULL,Option_t * = NULL, Int_t = NULL, Int_t = NULL);
void rootTreeR(Char_t * = NULL);
void Print();
  void Init();
  Header_t header;
  Int_t nEvWd;
  Int_t epicsN;
  Int_t scalerN;
  Int_t adcN;
  Int_t tdcN;
  Int_t evtNum;
  Int_t leda1flag;
  Int_t leda2flag;
UShort_t epicsList[65];
UShort_t scalerList[33];
UShort_t adcList[256];
UShort_t tdcList[256];
UShort_t epicsData[64];
UShort_t scalerData[32];
UShort_t adcData[256];
UShort_t tdcData[256];
  Double_t content;
  Double_t energy[512];
  Double_t e1;
  Double_t e2;
  Double_t etot1;
  Double_t etot;
  Double_t ec13;
  Double_t eb10;
  Double_t ebeam;
  Double_t b10;
  Double_t q;
  Double_t qa;
  Double_t theta;
  Double_t theta1;
  Double_t theta2;
  Double_t z;
  Double_t leda1_y;
  Double_t leda2_y;
  Short_t leda1;
  Short_t leda2;
  Short_t s2_1;
  Short_t s2_2;
  Short_t s2_1_strip;
  Short_t s2_2_strip;
  Short_t s2_1_sector;
  Short_t s2_2_sector;
  Double_t s2_1_y;
  Double_t s2_2_y;
  Short_t ledastrip1;
  Short_t ledastrip2;
  Short_t leda1_sector;
  Short_t leda2_sector;
private:
  Int_t fMap;
  Int_t fGain;
  Int_t fOffset;
  Int_t fNchan;
  Short_t fAdcMap[512];
  Short_t fTdcMap[512];
  Double_t fAdcOffset[512];
  Double_t fAdcGain[512];
};

#endif // TData.h
