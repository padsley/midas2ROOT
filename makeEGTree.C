#include "TROOT.h"

Int_t makeEGTree(Char_t *filein, Option_t *option,Int_t type, Int_t blocks)
{
  gROOT->LoadMacro("TData.cxx+");
  TData data;
  cout << "7" << endl;
  data.SetOffset("si_offsets3.dat");
  data.SetGain("si_gains3.dat");

//EGTreeW options: 	1: "cal" or "" (to read in gains and offsers as defined above and calibrate)
//			2. 1 or 0. 1= Solaris (TUDA) data 0=Linux DAQ data
//			3. No. of blocks to read in. 0 = all data

  data.EGTreeW(filein,option,type,blocks);

  strcat(filein,".root");
  //cout << "18" << endl;
  data.EGTreeR(filein);
  //cout << "20" << endl;
//  TH2F *h2=new TH2F("h2","h2",4096,0,4095,320,0,319);
  TH1F *h1=new TH1F("h1","h1",4096,0,4095);
  //cout << "23" << endl;
  //gROOT->ProcessLine(".q");//This is the line which causes the "free" error. Can ignore for the moment. 
  //cout << "25" << endl;
  return 1;
}

