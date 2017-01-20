// C++ headers
#include <fstream>
using namespace std;

// ROOT headers
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TROOT.h>
#include <TError.h>
#include <Riostream.h>
#include <Byteswap.h>
#include <TSystem.h>

// NPLib headers
//#include "/home/padsley/codes/nptool/NPLib/include/TW1Data.h"
//#include "/home/padsley/codes/nptool/NPLib/include/TSplitPoleData.h"


void Midas2Root(TString dirin, TString dirout, Int_t run, int section, Bool_t kOld = false)
{
  //cout << "start main loop" << endl;
   gErrorIgnoreLevel = kFatal;
   //gSystem->Load("/home/padsley/codes/nptool/NPLib/lib/libNPW1.so");
   //gSystem->Load("/home/padsley/codes/nptool/NPLib/lib/libNPSplitPole.so");
   //gROOT->ProcessLine(".L /home/padsley/codes/nptool/NPLib/Detectors/W1/TW1Data.cxx+");
   
   // build input and output file names
   TString f_in  = Form("%sR%d_%d",      dirin.Data(),  run, section);
   TString f_out;
   if(kOld)
     {
       f_out = Form("%sR%d_%d_old.root", dirout.Data(), run, section);
     }
   else
     {
       f_out = Form("%sR%d_%d.root", dirout.Data(), run, section);
     }
    cout << "formatted names" << endl;
    cout << f_in.Data() << endl;
   f_out;
   
   // declare variables
   // general
   Int_t evtNum, runNum = run;
   // dsssd
   Int_t adcN, tdcN, scalarN;
   UShort_t adcList[512], tdcList[512], adcData[512], tdcData[512];
   // scaler
   Int_t scalar[32];
   // split-pole
   UShort_t SPpos, SPde, SPwire, SPplasp, SPplasg;

   short SPposFromTime;

   short SPposStart, SPposStop;

   // open data file
   ifstream fin(f_in,ifstream::binary);

   // open output ROOT file
   TFile *fout = new TFile(f_out,"recreate");
   // create output TTree
   TTree *tout = new TTree("EGTree","EG Tree");

   // declare branches
   // general
   tout->Branch("RunNumber", &runNum, "RunNumber/I");
   // dsssd
   //TW1Data *fW1Data = new TW1Data();
   //tout->Branch("W1", "TW1Data", &fW1Data);
   // sp
   //TSplitPoleData *fSplitPoleData = new TSplitPoleData();
   //tout->Branch("SplitPole", "TSplitPoleData", &fSplitPoleData);

   if (kOld) {
      tout->Branch("evtNum",&evtNum,"evtNum/I");
      tout->Branch("scalarN",&scalarN,"scalarN/I");

      tout->Branch("adcN",&adcN,"adcN/I");
      tout->Branch("tdcN",&tdcN,"tdcN/I");
      tout->Branch("adcList",adcList,"adcList[adcN]/s");
      tout->Branch("adcData",adcData,"adcData[adcN]/s");
      tout->Branch("tdcList",tdcList,"tdcList[tdcN]/s");
      tout->Branch("tdcData",tdcData,"tdcData[tdcN]/s");

      tout->Branch("SPpos",&SPpos,"SPpos/s");
      tout->Branch("SPde",&SPde,"SPde/s");
      tout->Branch("SPwire",&SPwire,"SPwire/s");
      tout->Branch("SPplasp",&SPplasp,"SPplasp/s");
      tout->Branch("SPplasg",&SPplasg,"SPplasg/s");

      tout->Branch("SPposStart",&SPposStart,"SPposStart/s");
      tout->Branch("SPposStop",&SPposStop,"SPposStop/s");

      tout->Branch("SPposFromTime",&SPposFromTime,"SPposFromTime/s");
   }

   const Int_t blkSize = 16384;
   const Int_t nWdChar = blkSize / sizeof(Char_t);
   const Int_t nWdInt = blkSize / sizeof(Int_t);
   Char_t bytes[blkSize];
   Char_t *bytesPtr = bytes;

   Int_t   blocks, events, evt_len;
   Short_t group, item, address, *end_event;
   Short_t *half, *end_block, end_data;

   blocks = 0;
   events = 0;

   // read data file
   while (!fin.eof()) {
      if (fin.read(bytesPtr, nWdChar)) {
         blocks++;
         half = (Short_t *)bytesPtr;

         for (Int_t j = 0;  j < 4; j++)        half[j] = Rbswap_16(half[j]);
         for (Int_t j = 13; j < nWdInt*2; j++) half[j] = Rbswap_16(half[j]);

         half += 10;
         end_data  = *(half) / 2 + 13;
         end_block = half + end_data;
         half++;

         while (half < end_block){
            ++half;
            evt_len   = (*half-2) / 2;
            end_event = half + evt_len;
            if (evt_len > 1000) break;

            if (evt_len > 0) {
               // reset values
	      //          fW1Data->Clear();
	      //fSplitPoleData->Clear();

               scalarN = 0;
               evtNum=0;
               adcN = 0;
               tdcN = 0;

               SPpos = 0;
               SPde = 0;
               SPwire = 0;
               SPplasp = 0;
               SPplasg = 0;

	       SPposFromTime = 0;
	       SPposStop = 0;
	       SPposStart = 0;

               for (Int_t i = 0; i < 32; ++i) {
                  scalar[i] = 0;
               }

               (void) memset(adcList,0,sizeof(UShort_t) * 512);
               (void) memset(tdcList,0,sizeof(UShort_t) * 512);
               (void) memset(adcData,0,sizeof(UShort_t) * 512);
               (void) memset(tdcData,0,sizeof(UShort_t) * 512);

               ++half;

               // decode event
               while ((half < end_event)) {
                  // adc group and channel (item)
                  group = *half & 0x00ff;
                  item  = *half >> 8 & 0x003f;
		  //cout << group << endl;
                  // adc case
                  if (group > 0 && group < 20) {
                     address = 32 * (group - 1) + item;
                     ++half;

                     if (*half < 4095 && address < 384) {
                        // dsssd
                        adcList[adcN] = address;
                        adcData[adcN] = *half;

                        // sp
                        if (adcList[adcN] == 192) SPpos   = *half;
                        if (adcList[adcN] == 193) SPde    = *half;
                        if (adcList[adcN] == 194) SPwire  = *half;
                        if (adcList[adcN] == 195) SPplasp = *half;
                        if (adcList[adcN] == 196) SPplasg = *half;

                        if (group < 7) {  // dsssd case
                           if (item < 16) {  // p side (front)
                             // fW1Data->SetFrontE(group, item, *half);
                           }
                           else {   // n side (back)
			     //   fW1Data->SetBackE(group, item%16, *half);
                           }
//			   cout << group << "\t" << item << "\t" << *half << endl;
                        }
                        else {   // sp case
			  //           if (address == 192) fSplitPoleData->SetPosition(*half);
                          // if (address == 193) fSplitPoleData->SetDeltaE(*half);
                          // if (address == 194) fSplitPoleData->SetWire(*half);
                          // if (address == 195) fSplitPoleData->SetPlasticP(*half);
                          // if (address == 196) fSplitPoleData->SetPlasticG(*half);
                        }

                     }
                     ++adcN;
                  }

                  // tdc case
                  else if (group > 19 && group < 29){
                     address = 64*(group - 24) + item;
                     ++half;
                     if (address < 112) {    // dsssd's: remove TAC and trigger information
                        Int_t det   = address/16 + 1;
                        Int_t strip = address%16;
                        if (address < 96) {  // p-side (front)
                          // fW1Data->SetFrontT(det, strip, *half);
                        }
                        else {   // n-side (back)
                          // fW1Data->SetBackT(det-6, strip, *half);   // only D1 has timing for back signals
                        }
                     }

                     // SP
		     // if (address == 117) fSplitPoleData->SetTime1(*half);
		     // if (address == 126) fSplitPoleData->SetTime2(*half);

                     tdcList[tdcN] = address;
                     tdcData[tdcN] = *half;

                     if (address == 126 && *half < 800) {
//                        cout << tdcN << "\t" << group << "\t" << item << "\t" << address << "\t" << *half << "\t" << fSplitPoleData->GetTime2() << endl;
                     }

		     if(address==113 && SPposStop!=0)cout << "ACHTUNG! SPposStop" << endl;
		     if(address==114 && SPposStart!=0)cout << "ACHTUNG! SPposStart" << endl;


		     if(address==113)SPposStop = *half;
		     if(address==114)SPposStart = *half;
                     ++tdcN;
                  }
		  
		  

                  // scaler case
                  else if (group == 30){
                     for (Int_t i = 0; i < 3; i++) {
                        half++;
                        if (*half < 0) scalar[i] = 65536 + (*half);
                        else scalar[i] = *half;

                        half += 2;
                        scalar[i] += (*half) * 65536;
                        if (i == 2) {
                           scalarN = scalar[i];
			   // fSplitPoleData->SetTick(scalar[i]);
                        }
                        half++;
                     }
                     half--;
                  }

                  // ??
                  else if (group == 31){
                     half = end_event;
                  }
                  ++half;
               }
	       
	       SPposFromTime = SPposStop - SPposStart;

               events++;
               evtNum = events-1;
               tout->Fill();
               half = end_event;
            }
         }

         if(!(blocks%500)){
            cout << "\rProcessing Block:  " << blocks;
            cout.flush();
         }
      }
   }

   cout << "\r----------------------Processed All Blocks-----------------------" << endl;
   cout << "blocks " << blocks << endl;

   tout->Write();
   fout->Close();
   fin.close();

//   gROOT->ProcessLine(".q");
}
