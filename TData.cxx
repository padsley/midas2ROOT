#include "Riostream.h"
#include "TData.h"
#include "TTree.h"
#include "TFile.h"
#include "Bytes.h"
#include "Byteswap.h"

#define PI 3.14159265

#include <stdlib.h>

using std::cout;
using std::endl;

TData::TData()
{
  fMap = 0;
  fGain = 0;
  fOffset = 0;
  fNchan = 512;
  for (Int_t i = 0; i < fNchan; i++) {
    fAdcOffset[i] = 0;
    fAdcGain[i] = 1;
    fAdcMap[i] = i;
    fTdcMap[i] = i;
  }
  Init();
}

TData::~TData()
{
}

void TData::Init()
{

  adcN = 0;
  tdcN = 0;
  evtNum=0;

 (void) memset(adcList,0,sizeof(UShort_t) * 512);
 (void) memset(tdcList,0,sizeof(UShort_t) * 512);
 (void) memset(adcData,0,sizeof(UShort_t) * 512);
 (void) memset(tdcData,0,sizeof(UShort_t) * 512);
 (void) memset(energy,0,sizeof(Double_t) * 512);
}

Int_t TData::SetMap(const char *mapfName,Option_t *option)
{
  const int maxAddr = fNchan;
  const int maxChan = fNchan;
  Short_t addr;
  Short_t chan;
  Short_t *mapPtr;
  char buf[80];
  ifstream mapf;
  TString opt = option;

  opt.ToUpper();
  if (! opt.CompareTo("TDC"))
    mapPtr = fTdcMap;
  else
    mapPtr = fAdcMap;
  mapf.open(mapfName);
  if (! mapf.is_open())
    cout << "SetMap: Mapping file not found -"
	 << " using default correspondence." << endl;
  else {
    cout << "SetMap: Using map file \"" << mapfName << "\"." << endl;
    while (! mapf.eof()) {
      mapf.getline(buf,80);
        sscanf(buf,"chan%hd\t%hd",&chan,&addr);
//      sscanf(buf,"%hd\t%hd",&chan,&addr);
      if (chan < maxChan && addr < maxAddr)
	mapPtr[chan] = addr;
    }
    mapf.close();
    fMap = 1;
  }
  return fMap;
}

Int_t TData::SetGain(const char *gainfName)
{
  char buf[80];
  Short_t chan;
  Double_t gain;
  ifstream gainf;
  gainf.open(gainfName);

  if (! gainf.is_open())
    cout << "SetGain: Gains file not found." << endl;
  else {
    cout << "SetGain: Using gains file \"" << gainfName << "\"." << endl;
    while (! gainf.eof()) {
      gainf.getline(buf,80);  
      sscanf(buf,"chan%hd\t%lf\n",&chan,&gain);
      fAdcGain[chan] = gain;
    }
    gainf.close();
   for (Int_t i = 0; i < fNchan; i++){
     //if(fAdcGain[i]>0.9)fAdcGain[i]=0.004;
  }
   fGain = 1;
  }
  return fGain;

}

Int_t TData::SetOffset(const char *offsetfName)
{
  char buf[80];
  Short_t chan;
  Double_t offset;
  ifstream offsetf;
  offsetf.open(offsetfName);

  if (! offsetf.is_open())
    cout << "SetOffset: Offsets file not found." << endl;
  else {
    cout << "SetOffset: Using offset file \"" << offsetfName << "\"." << endl;
    while (! offsetf.eof()) {
      offsetf.getline(buf,80);
      sscanf(buf,"chan%hd\t%lf",&chan,&offset);
      fAdcOffset[chan] = offset;    }
    offsetf.close();
    fOffset = 1;
  }
  return fOffset;
}

void TData::GetAdcMap(Short_t *adcMap)
{
  for (Int_t i = 0; i < fNchan; i++)
    adcMap[i] = fAdcMap[i];
  return;
}

void TData::GetTdcMap(Short_t *tdcMap)
{
  for (Int_t i = 0; i < fNchan; i++)
    tdcMap[i] = fTdcMap[i];
  return;
}

void TData::GetGain(Double_t *gain)
{
  for (Int_t i = 0; i < fNchan; i++){
    gain[i] = fAdcGain[i];
    if(gain[i]>0.9)gain[i]=0.004;
  }
  return;
}

void TData::GetOffset(Double_t *offset)
{
  for (Int_t i = 0; i < fNchan; i++){
    offset[i] = fAdcOffset[i];
  }
  return;
}

void TData::rootTreeW(Char_t *dataf,Option_t *option, Int_t type, Int_t blocks_to_sort) {

  Int_t Solaris=type;

  const Int_t blkSize = 16384+16384*Solaris; // in bytes
  const Int_t nWdChar = blkSize / sizeof(Char_t);
  const Int_t nWdInt = blkSize / sizeof(Int_t);

  Char_t bytes[blkSize];
  Char_t *bytesPtr = bytes;


  Int_t blocks, events,block_events, evt_len;
  Int_t scaler_evt=0, epics_evt=0;
  Short_t print, group,item,address, *half,*end_event,*end_block,end_data;
  time_t startTime;
  time_t endTime;
  ifstream datafs;
  TString opt = option;
  TString rootf = dataf;
  rootf.Append(".root");
  TFile f1(rootf.Data(),"RECREATE");
  TTree *t1=new TTree("DATA","ROOT Tree");

  datafs.open(dataf,ifstream::binary);

  if (! datafs.is_open()){
    cout << "rootTreeW: Could not find file \"" << dataf << "\"." << endl;
  }
  else {

    //f1.Close();
    t1->Branch("tdcN",&tdcN,"tdcN/I");
    t1->Branch("adcN",&adcN,"adcN/I");
    //t1->Branch("tdcList",tdcList,"tdcList[tdcN]/s");
    //t1->Branch("adcList",adcList,"adcList[adcN]/s");
    t1->Branch("tdcList",tdcList,"tdcList[256]/s");
    t1->Branch("adcList",adcList,"adcList[256]/s");
    //t1->Branch("tdcData",tdcData,"tdcData[256]/s");
    t1->Branch("tdcData",tdcData,"tdcData[tdcN]/s");
    //t1->Branch("adcData",adcData,"adcData[256]/s");
    t1->Branch("adcData",adcData,"adcData[adcN]/s");
    t1->Branch("evtNum",&evtNum,"evtNum/I");
    opt.ToLower();
    if (! opt.CompareTo("cal")) {
      if (fGain && fOffset) {
	t1->Branch("energy",&energy,"energy[adcN]/D");
      }
      else {
	cout << "rootTreeW: Option \"cal\" given"
	     << " but calibration not set." << endl;
	cout << "rootTreeW: Check calibration files/filenames." << endl;
      }
    }
    else {
//       t1->Branch("tdcData",tdcData,"tdcData[tdcN]/s");
//       t1->Branch("adcData",adcData,"adcData[adcN]/s");
    }
	cout << "#################################################################" << endl;
	cout << "Reading UK MIDAS datafile: " << dataf << endl;
	startTime = time(NULL);
	blocks=0;
	events = 0;
	Int_t scaler[32];
	for (Int_t i=0;i<32;++i){
	 scaler[i]=0;
	}

	print=0;

	if(blocks_to_sort<1)blocks_to_sort=1000000;
    cout << "blocks_to_sort: \t" << blocks_to_sort << endl;
    while (!datafs.eof() && blocks < blocks_to_sort) {			//Loop until we reach the end-of-file....
        
      if (datafs.read(bytesPtr,nWdChar)) {//Read in a block of data.
        blocks++;
        block_events=0;
        if(print)printf("Block: %d\n",blocks);
        half=(Short_t *) bytesPtr;
        if(print)printf("Half: %x\n",*half);

          if(Solaris){
              //Byte swap the data.
              for(Int_t j=0;j < nWdInt*2;j++)half[j]=Rbswap_16(half[j]);
              half+=11;
              if(print)printf("Half: %x\n",*half);
              end_data=*(half)/2+12;
              if(print)printf("End_Data: %d\n",end_data);
          }
          else{
              //Byte swap the data.
              for(Int_t j=0;j < 4;j++)half[j]=Rbswap_16(half[j]);
              for(Int_t j=13;j < nWdInt*2;j++)half[j]=Rbswap_16(half[j]);
              
              half+=10;
              if(print)printf("Half: %x\n",*half);
              end_data=*(half)/2+13;
              if(print)printf("End_Data: %d\n",end_data);
          }

          end_block=half+end_data;//Set end_block to the end of event marker of the final event in the block
          half++;

          // This is the Event Loop
          
          while (half<end_block){	//While there are still words to read in the block.....
              ++half;//This will be the length of the current event.
              evt_len=(*half-2)/2;// -2 is to remove the words containing the event length and end event marker from the word count.
              end_event=half+evt_len;// Get a pointer to the end of the event.
              if (evt_len>1000)break;//evt_len>1000 happens when gone past the end of the valid data but before the end of the block.

              if(evt_len>0){
                  Init();
                  ++half;
                    while((half < end_event)){
                        group=*half & 0x00ff;
                        item=*half >> 8 & 0x003f;
                        if(group > 0 && group < 20){//ADC Data
                            address = 32 * (group - 1) + item;
                            ++half;
                            if(*half<4095 && address < 384){
                                adcData[adcN]=*half;
                                adcList[adcN]=address;
                                energy[adcN]=(Double_t) fAdcGain[address]*(adcData[adcN]-fAdcOffset[address]);
                                ++adcN;
                            }
                        }
                        else if (group > 23 && group < 29){//TDC Data
			  //address=64*(group - 24) + item -16;
			  address=64*(group - 24) + item;
                            if(address>368)address=address-16;
                            if(address>240)address=address-16;
                            if(address>128)address=address-16;
                            ++half;
                            tdcList[tdcN]=address;
                            tdcData[tdcN]=*half;
                            ++tdcN;
                        }
                        else if (group == 30){ //// Scaler Data	
                            half+=1;
                            ++scaler_evt;
                            for(Int_t i=0; i< 32;++i){
                            scaler[i]=(*half) *65536;
                            half+=2;
                            scaler[i] += *half;
                            if(scaler[i] < 0){
                                Int_t twos_c_number=scaler[i];
                                scaler[i]=~twos_c_number +1;
                            };
                            //printf("Scaler %d = %ld\n",i,scaler[i]);
                            half+=2;
                            }
                        }
                        else if (group == 31){ // EPICS Data	
                            //Not looking at the EPics Data for now so jump to the end of the EPICS event.
                            ++epics_evt;
                            half=end_event;
                        }
                            
                        ++half;
                    }

                  if(print && !(block_events%100) && !(blocks%100)){
                      printf("Event: %d ADC_Mult= %d", block_events,adcN);
                      for(Int_t ii=0;ii<adcN;++ii)printf(" Channel: %d Data: %d",adcList[ii],adcData[ii]);
                      printf("\n");
                  }    
                  events++;
                  ++block_events;
                  half=end_event;
		  evtNum=events-1;
                  t1->Fill();
              }//if(evt_len>0)
          } // while(half < end_block)
          
  if(!(blocks%1000)){
		//TFile f1(rootf.Data(),"UPDATE");
		//f1=t1->GetCurrentFile();
		//t1->Write();
		//f1.Close();
		printf("Block: %d\n",blocks);
}		

      }	// if datafs.read
    } 	// while datafs.eof

	t1->Write();
	f1.Close();
	endTime = time(NULL);
	datafs.close();
      cout << "#################################################################" << endl;

    cout << "Output ROOT file: " << rootf.Data() << endl;
	cout << "Blocks unpacked= " << blocks << endl;
	cout << "Time taken = " << endTime - startTime << " s." << endl;
      printf("Total Events = %d = %4.2e\n", events, (Float_t) events );
      cout << "#################################################################" << endl;

 }	// else datafs.is_open
  return;
}

void TData::rootTreeR(Char_t *rootfPtr)
{

  // we should change this, so that multiple calls don't memory leak
  // TFile handles it if we give a file that doesn't really exist
//   cout << "370" << endl;
  if (rootfPtr != NULL) {
    TFile *f1 = new TFile(rootfPtr);
    TTree *t1 = (TTree *)f1->Get("DATA");
    cout << "Tree pointer is \"" << t1->GetName() << "\"." << endl;
  }
  else
    cout << "rootTreeR: Specify a root file." << endl;
  return;
}

void TData::Print()
{
  cout << "Print: Information for " << fNchan << " channels:" << endl;
  cout << "Channel:" << "Gain:" << "Offset:" 
       << "ADC map:" << "TDC map" << endl;
  for (Int_t i = 0; i < fNchan; i++)
    cout << i << "\t" << fAdcGain[i] << "\t" << fAdcOffset[i] 
	 << "\t" << fAdcMap[i] << "\t" << fTdcMap[i] << endl;
  return;
}
