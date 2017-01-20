{
  gROOT->ProcessLine(".L Midas2Root.C+");
  for(int i=27;i<=27;i++)
    {
      gROOT->ProcessLine("Midas2Root(\"/home/splitpole/nsi91/raw/\",\"/home/splitpole/nsi91/sorted/\",i,0,false)");
      gROOT->ProcessLine("Midas2Root(\"/home/splitpole/nsi91/raw/\",\"/home/splitpole/nsi91/sorted/\",i,0,true)");
    }
  gROOT->ProcessLine(".q");
}
