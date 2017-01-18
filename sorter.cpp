{
  gROOT->ProcessLine(".L Midas2Root.C+");
  for(int i=10;i<=13;i++)
    {
      gROOT->ProcessLine("Midas2Root(\"\",\"\",i,0,false)");
      gROOT->ProcessLine("Midas2Root(\"\",\"\",i,0,true)");
    }
  gROOT->ProcessLine(".q");
}
