#include "data_interface.h"
#include "data_interface.cc"
#include "root_include.h"

void display(char *filename,int n=1, int canale=4)
{
  TCanvas *canvas;
  TH1F *h_scope;

  DataInterface *idata;
  acqPSDParam_t params;
  event_data evento; 

  idata = new DataInterface(filename);
  idata->SetChannel(canale);
  params = idata->GetParams();
  
  evento = idata->GetEntry(n);

  canvas = new TCanvas();
  
  h_scope = new TH1F("",Form("evento %i;channel",n),params.numsamples,0, params.numsamples);
  for (int i = 0; i < params.numsamples;i++)
  {
    h_scope->SetBinContent(i+1,evento.samples[i]);
  }

  canvas->cd();
  h_scope->Draw();
}
