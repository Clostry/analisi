#include "functions.cc"
#include "data_interface.h"
#include "data_interface.cc"
#include "root_include.h"

void psd(char *filename, int canale=4)
{
  TCanvas *canvas;
  TH2F *h_psd;

  DataInterface *idata;
  acqPSDParam_t params;
  event_data evento; 

  idata = new DataInterface(filename);
  idata->SetChannel(canale);
  params = idata->GetParams();
  
  h_psd = new TH2F("h_psd",Form("%s",filename),1000,0,20000,1000,0,1);

  int n = idata->GetEntries();

  for(int i = 0; i < n; i++)
  {
    evento = idata->GetEntry(n);
    float psd = (evento.qlong - evento.qshort) / (float) evento.qlong;
    h_psd->Fill(evento.qlong,psd);

  }

  canvas = new TCanvas();
  
  canvas->cd();
  h_psd->Draw("COLZ");
}
