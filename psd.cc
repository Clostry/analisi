#include "root_include.h"
#include "data_interface.h"
#include "functions.h"
#include "filtersparams.h"

#include <iostream>

using namespace std;

int main(int argc,char *argv[])
{

  TCanvas *canvas;
  TH2F *h_psd;
  TH1F *h;

  int dummyargc = 0;
  char **dummyargv = 0;
  char *filename = 0;
  bool calib = false;
  float m,q;

  int canale = 4;
  int bins = 1000;

  DataInterface *idata;
  acqPSDParam_t params;
  event_data evento; 

  TApplication *app;

  app = new TApplication("application", &dummyargc, &dummyargv[0]);

  if (argc < 2)
  {
    cout << "Uso :\n\t ./psd <filename>\n";
    return 0;
  }

  if (argc == 4)
  {
    calib = true;
    m = atof(argv[2]);
    q = atof(argv[3]);
    cout << "m: " << m << " q: " << q << endl;
  }

  filename = argv[1];

  idata = new DataInterface(filename);
  idata->SetChannel(canale);
  params = idata->GetParams();
  
  h_psd = new TH2F("h_psd",Form("%s",filename),bins,0,20000,bins,0,1);
  h = new TH1F("h",Form("%s",filename),bins,0,20000);

  int n = idata->GetEntries();
  float qlong,qshort,baseline;

  for(int i = 0; i < n; i++)
  {
    evento = idata->GetEntry(i);
    if (((!saturation(evento,params)))&&(!pileup(evento,params,70)))
    {
      GetIntegralFromScope(evento, params, &qlong, &qshort, &baseline, false);
      float diff = TMath::Abs( evento.baseline - baseline );
      if ( diff < baselineTH )
      {
        float psd = (qlong - qshort) / (float) qlong;
        h_psd->Fill(qlong,psd);
        h->Fill(qlong);
      }
    }
  }

  canvas = new TCanvas();
  if (calib)
  {
   h_psd->GetXaxis()->SetTitle("keVee");
   h_psd->GetXaxis()->Set(bins,q,q+20000*m);
   h->GetXaxis()->SetTitle("keVee");
   h->GetXaxis()->Set(bins,q,q+20000*m);
  }
  
  canvas->Divide(2,1);
  canvas->cd(1);
  canvas->SetLogz();
  h_psd->Draw("COLZ");
  canvas->cd(2);
  h->Draw();

  app->Run(kTRUE);

  return 0;
}
