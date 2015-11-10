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
  TH2F *h_psd_badtri;

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
  
  h_psd        = new TH2F("h_psd",Form("%s",filename),bins,0,20000,bins,0,1);
  h_psd_badtri = new TH2F("h_psd_badtri",Form("%s",filename),bins,0,20000,bins,0,1);

  float qlong,qshort,psd;

  int n = idata->GetEntries();

  for(int i = 0; i < n; i++)
  {
    evento = idata->GetEntry(i);
    qlong = evento.qlong;
    qshort = evento.qshort;
    psd = ( qlong - qshort ) / qlong;
    h_psd->Fill(qlong,psd);

    int diff = (int)params.pretrigger - trigger(evento,params);
    if ( diff > 5 )
      {
        h_psd_badtri->Fill(qlong,psd);
        if (qlong>4000)
        {
          cout << i << endl;
        }

      }
  }

  canvas = new TCanvas();
  if (calib)
  {
   h_psd->GetXaxis()->SetTitle("keVee");
   h_psd->GetXaxis()->Set(bins,q,q+20000*m);
   h_psd_badtri->GetXaxis()->SetTitle("keVee");
   h_psd_badtri->GetXaxis()->Set(bins,q,q+20000*m);
  }
  
  canvas->Divide(2,1);
  canvas->cd(1);
  canvas->SetLogz();
  h_psd->Draw("COLZ");
  canvas->cd(2);
  h_psd_badtri->Draw("COLZ");

  app->Run(kTRUE);

  return 0;
}
