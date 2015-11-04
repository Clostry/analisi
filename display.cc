#include "data_interface.h"
#include "root_include.h"
#include "functions.h"

#include <iostream>

using namespace std;

void display(event_data evento,acqPSDParam_t params,int n, TF1 *base, TF1 *base_online,TH1F *h_scope, TLegend *leg,float baseline)
{

    for (int i = 0; i < (int)params.numsamples;i++)
        {
              h_scope->SetBinContent(i+1,evento.samples[i]);
                }

      base->SetParameter(0,baseline);
        base_online->SetParameter(0,evento.baseline);

          base->SetLineColor(2);
            base_online->SetLineColor(3);

}

int main(int argc, char *argv[])
{

  char *filename,*graphic;
  int n;
  int canale = 4;

  int dummyargc = 0;
  char **dummyargv = 0;

  TApplication *app=0;

  TCanvas *canvas;
  TH1F *h_scope;

  DataInterface *idata;
  acqPSDParam_t params;
  event_data evento; 

  if (argc < 3) 
    return 1;

  if (argc < 4)
  {
    app = new TApplication("application", &dummyargc, &dummyargv[0]);
  }

  filename = argv[1];
  n = atoi(argv[2]);
  graphic = argv[3];

  cout << filename << " " << n << " " << graphic << endl;

  idata = new DataInterface(filename);
  idata->SetChannel(canale);
  params = idata->GetParams();
  
  evento = idata->GetEntry(n);

  canvas = new TCanvas();
  float psd = (evento.qlong - evento.qshort) / (float) evento.qlong;
  
  h_scope = new TH1F("",Form("evento %i ql %i psd %f;channel",n,evento.qlong,psd),params.numsamples,0, params.numsamples);
  for (int i = 0; i <(int)params.numsamples;i++)
  {
    h_scope->SetBinContent(i+1,evento.samples[i]);
  }

  TF1 *f = new TF1("baseline","[0]",0,params.numsamples);
  f->SetParameter(0,evento.baseline);
  f->SetLineColor(2);

  gStyle->SetOptStat(0);

  canvas->cd();
  h_scope->Draw();
  f->Draw("SAME");
   
  if (app) 
  {
    app->Run(kTRUE);
  }else
  {
    canvas->SaveAs(graphic);
  }

  return 0;
}
