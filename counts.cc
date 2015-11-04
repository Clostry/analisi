#include "root_include.h"
#include "data_interface.h"
#include "functions.h"

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc,char *argv[])
{

  TCanvas *canvas;
  TH2F *h_psd;

  TGraph *gr_gate,*gr_gate_min,*gr_gate_max;

  int dummyargc = 0;
  char **dummyargv = 0;
  char *filename = 0;
  char *gatefile = 0;

  int canale = 4;
  int qlongmax = 20000;

  float th_base = 2;
  int ncol = 14;

  fstream fgate;

  DataInterface *idata;
  acqPSDParam_t params;
  event_data evento; 

  TApplication *app;

  app = new TApplication("application", &dummyargc, &dummyargv[0]);

  canvas = new TCanvas("canvas","PSD",800,500);

  if (argc != 3)
  {
    cout << "Uso :\n\t ./counts <filename> <gatefile>\n";
    return 0;
  }

  filename = argv[1];
  gatefile = argv[2];

  idata = new DataInterface(filename);
  idata->SetChannel(canale);
  params = idata->GetParams();
  
  h_psd = new TH2F("h_psd",Form("%s",filename),1000,0,qlongmax,1000,0,1);

  int n = idata->GetEntries();
  float qlong,qshort,baseline;

  for(int i = 0; i < n; i++)
  {
    evento = idata->GetEntry(i);
    if ((!saturation(evento,params)))//&&(!pileup(evento,params,70)))
    {

      GetIntegralFromScope(evento, params, &qlong, &qshort, &baseline, false);
      float diff = TMath::Abs(evento.baseline - baseline);

      if ( diff < th_base )
      {
        float psd = (qlong - qshort) / (float) qlong;
        h_psd->Fill(qlong,psd);
      }
    }
  }

  vector<int> ql;
  vector<float> psd1,psd2,psd3;

  fgate.open(gatefile);

  float tmp;
  int count = 0;

  while(!fgate.eof())
  {
    fgate >> tmp;
    ql.push_back(tmp);
    fgate >> tmp;
    psd1.push_back(tmp);
    fgate >> tmp;
    psd2.push_back(tmp);
    fgate >> tmp;
    psd3.push_back(tmp);
    for(int i=0;i<ncol-4;i++)
    {
      fgate >> tmp;
    }

    count++;
  }

  gr_gate = new TGraph(count);
  gr_gate_min = new TGraph(count);
  gr_gate_max = new TGraph(count);

  for (int i=0;i<count-1;i++)
  {
    gr_gate_min->SetPoint(i,ql[i],psd1[i]);
    gr_gate->SetPoint(i,ql[i],psd2[i]);
    gr_gate_max->SetPoint(i,ql[i],psd3[i]);
  }

  gr_gate_min->SetPoint(count-1, qlongmax, psd1[count-2]);
  gr_gate->SetPoint(count-1 ,qlongmax ,psd2[count-2]);
  gr_gate_max->SetPoint(count-1 ,qlongmax ,psd3[count-2]);
  
  gr_gate_min->SetLineColor(2);
  gr_gate->SetLineColor(3);
  gr_gate_max->SetLineColor(1);

  gr_gate_min->SetLineWidth(2);
  gr_gate->SetLineWidth(2);
  gr_gate_max->SetLineWidth(2);

  canvas->cd();
  canvas->SetLogz();
  h_psd->Draw("COLZ");
  gr_gate_min->Draw("SAME");
  gr_gate->Draw("SAME");
  gr_gate_max->Draw("SAME");

  app->Run(kTRUE);

  return 0;
}
