#include "root_include.h"
#include "data_interface.h"
#include "functions.h"
#include "filtersparams.h"

#include <iostream>
#include <fstream>

const int n_eventi = 12;

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

	DataInterface *idata;
	acqPSDParam_t params;
	event_data evento;

  int dummyargc = 0;
  char **dummyargv = 0;

  char *filename,*graphic=0,*txtfile=0;
  int canale =4;
  int tot_events;
  float th_rel=0.06;
  int th_abs=pileupTH;
  bool th_float = false;
  bool pu;
  ofstream outfile;

	TCanvas *cnv;
  TH2F *h_psd, *h_psd_nopu;

  TApplication *app=0;


  if (argc < 2)
  {
    cout << "Filename required\n";
    cout << argv[0] << " filename [threshold] [graphic] [eventlist]\n";
    return 0;
  }

  cout << "Opening file " << argv[1] << endl;

  filename = argv[1];
  if (argc >= 3)
  {
    th_rel = atof(argv[2]);
    th_abs = atoi(argv[2]);

    if (th_rel == th_abs) 
      th_float = false;
  }
  if (argc >= 4) 
  { 
    graphic = argv[3];
  } else
  {
    app = new TApplication("application", &dummyargc, &dummyargv[0]);
  }
  if (argc >= 5) 
  {
    txtfile = argv[4];
    outfile.open(txtfile, ios::out | ios::trunc);
  }

  cout << "pileup threshold: " << th_rel << endl;
	
	idata = new DataInterface(filename);
	idata->SetChannel(canale);
	params = idata->GetParams();
	
	tot_events = idata->GetEntries();
	
	cout << "\n\t#####" << endl;
	cout << "channel " << params.channel << ": " << tot_events << " eventi" <<endl;
	cout << "threshold " << params.threshold << endl;
	cout << "pretrigger " << params.pretrigger << endl;
	cout << "pregate " << params.pregate << endl;
	cout << "shortgate " << params.shortgate << endl;
	cout << "longgate " << params.longgate << endl;
	cout << "numsamples " << params.numsamples << endl;
	cout << "\t#####\n" << endl;

	cnv = new TCanvas("cnv", "Pileup filter check", 0, 0, 1200, 500);
  cnv->Divide(2,1);

	h_psd = new TH2F("h_psd","eventi buoni",1000,0,20000,300,0,1);
	h_psd_nopu = new TH2F("h_psd_nopu","eventi filtrati",1000,0,20000,300,0,1);
	
  float qlong,qshort,baseline,psd;
	for (int i = 0; i < tot_events; i++) 
  {
		evento = idata->GetEntry(i);
		GetIntegralFromScope(evento, params, &qlong, &qshort, &baseline,false);
    psd = (qlong - qshort) / (float) qlong;
    if ( psd < 0.5)
    {
      if (th_float)
        pu = pileup(evento,params,th_rel);
      else
        pu = pileup(evento,params,th_abs);

      if (pu)
      {
        h_psd_nopu->Fill(qlong,psd);

        if ((qlong > 500)&&(outfile.is_open()))
          outfile << i << endl;
      }else
      {
        h_psd->Fill(qlong,psd);
      }
    }
	}

  cnv->cd(1);
  h_psd->Draw("COLZ");
  cnv->cd(2);
  h_psd_nopu->Draw("COLZ");

  if (graphic)
    cnv->SaveAs(graphic);
  else 
    app->Run(kTRUE);

  return 0;
	
}


