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
	int		tot_events, j = 0, j1 = 0, j2 = 0, j3 = 0, j4 = 0,
			alt=0, alt1=0, alt2=0;
	float	qlong, qshort, baseline, diff,
			th_base=baselineTH, th_qlong = 130, th_qshort = 30, th_time=500, th_trigger = 5; 

	DataInterface *idata;
	acqPSDParam_t params;
	event_data evento;

  int dummyargc = 0;
  char **dummyargv = 0;

  char *filename,*outfile=0;
  char title[20];
  int canale =4;

  ofstream eventlog;

	TCanvas *cnv_base, *cnv_qlong, *cnv_qshort, *cnv_time, *cnv_trigger, *cnv_event[n_eventi];
	TH1F	*h_base, *h_qlong, *h_qshort, *h_time, *h_trigger, *h_scope[n_eventi];
	TF1		*f_on[n_eventi], *f_off[n_eventi];
  TLegend *leg[n_eventi];

  TApplication *app;

  app = new TApplication("application", &dummyargc, &dummyargv[0]);

  if (argc < 2)
  {
    cout << "Filename required\n";
    return 0;
  }

  cout << "Opening file " << argv[1] << endl;

  filename = argv[1];
//  if (argc == 3) canale = atoi(argv[2]);
  if (argc == 3) outfile = argv[2];

  eventlog.open(outfile);
	
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

	cnv_base = new TCanvas("cnv_base", "Baseline check", 0, 0, 800, 500);
	cnv_qlong = new TCanvas("cnv_qlong", "Qlong check", 0, 0, 800, 500);
	cnv_qshort = new TCanvas("cnv_qshort", "Qshort check", 0, 0, 800, 500);
	cnv_time = new TCanvas("cnv_time", "Time check", 0, 0, 800, 500);
	cnv_trigger = new TCanvas("cnv_trigger", "Trigger check", 0, 0, 800, 500);

  for (int i=0; i < n_eventi; i++)
  {
    h_scope[i]=0;
  }

	h_base = new TH1F("h_base","Baseline check",100,-15,30);
	h_qlong = new TH1F("h_qlong","Qlong check",1000,-500,500);
	h_qshort = new TH1F("h_qshort","Qshort check",100,-40,40);
	h_time = new TH1F("h_time","Time check",1000,0,10000);
	h_trigger = new TH1F("h_trigger","Trigger check",100,-0.5,99.5);
	
  long int time=0;
  int timediff;
	for (int i = 0; i < tot_events; i++) {
		evento = idata->GetEntry(i);
		GetIntegralFromScope(evento, params, &qlong, &qshort, &baseline,true);
		timediff = evento.timetag - time;
		//baseline
		diff = TMath::Abs(evento.baseline - baseline);
		h_base->Fill( evento.baseline - baseline );
		
		if ( diff >= th_base ) {
      if (timediff > 8000 ) 
      {
        float psd = (evento.qlong - evento.qshort) / (float) evento.qlong;
        if ((evento.qlong>500))
        {
    //      cout << i << endl;
        }
      }
      if (eventlog.is_open())
      {
        eventlog << i << endl;
      }
			j++;
			if ( alt < 4 ) {
        strcpy(title,"baseline");
	      f_off[alt] = new TF1(Form("%s_base_off%i",title,i),"[0]",0,params.numsamples);
	      f_on[alt] = new TF1(Form("%s_base_on_%i",title,i),"[0]",0,params.numsamples);
        h_scope[alt] = new TH1F(Form("%s_%i;channel",title,i),Form("%s evento %i;channel",title,i),params.numsamples,0, params.numsamples);
				display(evento,params,i,f_off[alt],f_on[alt],h_scope[alt],leg[alt],baseline);
				alt++;
			}
		}
		
		//qlong
		diff = TMath::Abs(evento.qlong - qlong);
		h_qlong->Fill( evento.qlong - qlong);
		if ( diff >= th_qlong ) {
			j1++;
			if ( alt1 < 4 ) {
        strcpy(title,"qlong");
	      f_off[alt1+4] = new TF1(Form("%s_base_off%i",title,i),"[0]",0,params.numsamples);
	      f_on[alt1+4] = new TF1(Form("%s_base_on_%i",title,i),"[0]",0,params.numsamples);
        h_scope[alt1+4] = new TH1F(Form("%s_%i;channel",title,i),Form("%s evento %i;channel",title,i),params.numsamples,0, params.numsamples);
				display(evento,params,i,f_off[alt1+4],f_on[alt1+4],h_scope[alt1+4],leg[alt1+4],baseline);
				alt1++;
			}
		}
		
		//qshort
		diff = TMath::Abs(evento.qshort - qshort);
  	h_qshort->Fill( evento.qshort - qshort );
		if ( diff >= th_qshort ) {
			j2++;
			if ( alt2 < 4) {
        strcpy(title,"qshort");
	      f_off[alt2+8] = new TF1(Form("%s_base_off%i",title,i),"[0]",0,params.numsamples);
	      f_on[alt2+8] = new TF1(Form("%s_base_on_%i",title,i),"[0]",0,params.numsamples);
        h_scope[alt2+8] = new TH1F(Form("%s_%i;channel",title,i),Form("%s evento %i;channel",title,i),params.numsamples,0, params.numsamples);
				display(evento,params,i,f_off[alt2+8],f_on[alt2+8],h_scope[alt2+8],leg[alt2+8],baseline);
				alt2++;
			}	
		}

    h_time->Fill( timediff );
    if ( timediff < th_time ) 
    {
      j3++;
    }
    time = evento.timetag;
    int trig = trigger(evento,params);
		if (params.pretrigger - trig > th_trigger) 
    {
      j4++;
    }

    h_trigger->Fill(trig);

//		cout<<"\rLoad "<< (int)(100.0*(i)/(tot_events)) <<"%";
//		cout.flush();
	}
	cout << "\rLoad 100%" << endl;
	
	cout << "Eventi baseline diff.>=" << th_base << ": " << j << " (" << 100.0*j/tot_events << "%)" <<endl;
	cout << "Eventi qlong diff.>=" << th_qlong << ": " << j1 << " (" << 100.0*j1/tot_events << "%)" <<endl;
	cout << "Eventi qshort diff.>=" << th_qshort<< ": " << j2 << " (" << 100.0*j2/tot_events << "%)" <<endl;
	cout << "Eventi timest. diff.<" << th_time << ": " << j3 << " (" << 100.0*j3/tot_events << "%)" <<endl;
	cout << "Eventi trigger .<" << th_trigger << ": " << j4 << " (" << 100.0*j4/tot_events << "%)" <<endl;
		
	cnv_base->cd();
  cnv_base->SetLogy();
	h_base->Draw();

	cnv_qlong->cd();
  cnv_qlong->SetLogy();
	h_qlong->Draw();

	cnv_qshort->cd();
  cnv_qshort->SetLogy();
	h_qshort->Draw();

	cnv_time->cd();
  cnv_time->SetLogy();
	h_time->Draw();

	cnv_trigger->cd();
	h_trigger->Draw();

  for (int i=0; i < n_eventi; i++)
  {

    if (h_scope[i]) 
    {
      cnv_event[i] = new TCanvas();
      cnv_event[i]->cd();
      h_scope[i]->Draw();
      f_on[i]->Draw("SAME");
      f_off[i]->Draw("SAME");
      leg[i] = new TLegend(0.8,0.8,1,1);
      leg[i]->AddEntry(f_off[i],"base_off","l");
      leg[i]->AddEntry(f_on[i],"base_on","l");
      leg[i]->Draw();
    }
  }

  if (eventlog.is_open()) eventlog.close();

  app->Run(kTRUE);

  return 0;
	
}


