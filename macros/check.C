#include "data_interface.h"
#include "data_interface.cc"

void GetIntegralFromScope(event_data inc_data, acqPSDParam_t psd_params, float *qlong, float *qshort, float *baseline)
{

  float integral = 0;
  int max = 0;
  int j;

  //calcola baseline 
  *baseline = 0;
  for( j=0; j < (int) (psd_params.pretrigger-psd_params.pregate); j++)
  {
    *baseline += inc_data.samples[j];
  }
  //fai la media
  *baseline /= psd_params.pretrigger-psd_params.pregate;

  max = psd_params.pretrigger-psd_params.pregate+psd_params.shortgate;
  //non andare oltre il numero di sample
  if (max > (int)psd_params.numsamples)
    max = psd_params.numsamples;
  //calcola integrale qshort
  for( j = psd_params.pretrigger - psd_params.pregate; j < max; j++) 
  {
    integral += *baseline-inc_data.samples[j];
  }
  *qshort = integral;

  max = psd_params.pretrigger-psd_params.pregate+psd_params.longgate;
  //non andare oltre il numero di sample
  if (max > (int)psd_params.numsamples) 
    max = psd_params.numsamples;
  //calcola integrale qlong
  for( j = psd_params.pretrigger - psd_params.pregate + psd_params.shortgate; j < max; j++ )
  {
    integral += *baseline-inc_data.samples[j];
  }
  *qlong = integral;

}

void check(char *filename, int canale=4)
{
	int		tot_events, j = 0, j1 = 0, j2 = 0,
			alt=0, alt1=0, alt2=0;
	float	qlong, qshort, baseline, diff,
			th_base=4, th_qlong = 500, th_qshort = 90;
	
	DataInterface *idata;
	acqPSDParam_t params;
	event_data evento;

	TCanvas *cnv_base, *cnv_qlong, *cnv_qshort;
	TH1F	*h_base, *h_qlong, *h_qshort;
	TF1		*f1, *f2;
	
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
	h_base = new TH1F("h_base","Baseline check",200,0,10);
	cnv_qlong = new TCanvas("cnv_qlong", "Qlong check", 0, 0, 800, 500);
	h_qlong = new TH1F("h_qlong","Qlong check",10000,0,5000);
	cnv_qshort = new TCanvas("cnv_qshort", "Qshort check", 0, 0, 800, 500);
	h_qshort = new TH1F("h_qshort","Qshort check",10000,0,5000);
	
	f1 = new TF1("f1","[0]",0,params.numsamples);
	f2 = new TF1("f2","[0]",0,params.numsamples);	
	
	for (int i = 0; i < tot_events; i++) {
		evento = idata->GetEntry(i);
		GetIntegralFromScope(evento, params, &qlong, &qshort, &baseline);
		
		//baseline
		diff = TMath::Abs(evento.baseline - baseline);
		h_base->Fill( diff );
		
		if ( diff >= th_base ) {
			j++;
			if ( alt < 4 ) {
				f1->SetParameter(0,baseline);
				f2->SetParameter(0,evento.baseline);
				display(filename,i,canale,"Baseline",f1,f2);
				alt++;
			}
		}
		
		//qlong
		diff = TMath::Abs(evento.qlong - qlong);
		h_qlong->Fill( diff );
		if ( diff >= th_qlong ) {
			j1++;
			if ( alt1 < 4 ) {
				f1->SetParameter(0,baseline);
				f2->SetParameter(0,evento.baseline);
				display(filename,i,canale,"Qlong",f1,f2);
				alt1++;
			}
		}
		
		//qshort
		diff = TMath::Abs(evento.qshort - qshort);
		h_qshort->Fill( diff );
		if ( diff >= th_qshort ) {
			j2++;
			if ( alt2 < 4) {
				f1->SetParameter(0,baseline);
				f2->SetParameter(0,evento.baseline);
				display(filename,i,canale, "Qshort",f1,f2);
				alt2++;
			}	
		}
		
		
		
		
		
		cout<<"\rLoad "<< (int)(100.0*(i)/(tot_events)) <<"%";
		cout.flush();
	}
	cout << "\rLoad 100%" << endl;
	
	cout << "Eventi baseline diff.>=" << th_base << ": " << j << " (" << 100.0*j/tot_events << "%)" <<endl;
	cout << "Eventi qlong diff.>=" << th_qlong << ": " << j1 << " (" << 100.0*j1/tot_events << "%)" <<endl;
	cout << "Eventi qshort diff.>=" << th_qshort<< ": " << j2 << " (" << 100.0*j2/tot_events << "%)" <<endl;
		
	cnv_base->cd();
	h_base->Draw();
	cnv_qlong->cd();
	h_qlong->Draw();
	cnv_qshort->cd();
	h_qshort->Draw();
	
}

void display(char *filename,int n=1, int canale=4, char *title="", TF1 *base=0, TF1 *base_online=0)
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
  
  h_scope = new TH1F("",Form("%s evento %i;channel",title,n),params.numsamples,0, params.numsamples);
  for (int i = 0; i < params.numsamples;i++)
  {
    h_scope->SetBinContent(i+1,evento.samples[i]);
  }

  canvas->cd();
  h_scope->Draw();
  if (base)			{ base->SetLineColor(2); 		base->Draw("SAME"); }
  if (base_online)	{ base_online->SetLineColor(3); base_online->Draw("SAME"); }
}
