#include "root_include.h"
#include "data_interface.h"
#include "functions.h"
#include "filtersparams.h"

#include <iostream>

#define OFFLINE
#define FILTERS


using namespace std;

int scanrh(TH1F *h_psd, TH1F *h_bkg, int minbin, int maxbin, float stopratio, double &psd, double &ratio)
{
  int bin = maxbin;
  float ndata,nbkg;

  do
  {
    h_psd->GetXaxis()->SetRange(minbin,bin);
    h_bkg->GetXaxis()->SetRange(minbin,bin);
    //stesso binw non serve dividere 
    ndata = h_psd->Integral();
    nbkg = h_bkg->Integral();

    if((ndata>1)&&(nbkg>1))
    {
      ratio = nbkg/ndata;
    }else
    {
      ratio = 0;
    }
    psd = h_psd->GetBinCenter(bin);

    bin--;
  } while ( (bin > minbin) && (ratio < stopratio) );

  h_psd->GetXaxis()->SetRange(1, maxbin);
  h_bkg->GetXaxis()->SetRange(1, maxbin);

  return bin - 1;

}


int scanh(TH1F *h_psd, TH1F *h_bkg, int minbin, int maxbin, float stopratio, double &psd, double &ratio)
{
  int bin = minbin;
  float ndata,nbkg;

  do
  {
    h_psd->GetXaxis()->SetRange(bin,maxbin);
    h_bkg->GetXaxis()->SetRange(bin,maxbin);
    //stesso binw non serve dividere 
    ndata = h_psd->Integral();
    nbkg = h_bkg->Integral();

    ratio = nbkg/ndata;
    psd = h_psd->GetBinCenter(bin);

    bin++;
  } while ( (bin < maxbin) && (ratio > stopratio) );

  return bin - 1;

}

void scan(TF1 *f_data, TF1 *f_bkg, float xmin,float xmax, float step, float stopratio, double &psd, double &ratio)
{

  double ndata,nbkg;
  float x = xmin;

  do
  {
    ndata = f_data->Integral(x,xmax);
    nbkg = f_bkg->Integral(x,xmax);

    if ((ndata>0)&&(nbkg>0))
    {
      ratio = nbkg/ndata;
    }else
    {
      ratio = 1;
    }

    psd = x;
    x+=step;
    //cout << x << " " << ratio << " " << ndata <<  " " << nbkg <<  endl;

  } while ( (x < xmax) && (ratio > stopratio) );

  return;

}

void scanr(TF1 *f_data, TF1 *f_bkg, float xmin,float xmax, float step, float stopratio, double &psd, double &ratio)
{

  float ndata,nbkg;
  float x = xmax;

  do
  {
    ndata = f_data->Integral(xmin,x);
    nbkg =  f_bkg ->Integral(xmin,x);

    if ((ndata>0)&&(nbkg>0))
    {
      ratio = nbkg/ndata;
    }else
    {
      ratio = 1;
    }

    psd = x;
    x-=step;

  } while ( (x > xmin) && (ratio > stopratio) );

  return;

}

void fillh(TH1F *h, DataInterface *idata, int minql, int maxql, bool filters, bool offline)
{
  float qlong,qshort,baseline,diff,psd;
  bool pu,sat;

  acqPSDParam_t params;
	event_data evento;

  params = idata->GetParams();

  for (int i = 0; i < idata->GetEntries(); i++) 
  {
    evento = idata->GetEntry(i);

    qlong = evento.qlong;
    qshort = evento.qshort;
    baseline = evento.baseline;

    if ((qlong > minql)&&(qlong < maxql))
    {
      pu = false;
      sat = false;
      if (filters)
      {
        pu = pileup(evento,params,pileupTH);
        sat = saturation(evento,params);
      }
      {
        diff = 0;

        if( (offline)&&(!pu)&&(!sat) )
        {
          GetIntegralFromScope(evento, params, &qlong, &qshort, &baseline,false);
          diff = TMath::Abs((float)evento.baseline - baseline);
        }

        if ( ((diff<baselineTH)&&(!pu)&&(!sat)) || (!filters) )
        {
            psd = (qlong - qshort) / qlong;

            h->Fill(psd);
        }

      }
    }
  }

  return;
}
int main(int argc, char *argv[])
{

	DataInterface *idata,*ibkg;
	acqPSDParam_t params;

  TH1F *h_psd,*h_psd_bkg,*h_shifted, *h_diff, *h_diff_check;
  TF1 *f_g_cf, *f_g_bkg, *f_g_bkg_s, *f_n, *ftot, *ftot_n, *ftot_g;
  TCanvas *cnv;

  int dummyargc = 0;
  char **dummyargv = 0;

  float m = 0.206098;
  float q = 0.383650;

  int nbin = 300;
  float hmin = 0;
  float hmax = 0.5;

  //valori notevoli sigma gaussiana normalizzata
//  double f05pc = 2.807033768344;
  double f1pc = 2.575829303549;
  double f2pc = 2.326347874041;
  double f5pc = 1.959963984540;
  double f10pc = 1.644853626951;


  double binw;

  char *fbkg,*fdata;
  int canale = 4;
  int center, hwidth;
  int tot_events;
  char *graphfile=0;

  int neutroni, neutroni_gaus,
      neutroni_05pc, neutroni_1pc, neutroni_2pc, neutroni_5pc, neutroni_10pc;

  float eff05pc, eff1pc, eff2pc, eff5pc, eff10pc;
  double ratio05pc, ratio1pc, ratio2pc, ratio5pc, ratio10pc, ratio1pcn;

  float eff05pch, eff1pch, eff2pch, eff5pch, eff10pch;
  double ratio05pch, ratio1pch, ratio2pch, ratio5pch, ratio10pch, ratio1pcnh;

  double psd_max, psd_min, psd_mean;
  double psd05pc, psd1pc, psd2pc, psd5pc, psd10pc;

  double psd_maxh, psd_minh, psd_meanh;
  double psd05pch, psd1pch, psd2pch, psd5pch, psd10pch;

  int bin_05pc, bin_1pc, bin_2pc, bin_5pc, bin_10pc, bin_1pcn;
  
  float fom;
//  limite alla dx del quale ho solo il x% di gamma (sul totale di gamma)
  float g05pc,g1pc,g25pc,g5pc;
  float g05pcb,g1pcb,g25pcb,g5pcb;

  bool batch = true;
  bool fitdata = true;

  TApplication *app = 0;

  if (batch) 
  {
    dummyargc = 2;
    dummyargv = new char*[2];
    dummyargv[0] = new char[5];
    dummyargv[1] = new char[5];
    strcpy(dummyargv[0],"gate");
    strcpy(dummyargv[1],"-b");
  }	

  app = new TApplication("application", &dummyargc, &dummyargv[0]);

  if (argc < 5)
  {
    cout << "Filename required\n";
    cout << "./gate <datafile> <bkgfile> <center> <hwidth>\n";
    return 0;
  }


  fdata = argv[1];
  fbkg = argv[2];
  int energy = atoi(argv[3]); //keVee
  center = TMath::Nint((energy-q) / m);
  hwidth = TMath::Nint(atoi(argv[4]) / m);
  if (argc > 5) graphfile = argv[5];

  idata = new DataInterface(fdata);
  idata->SetChannel(canale);
  params = idata->GetParams();

  ibkg = new DataInterface(fbkg);
  ibkg->SetChannel(canale);

  tot_events = idata->GetEntries();

  if (!batch)
  {
    cout << "Opening file " << argv[1] << endl;
    cout << "\n\t#####" << endl;
    cout << "channel " << params.channel << ": " << tot_events << " eventi" <<endl;
    cout << "threshold " << params.threshold << endl;
    cout << "pretrigger " << params.pretrigger << endl;
    cout << "pregate " << params.pregate << endl;
    cout << "shortgate " << params.shortgate << endl;
    cout << "longgate " << params.longgate << endl;
    cout << "numsamples " << params.numsamples << endl;
    cout << "\t#####\n" << endl;
  }

  cnv = new TCanvas("cnv","",1200,1500);
  cnv->Divide(2,2);
  h_psd = new TH1F("h_psd","",nbin,hmin,hmax);
  h_psd_bkg = new TH1F("h_psd_bkg","",nbin,hmin,hmax);
  h_shifted = new TH1F("h_shifted","",nbin,hmin,hmax);
  h_diff = new TH1F("h_diff","",nbin,hmin,hmax);
  h_diff_check = new TH1F("h_diff_check","",nbin,hmin,hmax);

  binw = h_psd->GetBinWidth(1);

//void fillh(TH2F *h, datainterface *idata, int minql, int maxql, bool filters, bool offline)
  fillh(h_psd, idata, center-hwidth, center+hwidth, false, false);
  fillh(h_psd_bkg, ibkg, center-hwidth, center+hwidth, false, false);


  float b_mean, b_rms;

  b_mean = h_psd_bkg->GetMean();
  b_rms = h_psd_bkg->GetRMS();

  //fit gamma su californio
  f_g_cf = new TF1("f_g_cf", "gaus", b_mean-b_rms, b_mean+b_rms );
//  f_bkg = new TF1("f_bkg", "gaus", b_mean-2*b_rms, b_mean+2*b_rms );
  //fit gamma su file fondo
  f_g_bkg = new TF1("f_g_bkg", "gaus" );

  //faccio il fit del picco gamma
  if (batch)
  {
    //neutroni+gamma
    h_psd->Fit(f_g_cf,"R0Q");
    //fondo gamma
    h_psd_bkg->Fit(f_g_bkg,"0Q");
  }else
  {
    cout << "FIT f_g_cf\n";
    h_psd->Fit(f_g_cf,"R0");
    cout << "FIT f_g_bkg\n";
    h_psd_bkg->Fit(f_g_bkg,"0");

  }
  //  float scale = f_data->GetParameter(0)/f_bkg->GetParameter(0);
  float scale = f_g_cf->Integral(hmin,hmax) / ( h_psd_bkg->GetEntries() * binw );
  float shift = f_g_cf->GetParameter(1)-b_mean;
  int binshift = TMath::Nint(shift/binw);

  float rms = f_g_bkg->GetParameter(2);
  float mean = f_g_bkg->GetParameter(1); 

  g05pcb = mean + rms * f1pc;
  g1pcb  = mean + rms * f2pc;
  g25pcb = mean + rms * f5pc;
  g5pcb  = mean + rms * f10pc;

  if(!batch)
  {
    cout << "scale: " << scale << " shift: " << shift << endl;
    cout << "binshift: " << binshift << endl;
  }
//  h_psd_bkg->GetXaxis()->Set(nbin,hmin+shift,hmax+shift);

  //shift and scale bkg histo
  for(int i = 1; i <= nbin;i++)
  {
   int j = i + binshift;
   if( ( j >= 1 )&&( j <= nbin) )
   {
     h_shifted->SetBinContent( j, h_psd_bkg->GetBinContent( i ) * scale);
   }

  }

//  h_shifted->Scale(scale);
  h_shifted->ResetStats();

  f_g_bkg_s = new TF1("f_g_bkg_s","gaus");
  if (batch)
  {
    h_shifted->Fit(f_g_bkg_s,"0Q");
  }else
  {
    cout << "FIT f_g_bkg scaled\n";
    h_shifted->Fit(f_g_bkg_s,"0");
  }

  h_diff->Add(h_psd);
  h_diff->Add(f_g_bkg_s,-1);
  //h_diff->Add(h_shifted-1);
  h_diff->ResetStats();

  int bin;

  bin = h_psd->FindBin(b_mean+b_rms);

  bin_10pc = scanh(h_psd,h_shifted, bin, nbin, 0.1, psd10pch, ratio10pch);
  bin = bin_10pc;

  bin_5pc = scanh(h_psd,h_shifted, bin, nbin, 0.05, psd5pch, ratio5pch);
  bin = bin_5pc;

  bin_2pc = scanh(h_psd,h_shifted, bin, nbin, 0.02, psd2pch, ratio2pch);
  bin = bin_2pc;

  bin_1pc = scanh(h_psd,h_shifted, bin, nbin, 0.01, psd1pch, ratio1pch);
  bin = bin_1pc;

  bin_05pc = scanh(h_psd,h_shifted, bin, nbin, 0.005, psd05pch, ratio05pch);
  bin = bin_05pc;

  psd_maxh = psd1pch;

  //ora vogliamo stimare il punto in cui abbiamo solo gamma ~99%
  //ipotizziamo che il picco dei neutroni sia simmetrico
  //partiamo dal bin OK (fondo gamma<1%) e fittiamo fino a...boh.
  
  //no... voglio vedere in che punto perdo solo l'1% dei neutroni
//  bin_1pcn = scanrh(h_psd,h_shifted, 1, nbin, 0.99, psd_minh, ratio1pcnh);
//  ratio1pcnh = 1 - ratio1pcnh;

//  float n_mean = h_psd->GetMean();
//  float n_rms = h_psd->GetRMS();

  float n_mean = TMath::Abs(h_diff->GetMean());
  float n_rms = h_diff->GetRMS();
// per gaus: I = costant*sigma*sqrt(2*pi)
  float sqrt2pi = 2.50662827463;
  float n_const = h_diff->GetEntries()/(n_rms*sqrt2pi); 

  for(int i = 1; i<=nbin; i++)
  {
    if (h_diff->GetBinContent(i)<0)
    {
      h_diff->SetBinContent(i,0);
    }
  }

  if (!batch) cout << "m: " << n_mean << " rms: " << n_rms << endl;

  //ora possiamo cambiare i range
  h_psd->GetXaxis()->SetRange(hmin,nbin);
//  h_psd_bkg->GetXaxis()->SetRange(hmin,nbin);
  h_shifted->GetXaxis()->SetRange(hmin,nbin);

//  f_n = new TF1("f_n", "gaus", psd1pch, n_mean+3*n_rms );
//  f_n = new TF1("f_n", "gaus", n_mean-n_rms/2, hmax );
  f_n = new TF1("f_n", "gaus", n_mean-n_rms, hmax );


  f_n->SetParameter(0,n_const);
  f_n->SetParameter(1,n_mean);
  f_n->SetParameter(2,n_rms);

  f_n->SetParLimits(0,0,n_const*3);
//  float fitminmean = n_mean-3*n_rms;
//  if (fitminmean<0) fitminmean = 0;
//  f_n->SetParLimits(1,fitminmean,hmax);
  f_n->SetParLimits(1,b_mean,hmax);

  if (batch)
  {
//    h_psd->Fit(f_n,"R0Q");
    h_diff->Fit(f_n,"R0QL");
  } else
  {
//    h_psd->Fit(f_n,"R0");
    cout << "FIT f_n su h_diff\n";
    h_diff->Fit(f_n,"R0L");
  }

  ftot = new TF1("ftot","gaus(0)+gaus(3)",hmin+2*binw,hmax);
  double p[6];
  
  for(int i =0; i < 3; i++)
  {
    p[i] = f_g_cf->GetParameter(i);
    p[i+3] = f_n->GetParameter(i);
  }

  ftot->SetParameters(p);

  if (batch)
  {
    h_psd->Fit(ftot,"IMQ0");
  }else
  {
    cout << "FIT ftot\n";
    h_psd->Fit(ftot,"IM0");
  }

  mean = TMath::Abs(ftot->GetParameter(1));
  rms = TMath::Abs(ftot->GetParameter(2));
  n_mean = TMath::Abs(ftot->GetParameter(4));
  n_rms = TMath::Abs(ftot->GetParameter(5));
  fom = getFom(mean,rms,n_mean,n_rms);

  // soglia gamma da californio
  g05pc = mean + rms * f1pc;
  g1pc  = mean + rms * f2pc;
  g25pc = mean + rms * f5pc;
  g5pc  = mean + rms * f10pc;

  ftot_n = new TF1("ftot_n","gaus");
  ftot_g = new TF1("ftot_g","gaus");

  for(int i=0;i<3;i++)
  {
    ftot_g->SetParameter(i, TMath::Abs(ftot->GetParameter(i)));
    ftot_n->SetParameter(i, TMath::Abs(ftot->GetParameter(i+3)));
  }

  //fatto questo vogliamo trovare la valle.... facciamo lo scan
  //dei bin tra i 2 picchi e troviamo il minimo.
 
  float min = h_psd->GetMaximum();
  int y=0;

  for (int i = h_psd->FindBin(b_mean); i < h_psd->FindBin(n_mean); i++)
  {
    float x = h_psd->GetBinContent(i);
    if ( x < min )
    {
      min = x;
      y = i;
    }

  }

  psd_meanh = h_psd->GetBinCenter(y);
  psd_mean = ftot->GetMinimumX(b_mean,n_mean);

  //nonononoo e` sbagliato!!!
/*  double x = ftot->GetParameter(1);
  double s = ftot->GetParameter(2); 

  psd05pc = x + f05pc * s;
  psd1pc  = x + f1pc  * s;
  psd2pc  = x + f2pc  * s;
  psd5pc  = x + f5pc  * s;
  psd10pc = x + f10pc * s;
*/

  float step = 0.001;

  scan(ftot, ftot_g, hmin,  hmax, step, 0.1,   psd10pc, ratio10pc);
  scan(ftot, ftot_g, psd10pc, hmax, step, 0.05,  psd5pc,  ratio5pc);
  scan(ftot, ftot_g, psd5pc,  hmax, step, 0.02,  psd2pc,  ratio2pc);
  scan(ftot, ftot_g, psd2pc,  hmax, step, 0.01,  psd1pc,  ratio1pc);
  scan(ftot, ftot_g, psd1pc,  hmax, step, 0.005, psd05pc, ratio05pc);
  
  psd_max = psd1pc;

  //in che punto lasciamo a sx solo l'1% dei neutroni?
  psd_min = n_mean-f2pc*n_rms;
  //con l'istogramma non riesco a calcolarlo,serve sempre il fit...
  psd_minh = psd_min;

  neutroni_gaus = TMath::Nint( ftot_n->Integral(hmin,hmax)/binw );
  neutroni = h_diff->GetEntries();

  neutroni_05pc = TMath::Nint( ftot_n->Integral(psd05pc,hmax) / binw );
  neutroni_1pc = TMath::Nint( ftot_n->Integral(psd1pc,hmax) / binw );
  neutroni_2pc = TMath::Nint( ftot_n->Integral(psd2pc,hmax) / binw );
  neutroni_5pc = TMath::Nint( ftot_n->Integral(psd5pc,hmax) / binw );
  neutroni_10pc = TMath::Nint( ftot_n->Integral(psd10pc,hmax) / binw );

  //nota che queste non sono corrette ma sarebbero 0.5%, 0.99%, 0.196%, 4.8%, 9%
  //poiche` prende la percentuale rispetto al numero di neutroni e non di eventi totali
  if (!batch)
  {
    // per gaus: I = costant*sigma*sqrt(2*pi)
    cout << "Integrale gaussiana neutroni: " <<  ftot_n->Integral(hmin,hmax)/binw << endl;
    cout << "Conteggio neutroni da fit gaussiana: " << neutroni_gaus << endl;
    cout << "Conteggio neutroni sopra soglia 1% gamma: " << neutroni_1pc << endl;
    cout << "Conteggio neutroni da istogramma: " << neutroni << endl;
    cout << "Integrale somma gaussiane: " << ftot->Integral(hmin,hmax)/binw << endl;
  }

  //da rivedere
  /*
  float counts;
  counts = h_diff->Integral();
  h_diff->GetXaxis()->SetRange(bin_05pc,nbin);
  eff05pch = h_diff->Integral() / counts;
  h_diff->GetXaxis()->SetRange(bin_1pc,nbin);
  eff1pch = h_diff->Integral() / counts;
  h_diff->GetXaxis()->SetRange(bin_2pc,nbin);
  eff2pch = h_diff->Integral() / counts;
  h_diff->GetXaxis()->SetRange(bin_5pc,nbin);
  eff5pch = h_diff->Integral() / counts;
  h_diff->GetXaxis()->SetRange(bin_10pc,nbin);
  eff10pch = h_diff->Integral() / counts;
  */
  //ok rifacciamo pero` usando la gaussiana
  //il problema e` che ci sono le oscillazioni sull'istogramma che sballano tutto
 
  float counts;
  counts = ftot_n->Integral(hmin,hmax);

  eff05pc = ftot_n->Integral(psd05pc,hmax) / counts;
  eff1pc  = ftot_n->Integral(psd1pc,hmax)  / counts;
  eff2pc  = ftot_n->Integral(psd2pc,hmax)  / counts;
  eff5pc  = ftot_n->Integral(psd5pc,hmax)  / counts;
  eff10pc = ftot_n->Integral(psd10pc,hmax) / counts;
 
  eff05pch = ftot_n->Integral(psd05pch,hmax) / counts;
  eff1pch  = ftot_n->Integral(psd1pch,hmax)  / counts;
  eff2pch  = ftot_n->Integral(psd2pch,hmax)  / counts;
  eff5pch  = ftot_n->Integral(psd5pch,hmax)  / counts;
  eff10pch = ftot_n->Integral(psd10pch,hmax) / counts;

  if (!batch)
  {
    cout << "Bin width: " << binw << endl << endl;

    cout << "DA ISTOGRAMMA" << endl;
    cout << "0.5% ratio: " << ratio05pch << " bin: " << bin_05pc << " psd: " << psd05pch <<  " eff: " << eff05pch << endl;
    cout << "1% ratio: " << ratio1pch << " bin: " << bin_1pc << " psd: " << psd1pch <<  " eff: " << eff1pch << endl;
    cout << "2% ratio: " << ratio2pch << " bin: " << bin_2pc << " psd: " << psd2pch <<  " eff: " << eff2pch << endl;
    cout << "5% ratio: " << ratio5pch << " bin: " << bin_5pc << " psd: " << psd5pch <<  " eff: " << eff5pch << endl;
    cout << "10% ratio: " << ratio10pch << " bin: " << bin_10pc << " psd: " << psd10pch <<  " eff: " << eff10pch << endl;
    cout << "PSD oltre la quale perdo solo 1% neutroni : " << psd_minh << endl;
    cout << "PSD min, valle e max (1%): " << psd_minh << " " << psd_meanh << " " << psd_maxh << endl << endl;

    cout << "DA FIT" << endl;
    cout << "0.5% ratio: " << ratio05pc << " psd: " << psd05pc << " eff: " << eff05pc << endl;
    cout << "1% ratio: "   << ratio1pc  << " psd: " << psd1pc  << " eff: " << eff1pc  << endl;
    cout << "2% ratio: "   << ratio2pc  << " psd: " << psd2pc  << " eff: " << eff2pc  << endl;
    cout << "5% ratio: "   << ratio5pc  << " psd: " << psd5pc  << " eff: " << eff5pc  << endl;
    cout << "10% ratio: "  << ratio10pc << " psd: " << psd10pc << " eff: " << eff10pc << endl;
    cout << "PSD oltre la quale perdo solo 1% neutroni : " << psd_min << endl;
    cout << "PSD min, valle e max (1%): " << psd_min << " " << psd_mean << " " << psd_max << endl;
  }

  cnv->cd(2);
  //	h_psd_bkg->Draw();
  //  f_bkg->Draw("SAME");
  h_shifted->Draw();
  f_g_bkg_s->Draw("SAME");

  cnv->cd(3);
  h_diff->Draw();
  f_n->Draw("SAME");

  h_diff_check->Add(h_psd);
  h_diff_check->Add(ftot_g,-1);
  cnv->cd(4);
  h_diff_check->Draw();
  f_n->Draw("SAME");

  cnv->cd(1);
  h_psd->Draw();
  ftot_n->SetLineColor(1);
  ftot_n->Draw("SAME");
  ftot_g->SetLineColor(2);
  ftot_g->Draw("SAME");
  ftot->SetLineColor(3);
  ftot->Draw("SAME");


  cnv->SaveAs(graphfile);

  if (fitdata)
  {
    cout << energy << " " << psd_min << " " << psd_mean << " " << psd_max << " " 
      << psd05pc   << " " << psd1pc  << " " << psd2pc   << " " << psd5pc  << " " << psd10pc << " " 
      << eff05pc   << " " << eff1pc  << " " << eff2pc   << " " << eff5pc  << " " << eff10pc << " "
      << g05pc     << " " << g05pcb  << " " << g1pc     << " " << g1pcb   << " " 
      << g25pc     << " " << g25pcb  << " " << g5pc     << " " << g5pcb   << " " << fom << endl;
  }else
  {
    cout << center << " " << psd_minh << " " << psd_meanh << " " << psd_maxh << " " 
      << psd05pch << " " << psd1pch << " " << psd2pch << " " << psd5pch << " " << psd10pch << " " 
      << eff05pch << " " << eff1pch << " " << eff2pch << " " << eff5pch << " " << eff10pch << endl;
  }

  if (!batch) app->Run(kTRUE);

  return 0;
	
}


