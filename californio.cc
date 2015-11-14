#include <iostream>
#include <sstream>
#include <fstream>
#include <typeinfo>

#include <root_include.h>

#include <boost/program_options.hpp>

#include <data_interface.h>
#include <functions.h>

namespace po = boost::program_options;

using namespace std;

template<class T> float getFOM(T *h, float &prob = 0)
{

  float delta;
  float fom;
  double fwhm1,fwhm2;
  double mean1,mean2;
  double max[3];
  double tmp;
  int peaks;
  int entries;

  TSpectrum *spec;
  TF1 *f1,*f2,*ftot;

  spec = new TSpectrum(3);
/*
  threshold: threshold value in % for selected peaks, peaks with amplitude less than threshold*highest_peak/100 are ignored, see manual.o
  Int_t Search(const TH1* hist, Double_t sigma = 2, Option_t* option = "", Double_t threshold = 0.050000000000000003)
  */
  //prendi per buoni il primo e l'ultimo picco
//  peaks = spec->Search(h->Rebin(4),0.0002,"nobackground,noMarkov",0.01);
  h->SetAxisRange(0,0.5);
  //peaks = spec->Search(h->Rebin(4),0.02,"nobackground,noMarkov",0.001);
//  h->Rebin(2);
//  h->SetBinContent(1,0);
  entries = h->GetEntries();
  if (entries > 5000)
  {
    T *rebinned =(T*) h->Rebin(2,"rebinned");
    peaks = spec->Search(rebinned,1,"nobackground,noMarkov",0.05);
  } else
  {
    peaks = spec->Search(h,1,"nobackground,noMarkov",0.001);
  }
  //peaks = spec->Search(h->Rebin(4),0.01,"nobackground,noMarkov",0.0005);
  
  for (int i = 0; i< peaks;i++)
  {
    max[i] = spec->GetPositionX()[i];
  }

  if (peaks < 2)
  {
    return 0;
  }

  //ordina: funziona con 2 o tre picchi
 
  for (int i=0; i < peaks-1; i++)
  {
    if (max[i]>max[i+1]) 
    {
      tmp = max[i+1];
      max[i+1] = max[i];
      max[i] = tmp;
    }
  }

  if (max[0]>max[1])
  {
    tmp = max[1];
    max[1] = max[0];
    max[0] = tmp;
  }
  
    
  float max1 = h->GetBinCenter(h->GetMaximumBin());
  float max2 = max[peaks-1];
  delta = (max2-max1)/4;
  //picco gamma
  f1 = new TF1("f1", "gaus", max1-delta*1.5, max1+delta );
  //picco neutroni
  f2 = new TF1("f2", "gaus", max2-delta, max2+4*delta );
  //ftot = new TF1("ftot","gaus(0)+gaus(3)",par[4]-2*par[5],par[1]+3*par[2]);
  ftot = new TF1("ftot","gaus(0)+gaus(3)",max1-delta*1.5,max2+4*delta);

  //rosso gamma
  f1->SetLineColor(2);
  //verde neutroni
  f2->SetLineColor(3);
  //blu totale
  ftot->SetLineColor(4);

  f1->SetParameter(1, max1);
  f2->SetParameter(1, max2);
/*
  h->Fit(f1,"QLR" );
  double m = f1->GetParameter(1);
  double s = f1->GetParameter(2);
  f1->SetRange( m-3*s , m+3*s);
  f2->SetRange( m+2*s , max[peaks-1]+4*delta);
  h->Add(f1,-1);
  h->Fit(f2,"QLR+");
*/
  double par[6];
  h->Fit(f1,"QLR" );

  if (entries<5000)
  {
    h->Fit(f2,"QLLR+");
  }else
  {
    h->Fit(f2,"QLR+");
  }


  f1->GetParameters(&par[0]);
  f2->GetParameters(&par[3]);
  
  //ftot = new TF1("ftot","gaus(0)+gaus(3)");


  ftot->SetParameters(par);
	
//  double xmin = par[1]-3*par[2];
//  if (xmin > par[4]-3*par[5]) xmin = par[4]-3*par[5]
  ftot->SetRange(par[1]-3*par[2],par[4]+4*par[5]);
  if (entries<20000)
  {
    /*
    if (entries<5000)
    {
      h->Fit(ftot,"QLLR+");
    }else
    {
      h->Fit(ftot,"QLR+");
    }
    */
    h->Fit(ftot,"QLR+");
  } else
  {
    h->Fit(ftot,"QR+");
  }

  ftot->GetParameters(&par[0]);

  //leggi parametri del fit
  mean1 = par[1];
  fwhm1 = par[2]*2.3548;
  mean2 = par[4];
  fwhm2 = par[5]*2.3548;
  

  //calcola fom
  fom = ( mean2 - mean1 ) / ( fwhm1 + fwhm2 );

  prob = ftot->GetProb();

  delete spec;
  delete f1;
  delete f2;
  delete ftot;

  return fom;
}

int main(int argc, char *argv[])
{

  int dummyargc=0;
  char **dummyargv=0;

  bool draw_psd = true;

  int psd_nbin = 1000;
  int psd_maxqlong = 38000;
  float psd_max_value = 1;

  int pregate, longgate, shortgate;
  vector<int> threshold(2);
  int event;
  int channel;
  string root_file, prefix, extension, str_th, gate_file, str_optim;
  ofstream optim_file;
//  string config_file

  acqPSDParam_t psd_params;
  DataInterface *idata = 0;

  TApplication *application = NULL;
  TH2F *h_psd_q, *h_psd_q_online;
  TH1F *h_scope = 0;
  TH1D *h_fom, *h_fom_online,*h_qlong_online,*h_qlong;
  TGraph *gr_fom_th = 0,
         *gr_fom_th_online = 0,
         *gr_fom_pg = 0,
         *gr_fom_sg = 0,
         *gr_fom_lg = 0,
         *gr_gate = 0;
  TCanvas *cnv_psd, *cnv_psd_online, *cnv_fom_th, *cnv_fom_th_online, 
          *cnv_fom_pg, *cnv_fom_sg, *cnv_fom_lg, *cnv_scope, *cnv_fom, *cnv_fom_online, *cnv_dummy,
          *cnv_qlong, *cnv_qlong_online;
  string str_psd, str_psd_online, str_fom_th, str_fom_th_online, str_fom_pg, str_fom_sg, 
         str_fom_lg, str_scope, str_fom, str_fom_online;

  //program options from boost library
  try 
  {
    po::options_description generic("Generic options");
    generic.add_options()
      ("help,h", "produce help message")        
//      ("config,c", po::value<string>(&config_file), "load config file")
      ;

    po::options_description config("Configuration");
    config.add_options()
      ("root-file,f",  po::value<string>(&root_file),                       "root file to analyze")
      ("prefix,e",     po::value<string>(&prefix)->default_value(""),       "prefix for output graphs filenames")
      ("extension,x",  po::value<string>(&extension)->default_value("png"), "extension type for graphics")
      ("channel,c",    po::value<int>(&channel)->default_value(4),                            "read specific channel")
      ("event,n",      po::value<int>(&event),                              "display event n. arg 1-nevents")
      ("pre-gate,p",   po::value<int>(&pregate),   "set pregate")
      ("short-gate,s", po::value<int>(&shortgate), "set short gate")
      ("long-gate,l",  po::value<int>(&longgate),  "set long gate")
      //("threshold,t",  po::value<int>(&threshold)->default_value(0),  "discard qlong<threshold")
      ("threshold,t",  po::value< vector<int> >()->multitoken(),      "discard qlong<threshold (min and max), 0 means no threshold")
      ("scan-th",      po::value< vector<int> >()->multitoken(),      "test different qlong threshold values between specified range")
      ("scan-pg",      po::value< vector<int> >()->multitoken(),      "test different pre gate values between specified range (step=1bin)")
      ("scan-sg",      po::value< vector<int> >()->multitoken(),      "test different short gate values between specified range (step=1bin)")
      ("scan-lg",      po::value< vector<int> >()->multitoken(),      "test different long gate values between specified range (step=1bin)")
      ("calibration,a",po::value< vector<float> >()->multitoken(),    "requires slope and intercept")
      ("gatepoints,g", po::value<string>(&gate_file),                 "2 column gate points") 
      ("non-interactive,q", "run and print graphs without loading the GUI")
//      ("other-format,o", "use other data format")
      ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);

    po::options_description config_file_options;
    config_file_options.add(config);

    po::variables_map vm;
    store(po::parse_command_line(argc, argv, cmdline_options), vm);
    notify(vm);
/*
    if (vm.count("config"))
    {
      ifstream ifs(config_file.c_str());
      if (!ifs)
      {
        cout << "can not open config file: " << config_file << "\n";
        return 0;
      }
      else
      {
        store(parse_config_file(ifs, config_file_options), vm);
        notify(vm);
      }
    }
*/
    if (vm.count("help"))
    {
      cout << cmdline_options << "\n";
      return 1;
    }

    if (vm.count("non-interactive"))
    {
      dummyargc = 2;
      dummyargv = new char*[2];
      dummyargv[0] = new char[5];
      dummyargv[1] = new char[5];
      strcpy(dummyargv[0],argv[0]);
      strcpy(dummyargv[1],"-b");
    }

    application = new TApplication("application", &dummyargc, &dummyargv[0]);

    
    idata = new DataInterface(root_file.c_str());
    idata->SetChannel(channel);
    psd_params = idata->GetParams();

/*    if (! vm.count("other-format"))
    {
      tree   = (TTree*)infile->Get("acq_tree_0");
      branch = tree->GetBranch("acq_ch4");
      psd_params.numsamples = 400;
    }
    else
    {
      tree   = (TTree*)infile->Get("datatree");
      branch = tree->GetBranch("ACQ_ch2");
      psd_params.numsamples = 80;
    }
*/

    if (vm.count("pre-gate"))
    {
      psd_params.pregate = pregate;
    }
    if (vm.count("short-gate"))
    {
      psd_params.shortgate = shortgate;
    }
    if (vm.count("long-gate"))
    {
      psd_params.longgate = longgate;
    }

    cout << "****PSD parameters****" << endl;
    cout << "channel " << (psd_params.channel) << endl;
    cout << "threshold " << (psd_params.threshold) << endl; //trigger threshold, unused
    cout << "pretrigger " << (psd_params.pretrigger) << endl;
    cout << "pregate " << (psd_params.pregate) << endl;
    cout << "shortgate " << (psd_params.shortgate)  << endl;
    cout << "longgate " << (psd_params.longgate) << endl;
    cout << "numsamples " << (psd_params.numsamples) << endl << endl;
    /*cout << "pregate " << (psd_params.pregate = pregate) << endl;
    cout << "shortgate " << (psd_params.shortgate = shortgate)  << endl;
    cout << "longgate " << (psd_params.longgate = longgate) << endl;
*/


    // Adjust output filenames
    str_psd =           prefix + "_psd."            + extension;
    str_psd_online =    prefix + "_psd_online."     + extension;
    str_fom_th =        prefix + "_fom_th."         + extension;
    str_fom_th_online = prefix + "_fom_th_online."  + extension;
    str_fom_pg =        prefix + "_fom_pg."         + extension;
    str_fom_sg =        prefix + "_fom_sg."         + extension;
    str_fom_lg =        prefix + "_fom_lg."         + extension;
    str_fom =           prefix + "_fom."            + extension;
    str_fom_online =    prefix + "_fom_online."     + extension;
    str_scope =         prefix + "_scope."          + extension;
    str_optim =         prefix + "_optim.txt";


    // root style tricks
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat("en");
    //gStyle->SetOptStat("enmr");
    //  gStyle->SetOptFit(1111);
    gStyle->SetPalette(1);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetLabelFont(52, "xyz");
    gStyle->SetTitleFont(72, "xyz");
    gStyle->SetLabelSize(0.03,"xyz");
    gStyle->SetTextSize(0.03);
    gStyle->SetLineScalePS(1);
    //  gStyle->SetHistLineColor(kOrange+10);
    //  gStyle->SetHistFillColor(kOrange+10);
    gStyle->SetFuncColor(kBlue);
    gStyle->SetFuncWidth(2);
    
    int tot_events = idata->GetEntries();
    cout << tot_events << " Eventi totali" << endl << endl;

    //h_qlong = new TH1F("h_qlong", "qlong", 3000, 0, 32000 );
    //h_qlong_online = new TH1F("h_qlong_online", "qlong", 3000, 0, 32000 );

    /*
    int saturated = 0;
    int saturated_min = 0;
    int saturated_max = 0;
    cout << "Eventi saturati: ";
    for( int i = 0; i < tot_events; i++ )
    {
      event_data inc_data = idata->GetEntry(i);
      //h_qlong_online->Fill(inc_data.qlong);
      if (!saturation(inc_data,psd_params))
      {
      } 
      else
      {
        if (saturated == 0)
        {
          saturated_min = saturated_max = inc_data.qlong;
        } else
        {
          if (inc_data.qlong > saturated_max)
          {

            cout << i << " ";
            saturated_max = inc_data.qlong;
          }
          if (inc_data.qlong < saturated_min)
          { 
            cout << i << " ";
            saturated_min = inc_data.qlong;
          }
        }
        saturated++;
      }
    }
    cout << endl << saturated << " eventi sono saturati\n";
    cout << "qlong min " << saturated_min << endl;
    cout << "qlong max " << saturated_max << endl;
    */

    //root wants it....
    cnv_dummy = new TCanvas("dummy", "dummy", 0, 0, 100, 100);

    if (vm.count("event"))
    {
      if ((event >= 0)&&(event < tot_events ))
      {
        event_data inc_data = idata->GetEntry(event);
        h_scope = new TH1F("h_scope","Scope",psd_params.numsamples,0,psd_params.numsamples);
        for ( int i = 0; i < (int) psd_params.numsamples; i++ )
        {
          h_scope->SetBinContent(i + 1,inc_data.samples[i]);
        }
      }
      else
      {
        cout << "Invalid event number.\n";
      }

    }

    threshold[0] =  0;
    threshold[1] =  0;

    if (vm.count("threshold"))
    {
      threshold = vm["threshold"].as<vector<int> >();
      if (threshold.size() != 2)
      {
        cout << "--threshold requires 2 values\n";
        return 1;
      }
    }

    cout << "Lower threshold: " << threshold[0] << " Upper threshold: " << threshold[1] << endl;

    // compute FoM vs pregate 
    if (vm.count("scan-pg"))
    {
      draw_psd = false;
      
      vector<int> range = vm["scan-pg"].as<vector<int> >();
      if (range.size() != 2)
      {
        cout << "--scan-pg requires min and max value" << endl;
        return 1;
      }

      int npoints = range[1]-range[0]+1;
      int min = range[0];
      int max = range[1];
      gr_fom_pg = new TGraph(npoints);

       
      cout << "Scan pre-gate " << min << "-" << max << ": ";
      
      gr_fom_pg = new TGraph(npoints);

      acqPSDParam_t tmp_params = psd_params;
  
      if (!optim_file.is_open())
      {
        optim_file.open(str_optim.c_str(),ofstream::out|ofstream::trunc);
      }

      TH1::AddDirectory(kFALSE);

      for ( int i = min; i <= max; i++ )
      {

        tmp_params.pregate = i;
        cout << i << " "; 
        cout.flush();
        TH2F *h = new TH2F("","", psd_nbin, 0,psd_maxqlong ,300 ,0 ,psd_max_value  );
        fill_psd_q(h, idata, tmp_params, threshold, true  );
         
        TH1D *h_p = h->ProjectionY("");

        float fom,prob;
        fom = getFOM(h_p,prob);
        optim_file << tmp_params.longgate << " " << tmp_params.shortgate << " " << tmp_params.pregate << " " << fom << " " << prob << endl;
        optim_file.flush();

        gr_fom_pg->SetPoint(i-min,i,fom);

        delete h;
        delete h_p;
      
      }

      cout << endl;

      TH1::AddDirectory(kTRUE);

      gr_fom_pg->SetMarkerStyle(20);
      gr_fom_pg->SetMarkerColor(kOrange+10);
      gr_fom_pg->GetXaxis()->SetTitle("pre-gate");
      gr_fom_pg->GetYaxis()->SetTitle("FoM");
      gr_fom_pg->SetTitle("FoM vs pre-gate");
      
    }

    // compute FoM vs short gate
    if (vm.count("scan-sg"))
    {

      draw_psd = false;
      vector<int> range = vm["scan-sg"].as<vector<int> >();
      if (range.size() != 2)
      {
        cout << "--scan-sg requires min and max value" << endl;
        return 1;
      }

      int npoints = range[1]-range[0]+1;
      int min = range[0];
      int max = range[1];
      
      cout << "Scan short-gate " << min << "-" << max << ": ";
      
      gr_fom_sg = new TGraph(npoints);

      acqPSDParam_t tmp_params = psd_params;

      if (!optim_file.is_open())
      {
        optim_file.open(str_optim.c_str(),ofstream::out|ofstream::trunc);
      }

      TH1::AddDirectory(kFALSE);

      for ( int i = min; i <= max; i++ )
      {

        tmp_params.shortgate = i;
        cout << i << " "; 
        cout.flush();

        TH2F *h = new TH2F("","", psd_nbin, 0,psd_maxqlong ,300 ,0 ,psd_max_value  );
        fill_psd_q(h, idata, tmp_params, threshold, true  );
         
        TH1D *h_p = h->ProjectionY("");

        float fom,prob;
        fom = getFOM(h_p,prob);
        optim_file << tmp_params.longgate << " " << tmp_params.shortgate << " " << tmp_params.pregate << " " << fom << " " << prob << endl;
        optim_file.flush();

        gr_fom_sg->SetPoint(i-min,i,fom);

        delete h;
        delete h_p;
      
      }
      
      cout << endl;
      
      TH1::AddDirectory(kTRUE);

      gr_fom_sg->SetMarkerStyle(20);
      gr_fom_sg->SetMarkerColor(kOrange+10);
      gr_fom_sg->GetXaxis()->SetTitle("short-gate");
      gr_fom_sg->GetYaxis()->SetTitle("FoM");
      gr_fom_sg->SetTitle("FoM vs short-gate");
    }

    // compute FoM vs long gate
    if (vm.count("scan-lg"))
    {

      draw_psd = false;
      vector<int> range = vm["scan-lg"].as<vector<int> >();
      if (range.size() != 2)
      {
        cout << "--scan-lg requires min and max value" << endl;
        return 1;
      }

      int inc = 10;
      int npoints = (range[1]-range[0])/inc+1;
      int min = range[0];
      int max = range[1];
      gr_fom_lg = new TGraph(npoints);

 
      cout << "Scan long-gate " << min << "-" << max << ": ";
      
      gr_fom_lg = new TGraph(npoints);

      acqPSDParam_t tmp_params = psd_params;

      if (!optim_file.is_open())
      {
        optim_file.open(str_optim.c_str(),ofstream::out|ofstream::trunc);
      }

      TH1::AddDirectory(kFALSE);

      for ( int i = 0; i < npoints; i++ )
      {

        tmp_params.longgate = min + i*inc;
        cout << tmp_params.longgate << " ";
        cout.flush();
        TH2F *h = new TH2F("","", psd_nbin, 0,psd_maxqlong ,300 ,0 ,psd_max_value  );
        fill_psd_q(h, idata, tmp_params, threshold, true  );
         
        TH1D *h_p = h->ProjectionY("");

        float fom,prob;
        fom = getFOM(h_p,prob);
        optim_file << tmp_params.longgate << " " << tmp_params.shortgate << " " << tmp_params.pregate << " " << fom << " " << prob << endl;
        optim_file.flush();

        gr_fom_lg->SetPoint(i,tmp_params.longgate,fom);

        delete h;
        delete h_p;
      
      }

      cout << endl;

      TH1::AddDirectory(kTRUE);

      gr_fom_lg->SetMarkerStyle(20);
      gr_fom_lg->SetMarkerColor(kOrange+10);
      gr_fom_lg->GetXaxis()->SetTitle("long-gate");
      gr_fom_lg->GetYaxis()->SetTitle("FoM");
      gr_fom_lg->SetTitle("FoM vs long-gate");

    }

     
    if (draw_psd)
    {
      h_psd_q_online =  new TH2F("h_psd_q_online", "PSD vs qlong Online;qlong;PSD", psd_nbin, 0,psd_maxqlong ,300 ,0 ,psd_max_value  );
      h_psd_q        =  new TH2F("h_psd_q",        "PSD vs qlong Online;qlong;PSD", psd_nbin, 0,psd_maxqlong ,300 ,0 ,psd_max_value  );

      vector<int> th(2);

      th[0] = 0;
      th[1] = 0;

      fill_psd_q(h_psd_q_online, idata, psd_params, th, false );
      fill_psd_q(h_psd_q,        idata, psd_params, th, true  );

      h_qlong_online = h_psd_q_online->ProjectionX("h_qlong_online");
      h_qlong = h_psd_q->ProjectionX("h_qlong");

      if (threshold[0]!=0)
      {
        if ( threshold[1]==0 )
        {
          h_fom = h_psd_q->ProjectionY("h_fom",threshold[0]*psd_nbin/psd_maxqlong);
          h_fom_online = h_psd_q_online->ProjectionY("h_fom_online",threshold[0]*psd_nbin/psd_maxqlong);
        }
        else
        {
          h_fom = h_psd_q->ProjectionY("h_fom",threshold[0]*psd_nbin/psd_maxqlong,threshold[1]*psd_nbin/psd_maxqlong);
          h_fom_online = h_psd_q_online->ProjectionY("h_fom_online",threshold[0]*psd_nbin/psd_maxqlong,threshold[1]*psd_nbin/psd_maxqlong);
        }
      }
      else
      {
        h_fom = h_psd_q->ProjectionY("h_fom");
        h_fom_online = h_psd_q_online->ProjectionY("h_fom_online");
      }

      float prob;
      cout << "FoM from online data: " << getFOM(h_fom_online,prob) << endl;
      cout << "FoM from data: " << getFOM(h_fom,prob) << endl;
      
      


      // compute FoM vs energy threshold
      if (vm.count("scan-th")) 
      {

        vector<int> range = vm["scan-th"].as<vector<int> >();
        if (range.size() != 2) 
        {
          cout << "--scan-th requires min and max value" << endl;
          return 1;
        }

        int npoints = 20;
        int min_bin = (range[0]*psd_nbin)/psd_maxqlong;
        int max_bin = (range[1]*psd_nbin)/psd_maxqlong;
        int th_bin;
        gr_fom_th = new TGraph(npoints);
        gr_fom_th_online = new TGraph(npoints);

        TH1::AddDirectory(kFALSE);

        for ( int i = 1; i <= npoints; i++ )
        {

          float p;
          th_bin = min_bin + (max_bin-min_bin)*i/npoints;
          TH1D *h = h_psd_q_online->ProjectionY("",th_bin);
          gr_fom_th_online->SetPoint(i-1,th_bin*psd_maxqlong/psd_nbin,getFOM(h,p));

          delete h;

          h = h_psd_q->ProjectionY("",th_bin);
          gr_fom_th->SetPoint(i-1,th_bin*psd_maxqlong/psd_nbin,getFOM(h,p));

          delete h;

        } 

        TH1::AddDirectory(kTRUE);

        gr_fom_th_online->SetMarkerStyle(20);
        gr_fom_th_online->SetMarkerColor(kOrange+10);
        gr_fom_th_online->GetXaxis()->SetTitle("qlong threshold");
        gr_fom_th_online->GetYaxis()->SetTitle("FoM");
        gr_fom_th_online->SetTitle("FoM vs threshold (online)");

        gr_fom_th->SetMarkerStyle(20);
        gr_fom_th->SetMarkerColor(kOrange+10);
        gr_fom_th->GetXaxis()->SetTitle("qlong threshold");
        gr_fom_th->GetYaxis()->SetTitle("FoM");
        gr_fom_th->SetTitle("FoM vs threshold");

      }

      if (vm.count("calibration"))
      {
        vector<float> x = vm["calibration"].as<vector<float> >();

        h_psd_q->GetXaxis()->SetTitle("qlong (keVee)");
        h_psd_q->GetXaxis()->Set(psd_nbin,x[1],x[1]+psd_maxqlong*x[0]);

        h_psd_q_online->GetXaxis()->SetTitle("qlong (keVee)");
        h_psd_q_online->GetXaxis()->Set(psd_nbin,x[1],x[1]+psd_maxqlong*x[0]);

        h_qlong_online->GetXaxis()->SetTitle("qlong (keVee)");
        h_qlong_online->GetXaxis()->Set(psd_nbin,x[1],x[1]+psd_maxqlong*x[0]);

        h_qlong->GetXaxis()->SetTitle("qlong (keVee)");
        h_qlong->GetXaxis()->Set(psd_nbin,x[1],x[1]+psd_maxqlong*x[0]);
      }

      if (vm.count("gatepoints"))
      {
        ifstream p_file;

        p_file.open(gate_file.c_str());
        if (p_file.is_open())
        {
          vector<float> ql,psd;
          float tmp;
          int count =0;

          while(!p_file.eof())
          {
            p_file >> tmp;
            ql.push_back(tmp);
            p_file >> tmp;
            psd.push_back(tmp);
            count++;
          }

          count--;
          gr_gate = new TGraph(count);
          for (int i=0;i<count;i++)
          {
            gr_gate->SetPoint(i,ql[i],psd[i]);
          }
          gr_gate->SetMarkerStyle(20);
          gr_gate->SetMarkerColor(kOrange+10);
//          gr_gate->Fit("pol2");
            
        }else
        {
          cout << "Error opening gate_file\n";
          return 1;
        }
      }

      cnv_psd = new TCanvas("cnv_psd", "PSD", 0, 0, 800, 500);
      cnv_psd->cd();
      h_psd_q->Draw("COLZ");
      cnv_psd->SetLogz();
      cnv_psd->SaveAs(str_psd.c_str());

      cnv_psd_online = new TCanvas("cnv_psd_online", "PSD online", 0, 0, 800, 500);
      cnv_psd_online->cd();
      h_psd_q_online->Draw("COLZ");
      if (gr_gate) gr_gate->Draw("SAME");
      cnv_psd_online->SetLogz();
      cnv_psd_online->SaveAs(str_psd_online.c_str());

      cnv_fom_online = new TCanvas("cnv_fom_online", "FOM online", 0, 0, 800, 500);
      cnv_fom_online->cd();
      h_fom_online->Draw();
      cnv_fom_online->SaveAs(str_fom_online.c_str());

      cnv_fom = new TCanvas("cnv_fom", "FOM", 0, 0, 800, 500);
      cnv_fom->cd();
      h_fom->Draw();
      cnv_fom->SaveAs(str_fom.c_str());

      cnv_qlong_online = new TCanvas("cnv_qlong_online", "qlong_online", 0, 0, 800, 500);
      cnv_qlong_online->cd();
      h_qlong_online->Draw();

      cnv_qlong = new TCanvas("cnv_qlong", "qlong", 0, 0, 800, 500);
      cnv_qlong->cd();
      h_qlong->Draw();
    }

    if (gr_fom_th)
    {
      cnv_fom_th_online = new TCanvas("cnv_fom_th_online", "FoM vs threshold - online", 0, 0, 800, 500);
      cnv_fom_th_online->cd();
      gr_fom_th_online->Draw("ALP");
      cnv_fom_th_online->SaveAs(str_fom_th_online.c_str());

      cnv_fom_th = new TCanvas("cnv_fom_th", "FoM vs threshold", 0, 0, 800, 500);
      cnv_fom_th->cd();
      gr_fom_th->Draw("ALP");
      cnv_fom_th->SaveAs(str_fom_th.c_str());
    }

    if (gr_fom_pg)
    {
      cnv_fom_pg = new TCanvas("cnv_fom_pg", "FoM vs pregate", 0, 0, 800, 500);
      cnv_fom_pg->cd();
      gr_fom_pg->Draw("ALP");
      cnv_fom_pg->SaveAs(str_fom_pg.c_str());
    }

    if (gr_fom_sg)
    { 
      cnv_fom_sg = new TCanvas("cnv_fom_sg", "FoM vs shortgate", 0, 0, 800, 500);
      cnv_fom_sg->cd();
      gr_fom_sg->Draw("ALP");
      cnv_fom_sg->SaveAs(str_fom_sg.c_str());
    }

    if (gr_fom_lg)
    {
      cnv_fom_lg = new TCanvas("cnv_fom_lg", "FoM vs longgate", 0, 0, 800, 500);
      cnv_fom_lg->cd();
      gr_fom_lg->Draw("ALP");
      cnv_fom_lg->SaveAs(str_fom_lg.c_str());
    }

    if (h_scope)
    {
      cnv_scope = new TCanvas("cnv_scope","Scope", 0, 0, 800, 500);
      cnv_scope->cd();
      h_scope->Draw();
      cnv_scope->SaveAs(str_scope.c_str());
    }

    cnv_dummy->Close();

    if (!vm.count("non-interactive"))
    {
      cout << "Running Application\n";
      application->Run(kTRUE);
    }

    delete application;

    if (optim_file.is_open())
    {
      optim_file.close();
    }


  }
  catch(exception& e)
  {
    cout << e.what() << "\n";
    return 1;
  }

  return 0;

}
