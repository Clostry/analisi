#include "StructDefs.h"
#include "TSpectrum.h"

#define NSMEARINGS 200
#define NBINS 512
//#define NBINS 256
#define MAXHISTONRG 2000
#define HW_NBITS 16
#define MAXHISTOCHN (0x1<<HW_NBITS)
#define KN_NORM 5000
#define THRSHIFT 10
#define NUMENERGIES 2

char *branchname[8];

Char_t		smoothName[25];
TFile			*f_smearings;
TH1F			*h_ideal, *h_calib, *h_smooth;
Float_t		t_shift, max_calib, max_smooth;

Double_t f_profile(Double_t *x, Double_t *par) {
	Double_t	x_effective, value, norm;
	Int_t		bin;
	x_effective = x[0]-t_shift;
	bin = (h_smooth->FindBin(x_effective));
	norm = (max_calib-par[0])/max_smooth;
	value = norm*(h_smooth->GetBinContent(bin))+par[0]; 
	//printf("x[0]: %f\tx_eff: %f\tt_sh: %f\tbin: %d, MC: %f\tMS: %f\n",x[0],x_effective,t_shift,bin,max_calib,max_smooth);
	return value;
}

void calibraPlastico(char* filename, int channel=4) {

	FILE *outfile[2];
	outfile[0] = fopen("Chi2_511","w");
	outfile[1] = fopen("Chi2_1275","w");
	Int_t			number_of_loop=0;


	Int_t i, j, k, i_sm, rsen[2];
	Float_t r, alpha, energia;
	Int_t b_altezza;
			
	// energies
	Float_t	E_peak[NUMENERGIES];
	Float_t	E_compton[NUMENERGIES];
	E_peak[0] = 511.;
	E_peak[1] = 1275.;
	for (i=0; i<NUMENERGIES; i++) {
		E_compton[i] = 2*E_peak[i]*E_peak[i]/(511+2*E_peak[i]);
		printf("E_compton[%d] = %f;\n",i,E_compton[i]);
	}

  TTimer  *timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
	TCanvas *c0 = new TCanvas("c0");
	c0->cd();

	h_ideal = new TH1F("h_ideal","Compton ideale",NBINS,0,MAXHISTONRG);
	// check file existance
	f_smearings = new TFile("smearings.root","UPDATE");
	// check smearing samples existance
	for(i=0; i<NUMENERGIES; i++) {
		sprintf(smoothName,"smooth_%.1f_%d;1",E_peak[i],NSMEARINGS-2); 
		if (!f_smearings->Get(smoothName)) {
			cout << smoothName << " " << f_smearings->FindObject(smoothName) << endl;
			// smearing for that energy do not exist
			printf("Non esistono.\n");
			// ideal compton histogram
			for (j=0; j<NBINS; j++) {
				if (j>h_ideal->FindBin(50) && j<h_ideal->FindBin(E_compton[i])) {
					r = h_ideal->GetBinCenter(j)/E_peak[i];
					alpha = E_peak[i]/511.0;
					energia = KN_NORM * (2+r*r/(alpha*alpha*(1-r)*(1-r))+r/(1-r) * (r-2/alpha));
					h_ideal->SetBinContent(j,energia);
				} else { h_ideal->SetBinContent(j,0); }
			}
			h_ideal->Draw();
      c0->Update();
      timer->TurnOn();
      timer->Reset();
      timer->TurnOff();
			// creazione spettri smussati
			for (i_sm=1; i_sm<NSMEARINGS; i_sm++) {
				sprintf(smoothName,"smooth_%.1f_%d",E_peak[i],i_sm); 
				printf("Creo '%s': ",smoothName);
				h_smooth = new TH1F(smoothName,"Smooth",NBINS,0,MAXHISTONRG); // istogramma in energia
				for (j=1; j<h_ideal->FindBin(E_compton[i]); j++) {
					b_altezza = h_ideal->GetBinContent(j);
					for (k=1; k<b_altezza; k++){ h_smooth->Fill(gRandom->Gaus(h_ideal->GetBinCenter(j),i_sm)); } printf(".");
				} printf("\n");
				h_smooth->Write();
			}
		}
	}
	// ok, we've got the smearings!	

	f_smearings->Close();
	f_smearings = new TFile("smearings.root","READ");


	// ----------------------------------
	// retrieving "raw" histogram
	// ----------------------------------
	TFile *infile = new TFile(filename);
	TTree *tree= (TTree*)infile->Get("acq_tree_0");
	TBranch *branch = tree->GetBranch(Form("acq_ch%d",channel));
	branch->SetAddress(&inc_data.timetag);	

	TH1F *h_raw = new TH1F("h_raw","Acquisizione",NBINS,0,MAXHISTOCHN);
	UInt_t toentry=branch->GetEntries();
	printf("getHistoFromFile: There are %d entries in the branch\n",toentry);
	for(i=0; i<toentry; i++) {
		branch->GetEntry(i);
		h_raw->Fill(inc_data.qlong);
	}

  h_raw->Draw();
    
	TSpectrum	*s = new TSpectrum(10);
	Int_t		e, nPeaks, bTemp, bFirstPeak = 9999;
	Float_t	*xPeaks;
    
	// trovo il primo picco
	nPeaks = s->Search(h_raw->Rebin(2, "h_raw_rebinned"));
	if (nPeaks>0) {
		xPeaks = s->GetPositionX();
		// loop sui picchi per trovare il primo
		for (i=0;i<nPeaks;i++) {
			bTemp = h_raw->GetXaxis()->FindBin(xPeaks[i]);
			if (bTemp<bFirstPeak) { bFirstPeak = bTemp; }
		}
	} else { bFirstPeak = 0; }
    // sottraggo il fondo Compton
	Float_t 		*bgBins = new Float_t[NBINS];
	TH1F			*h_filtered = new TH1F("h_filtered","Picchi",NBINS,0,MAXHISTOCHN);
	/*for (i = 0; i < bFirstPeak; i++) { bgBins[i]=h_raw->GetBinContent(i+1); }
	s->Background(bgBins,bFirstPeak,20,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kFALSE,TSpectrum::kBackSmoothing15,kFALSE);
	for (i = 0; i < bFirstPeak; i++) { h_filtered->SetBinContent(i+1, (h_raw->GetBinContent(i+1)-bgBins[i])); }
	for (i = 0; i < NBINS-bFirstPeak; i++) { bgBins[i]=h_raw->GetBinContent(bFirstPeak+i+1); }
	s->Background(bgBins,NBINS-bFirstPeak,20,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,kFALSE,TSpectrum::kBackSmoothing15,kFALSE);
	for (i = bFirstPeak; i < NBINS; i++) { h_filtered->SetBinContent(i+1, (h_raw->GetBinContent(i+1)-bgBins[i-bFirstPeak])); }
	//TCanvas * background = new TCanvas("background","Estimation of bg",10,10,1000,700);*/
	
  for (i = 0; i < NBINS; i++) { bgBins[i]=h_raw->GetBinContent(i+1); }

	s->Background(bgBins, NBINS, 20, TSpectrum::kBackDecreasingWindow, TSpectrum::kBackOrder2, kFALSE, TSpectrum::kBackSmoothing15, kFALSE);

	for (i = 0; i < NBINS; i++) { h_filtered->SetBinContent(i+1, (h_raw->GetBinContent(i+1)-bgBins[i])); }
	
  h_raw->Draw("L");
	h_filtered->SetLineColor(kRed);
	h_filtered->Draw("SAME L");
  c0->Update();
		// trovo i picchi e calibro
    
    
	nPeaks = s->Search(h_filtered->Rebin(2, "h_filtered_rebinned"),2,"",0.05);
    xPeaks = s->GetPositionX();
	//	sigma_fall[0] = 100;


  timer->TurnOn();
  timer->Reset();
  timer->TurnOff();



	if (nPeaks<NUMENERGIES) {
		// trovati troppi pochi picchi
		printf("EPIC FAIL while calibrating - too few peaks!\n");
	} else {
		// possiamo calibrare
		TGraphErrors	*graphErr;
		TF1				*fitMeasPeaks, *fitfun;
		Float_t			nrg[NUMENERGIES], shift[NUMENERGIES], sigma_fall[NUMENERGIES], a, b, bPeak, cPeak, nrgPeak, chi2, fondo, xmin, xmax, chi2min;
		// fitto i due picchi, mi servono le sigma...
		for (e=0; e<NUMENERGIES; e++) {
			xmin = xPeaks[e] - (0.5*xPeaks[e]);	// fit left margin
			xmax = xPeaks[e] + (0.5*xPeaks[e]);	// fit right margin
			fitMeasPeaks = new TF1("fitMeasPeaks","gaus",xmin,xmax);
			fitMeasPeaks->SetNpx(1000);
			fitMeasPeaks->SetParameters(1,xPeaks[e]);
			h_filtered->Fit("fitMeasPeaks","QNO","",xmin,xmax);
			sigma_fall[e] = fitMeasPeaks->GetParameter(2);
			printf("%f - %f - %f Sigma_fall[%d]: %f\n",xmin,xPeaks[e],xmax,e,sigma_fall[e]);
			// inizializzazione
			shift[e]=0; nrg[e]=E_compton[e];
		}
		fitfun = new TF1("calfitfun","pol1",xPeaks[0],xPeaks[NUMENERGIES-1]);
		// costruiamo uno spettro sperimentale calibrato
		h_calib = (TH1F*)h_raw->Clone();
		h_calib->SetNameTitle("h_calib","Acquisiz. calibrata");
		Int_t debug = 0;
		Int_t	loop_flag = 1;
		// do..while shift begin
		while(loop_flag) {
			number_of_loop++;
			loop_flag = 0;
			for (e=0; e<NUMENERGIES; e++) {
				// energie aggiornate
				nrg[e] = nrg[e]-shift[e];
				printf("xPeaks[%d] = %f\tshift[%d] = %f\tnrg[%d] = %f\n",e,xPeaks[e],e,shift[e],e,nrg[e]);
			}
			// calibrazione
			fitfun->SetParLimits(0,-100,100); fitfun->SetParLimits(1,0.,1.);
			graphErr = new TGraphErrors(NUMENERGIES,xPeaks,nrg);	graphErr->Fit(fitfun,"RN");
			a = fitfun->GetParameter(1); b = fitfun->GetParameter(0);// graphErr->Delete();
			printf("intercetta=%f pendenza=%f\n",b,a);
			h_calib->GetXaxis()->SetLimits(b,h_raw->GetBinCenter(NBINS)*a+b);
			h_calib->Draw();
      c0->Update();
      timer->TurnOn();
      timer->Reset();
      timer->TurnOff();
			for (e=0; e<NUMENERGIES; e++) {
				chi2min = 9999999;
				printf("Looppo sull'energia %.1f\n",nrg[e]);
				// loop on smearings
				for (i_sm=19; i_sm<NSMEARINGS; i_sm++) {
					sprintf(smoothName,"smooth_%.1f_%d",E_peak[e],i_sm);
					h_smooth = (TH1F*)f_smearings->Get(smoothName);				if (debug) printf("Recupero l'histo %s\n",smoothName);
					//h_smooth->Draw();
					// momentaneamente assumiamo che il CE sia:
					bPeak = h_smooth->GetMaximumBin();								if (debug) printf("bPeak = %f\n",bPeak);
					nrgPeak = bPeak*MAXHISTONRG/NBINS;								if (debug) printf("nrgPeak = %f\n",nrgPeak);
					t_shift = nrg[e]-nrgPeak;											if (debug) printf("t_shift = %f\n",t_shift);
					cPeak = h_calib->FindBin(a*xPeaks[e]+b);						if (debug) printf("peak calib: %.1f\n",cPeak);
					max_calib = h_calib->GetBinContent(cPeak);					if (debug) printf("max_calib = %f\n",max_calib);
					max_smooth = h_smooth->GetBinContent(bPeak);					if (debug) printf("max_smooth = %f\n",max_smooth);
					if (debug) printf("sigma_fall[%d] = %f\n",e,sigma_fall[e]);	if (debug) printf("a*sigma_fall[%d] = %f\n",e,a*sigma_fall[e]);
					// definisco la funzione per il fit
					TF1 *f_smear_profile = new TF1("f_smear_profile",f_profile,nrgPeak-a*sigma_fall[e],nrgPeak+a*sigma_fall[e]*3,1);
					f_smear_profile->SetParameters(0,300);
					h_calib->Fit("f_smear_profile","QON","",a*xPeaks[e]+b-a*sigma_fall[e],a*xPeaks[e]+b+a*sigma_fall[e]*3);
					chi2 = f_smear_profile->GetChisquare();						if (debug) printf("i_sm: %i\tchiquadro: %f\n",i_sm,chi2);
					fondo = f_smear_profile->GetParameter(0);						if (debug) printf("fondo = %f\n",fondo);
					if (debug) printf("i_sm: %d\tchi2min vs. chi2: %f %f\n",i_sm,chi2min,chi2);
					if (number_of_loop==1) {fprintf(outfile[e],"%i\t%f\t%f\n",i_sm,chi2, fondo);}
					if (chi2<chi2min) { chi2min = chi2; shift[e] = t_shift; rsen[e]=i_sm;}
				} // loop on smearings ends here
				if (debug) printf("shift[%d]: %f\n",e,t_shift);
				if (shift[e]>THRSHIFT) { loop_flag++; }
				if (number_of_loop==1) {fclose(outfile[e]);}
			}	// end loop on energies
		} // while shift ends here
		printf("Esco con shift[0]=%.1f e shift[1]=%.1f\n",shift[0],shift[1]);
		cout << "Esco con risoluzione[0] = " << rsen[0] << " keV (" << rsen[0]*100/E_compton[0] << "%)";
		cout << " e risoluzione[1] = " << rsen[1] << " keV (" << rsen[1]*100/E_compton[1] << "%)" << endl;
		printf("\nRESULTS: nrg1 = %f\tnrg2 = %f\n",nrg[0],nrg[1]);
		printf("RESULTS:    m = %f\t   q = %f\n",a,b);
	} // end if enough peaks
    infile->Close();
	//f_smearings->Close();
}
