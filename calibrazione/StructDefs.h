#ifndef _STRUCTDEF_H_
#define _STRUCTDEF_H_

#include <TGraph.h>
#include <TBranch.h>
#include <TSpectrum.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TRandom1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

/*struct acqEventSTD_t {
	ULong64_t	timetag;
	UInt_t	counter;
	UShort_t	samples[4096];
} input_data;*/

struct acqEventSTD_t {
	ULong64_t	timetag;
	UInt_t	    baseline;
	UShort_t	qshort;
	UShort_t	qlong;
	UShort_t	pur;
	UShort_t	samples[4096];
} input_data;

struct acqEventPSD_t {
	ULong64_t	timetag;
	UInt_t		baseline;
	UShort_t		qshort;
	UShort_t		qlong;
	UShort_t		pur;
	UShort_t		samples[4096];
} inc_data, output_data;

/*struct acqEventPSD_t {
	ULong64_t	timetag;
	UInt_t		baseline;
	UInt_t      qshort;
	UInt_t  	qlong;
	UShort_t  	pur;
	UShort_t	samples[4096];
} output_data, inc_data;*/

struct acqParams_t {
	UInt_t      acqTimestamp;
	UShort_t	dppVersion;
	UShort_t	nsPerSample;
	UShort_t	nsPerTimetag;
	UShort_t	numChannels;
	UShort_t	numSamples[32];
} inc_params;

struct acqParamsCal_t {
	UInt_t      acqTimestamp;
    Float_t     cal_m;
    Float_t     cal_q;
	UShort_t	dppVersion;
	UShort_t	nsPerSample;
	UShort_t	nsPerTimetag;
	UShort_t	numChannels;
	UShort_t	numSamples[32];
} out_params, params;

#endif