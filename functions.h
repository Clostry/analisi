//functions.h
//
//

int trigger(event_data data, acqPSDParam_t params);

bool saturation(event_data data, acqPSDParam_t params, int *lastbin, int *nbin );

bool saturation(event_data data, acqPSDParam_t params );

void GetIntegralFromScope(event_data inc_data, acqPSDParam_t psd_params, float *qlong, float *qshort, float *baseline, bool baseline_from_scope = true );
 
bool pileup(event_data data, acqPSDParam_t params,int threshold);

void fill_psd_q(TH2F *h, DataInterface *idata, acqPSDParam_t psd_params, std::vector<int> th, bool get_from_scope);
