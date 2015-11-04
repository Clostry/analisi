#include <data_interface.h>
#include "root_include.h"

using namespace std;

int trigger(event_data data, acqPSDParam_t params)
{
  for(int i=0;i<(int)params.numsamples;i++)
  {
    if ((data.baseline - data.samples[i]) > params.threshold)
    {
      return i;
    }
  }
  return -1;
}

bool saturation(event_data data, acqPSDParam_t params, int *lastbin , int *nbin )
{

  int max = params.pretrigger-params.pregate+params.longgate;
  //non andare oltre il numero di sample
  if (max > (int)params.numsamples)
    max = params.numsamples;

  if (nbin) *nbin = 0;
  for( int i = 0; i < max; i++ )
  {
    if (data.samples[i] == 0 )
    {
      if (lastbin) *lastbin = i;
      if (nbin) (*nbin)++;
      return true;
    }
  }

  return false;

}

bool saturation(event_data data, acqPSDParam_t params)
{

  int max = params.pretrigger-params.pregate+params.longgate;
  //non andare oltre il numero di sample
  if (max > (int)params.numsamples)
    max = params.numsamples;

  for( int i = 0; i < max; i++ )
  {
    if (data.samples[i] == 0 )
    {
      return true;
    }
  }

  return false;

}



void GetIntegralFromScope(event_data inc_data, acqPSDParam_t psd_params, float *qlong, float *qshort, float *baseline, bool baseline_from_scope = true )
{

  float integral = 0;
  int max = 0;
  int j;
  if (baseline_from_scope)
  {
    //calcola baseline 
    *baseline = 0;
    for( j=0; j < (int) (psd_params.pretrigger-psd_params.pregate); j++)
    {
      *baseline += inc_data.samples[j];
    }
    //fai la media
    *baseline /= psd_params.pretrigger-psd_params.pregate;
  } else
  { 
    // usa baseline scheda
    *baseline = inc_data.baseline;
  }
  max = psd_params.pretrigger-psd_params.pregate+psd_params.shortgate+1;
  //non andare oltre il numero di sample
  if (max > (int)psd_params.numsamples)
    max = psd_params.numsamples;
  //calcola integrale qshort
  for( j = psd_params.pretrigger - psd_params.pregate + 1; j < max; j++) 
  {
    integral += *baseline-inc_data.samples[j];
  }
  *qshort = integral;

  max = psd_params.pretrigger-psd_params.pregate+psd_params.longgate + 1;
  //non andare oltre il numero di sample
  if (max > (int)psd_params.numsamples) 
    max = psd_params.numsamples;
  //calcola integrale qlong
  for( j = psd_params.pretrigger - psd_params.pregate + psd_params.shortgate + 1 ; j < max; j++ )
  {
    integral += *baseline-inc_data.samples[j];
  }
  *qlong = integral;

}

bool pileup(event_data data, acqPSDParam_t params, int threshold)
{
  int maxbin;
  int max;
  maxbin = params.pretrigger - params.pregate+params.longgate + 1;
  if (maxbin > (int)params.numsamples)
        maxbin = params.numsamples;

  int min;
  int minbin = params.pretrigger - 2;

  do
  {

    min = data.samples[minbin];
    minbin++;

  }while((data.samples[minbin] < min) && (minbin < maxbin));

  max = data.samples[minbin];

  for( int i = minbin + 1; i < maxbin; i++)
  {
    if ( (max - data.samples[i])  > threshold )
    {
      return true;
    }
    if (max < data.samples[i])
    {
      max = data.samples[i];
    }

  }

  return false;
}

/*
bool pileup(event_data data, acqPSDParam_t params, int threshold)
{
  int maxbin;
  maxbin = params.pretrigger - params.pregate+params.longgate + 1;
  if (maxbin > (int)params.numsamples)
        maxbin = params.numsamples;

  int min;
  int minbin = params.pretrigger - 2;

  do
  {

    min = data.samples[minbin];
    minbin++;

  }while((data.samples[minbin] < min) && (minbin < maxbin));

  for( int i = minbin + 1; i < maxbin; i++)
  {
    if ( (data.samples[i-1] - data.samples[i])  > threshold )
    {
      return true;
    }

    if ( (data.samples[i-2] - data.samples[i])  > threshold )
    {
      return true;
    }

  }

  return false;
}*/


void fill_psd_q(TH2F *h, DataInterface *idata, acqPSDParam_t psd_params, vector<int> th, bool get_from_scope)
{
  int q_scaling = 0; //scaling da implementare
  float psd;
  float baseline;
  int n = idata->GetEntries();
  float qlong;
  float qshort;

  if (get_from_scope)
  {

    for (int i = 0; i < n; i++)
    {
      event_data inc_data = idata->GetEntry(i);
      //if (!saturation(inc_data,psd_params))
      //if (!pileup(inc_data,psd_params)
      {
        GetIntegralFromScope(inc_data, psd_params, &qlong, &qshort, &baseline, false );
        if ( (qlong > th[0])&&((th[1]==0)||(qlong < th[1])))
        {
          psd = (qlong - qshort) / qlong; 
          h->Fill( ( ((int)qlong) >> q_scaling ),psd );
        }
      }
    }
  }
  else
  {
    for (int i = 0; i < n; i++)
    {
      event_data inc_data = idata->GetEntry(i);
      //if (( inc_data.qlong > th )&&( !pileup(inc_data,psd_params)))
      if (( inc_data.qlong > th[0] )&&((th[1]==0)||(inc_data.qlong < th[1])))
      {
        psd = (inc_data.qlong - inc_data.qshort) / (float) inc_data.qlong;
        h->Fill( ( ((int)inc_data.qlong) >> q_scaling ),psd );
      }
    }
  }

  return;

}
