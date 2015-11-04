#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <sstream>

struct event_data 
{
  ULong64_t timetag;
  UInt_t    baseline;
  UShort_t  qshort;
  UShort_t  qlong;
  UShort_t  pur;
  UShort_t  samples[1024];
}; 

struct acqEventPSD_old_t 
{
  UInt_t timetag;
  UInt_t qshort;
  UInt_t qlong;
  UInt_t baseline;
  UInt_t samples[100];
}; 

struct acqPSDParam_t 
{
  UInt_t channel;
  UInt_t threshold;
  UInt_t pretrigger;
  UInt_t pregate;
  UInt_t shortgate;
  UInt_t longgate;
  UInt_t numsamples;
};

class DataInterface {
  public:
    DataInterface(const char *filename);
    ~DataInterface();
    int SetChannel(int channel);
    int GetEntries();
    event_data GetEntry(int n);
    acqPSDParam_t GetParams();

  private:
    TFile *f;
    TTree *t;
    TBranch *b;
    int entries;
    event_data entry;
    acqPSDParam_t params;
};  
