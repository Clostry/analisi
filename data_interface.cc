// data_interface.cc
#include <data_interface.h>

DataInterface::DataInterface(const char *filename)
{
  f = new TFile(filename);
  t = (TTree*)f->Get("acq_tree_0");
  b = 0;
  SetChannel(0);

}


DataInterface::~DataInterface()
{
  delete b;
  delete t;
  delete f;
}


int DataInterface::SetChannel(int channel)
{
  if ((channel<0)||(channel>7))
  {
    return 1;
  }

  b = t->GetBranch("psd_params");
  b->SetAddress(&params);
  b->GetEntry(channel);
//  std::string tmp = "acq_ch" + std::to_string(channel);
//  b = t->GetBranch(tmp.c_str());
  b = t->GetBranch(Form("acq_ch%i",channel));
  b->SetAddress(&entry);
  entries = b->GetEntries();

  return 0;
}

int DataInterface::GetEntries()
{
  return entries;
}

event_data DataInterface::GetEntry(int n)
{
  b->GetEntry(n);
  return entry;
}

acqPSDParam_t DataInterface::GetParams()
{
  return params;
}
