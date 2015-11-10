#include "../data_interface.h"
#include "../data_interface.cc"

void dump(char *filename, int canale=4)
{
	DataInterface *idata;
	acqPSDParam_t params;
	event_data evento;

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
}
