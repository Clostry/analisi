void effgraph(char *filename)
{
  ifstream file;
  TCanvas *cnv;
  float readings[5][200];
  int energies[200];
  TGraph *gr[5];
  float tmp;

/*
  float m = 0.2085;
  float q = 7.16;
*/
  float m = 0.254822;
  float q = -0.585938;

  cnv = new TCanvas("cnv","Efficienza",0,0,800,500);
  file.open(filename);
  int lines=0;
  file >> tmp;

  while(!file.eof())
  {
    energies[lines] = tmp;
    for(int i=0;i<8;i++)
    {
      file >> tmp;
    }
    for(int i=0;i<5;i++)
    {
      file >> tmp;
      readings[i][lines] = tmp;
    }
    lines++;
    file >> tmp;
  }

  file.close();
  cout << "file closed" << endl;

  cnv->cd();
  for(int i=0;i<5;i++)
  {
    cout << "doing graph " << i+1 << endl;
    gr[i]=new TGraph(lines);
    for(int j=0;j<lines;j++)
    {
      cout << j << endl;
      gr[i]->SetPoint(j,energies[j]*m+q,readings[i][j]);
    }
    cout << "got points\n";
    gr[i]->SetLineColor(5-i);
  }
  TLegend *leg = new TLegend(0.7,0.3,0.9,0.1);
  leg->AddEntry(gr[0],"0.5%","l");
  leg->AddEntry(gr[1],"1%","l");
  leg->AddEntry(gr[2],"2%","l");
  leg->AddEntry(gr[3],"5%","l");
  leg->AddEntry(gr[4],"10%","l");

  gr[0]->SetTitle("Efficienza");
  gr[0]->GetXaxis()->SetTitle("keVee");
  gr[0]->Draw();
  gr[1]->Draw("SAME");
  gr[2]->Draw("SAME");
  gr[3]->Draw("SAME");
  gr[4]->Draw("SAME");
  leg->Draw("SAME");

  float min = gr[0]->GetXaxis()->GetXmin();
  float max = gr[0]->GetXaxis()->GetXmax();
  TF1 *f = new TF1("f","[0]",min, max);
  f->SetLineColor(2);
  f->SetLineWidth(1);
  f->SetParameter(0,1);
  f->Draw("SAME");
  
}
