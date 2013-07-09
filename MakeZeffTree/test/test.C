void test(double mean, int denom, int cycles) {

  EfficiencyStatistics es(mean,denom);
  TH1* h1=new TH1F("Rand","Rand",102,0,1.00);
  TH1* h2=new TH1F("Lyt","Lyt",102,0,1.00);

  for (int i=0; i<cycles; ++i) {
    h1->Fill(es.rand());
  }

  for (int j=1; j<=h2->GetNbinsX(); ++j) {
    h2->Fill(h2->GetBinCenter(j),es.prob(h2->GetBinCenter(j)));
  }

  TCanvas* c1=new TCanvas("c1","c1");
  c1->Divide(1,2);
  c1->cd(1);
  h1->Draw();
  c1->cd(2);
  h2->Draw();

}
