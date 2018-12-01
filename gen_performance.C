void gen_performance() {

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);

  TCanvas *c_gen_perform = new TCanvas("c_gen_perform","c_gen_perform",600,400);

  double exp[10]={100,500,1000,5000,10000,50000,100000,500000,1000000,5000000};
  double cput[10]={0.68,0.72,0.78,1.11,1.5,5.2,9.6,47,89,490};
  double runt[10]={0.71,0.76,0.82,1.16,1.6,5.3,9.9,49,92,510};
  double fsize[10]={0.035156,0.169922,0.353516,1.3,2.8,16,28,189,249,1945.6};

  TPad *pad = new TPad("pad","",0,0,1,1);
  pad->SetFillColor(0);
  pad->Draw();
  pad->SetLogx();
  pad->SetLogy();
  pad->cd();

  TGraph *p_cpu=new TGraph(10,exp,cput);
  p_cpu->SetLineColor(kBlue);
  p_cpu->SetLineWidth(2);
  TGraph *p_run=new TGraph(10,exp,runt);
  p_run->SetLineColor(kRed);
  p_run->SetLineWidth(2);

  TMultiGraph *p_time=new TMultiGraph();
  p_time->Add(p_cpu);
  p_time->Add(p_run);
  p_time->Draw("AL");
  p_time->SetTitle("TRACKS generation - performances;# collisions;t [s]");
  p_time->GetXaxis()->SetLabelSize(0.045);
  p_time->GetYaxis()->SetLabelSize(0.045);
  p_time->GetXaxis()->SetTitleSize(0.05);
  p_time->GetYaxis()->SetTitleSize(0.05);
  p_time->GetYaxis()->SetTitleOffset(0.5);

  auto legt = new TLegend(0.15,0.65,0.3,0.85);
  legt->SetHeader("#varepsilon_{CPU} > 96%","");
  legt->AddEntry(p_cpu,"CPU time","l");
  legt->AddEntry(p_run,"RUN time","l");
  legt->Draw();

  TGraph *p_size=new TGraph(10,exp,fsize);
  p_size->SetLineColor(kGreen+2);
  p_size->SetLineWidth(2);
  p_size->Draw("SAMEL");

  //gPad->Update();

  TPad *overlay = new TPad("overlay","",0,0,1,1);
  overlay->SetFillStyle(4000);
  overlay->SetFillColor(0);
  overlay->SetFrameFillStyle(4000);
  overlay->Draw();
  overlay->cd();
  Double_t xmin = pad->GetUxmin();
  Double_t ymin = 0.25;
  Double_t xmax = pad->GetUxmax();
  Double_t ymax = pad->GetUymax();
  TH1F *hframe = overlay->DrawFrame(xmin,10,xmax,ymax);
  hframe->GetXaxis()->SetNdivisions(0);
  hframe->GetYaxis()->SetNdivisions(0);

  TGaxis *axis = new TGaxis(xmax,ymin,xmax,5.8,0.1,2000,505,"G+L");
  axis->SetLineColor(kGreen+2);
  axis->SetLabelColor(kGreen+2);
  axis->SetTitleColor(kGreen+2);
  axis->SetTitle("gen.root size [MB]");
  axis->SetLabelSize(0.045);
  axis->SetTitleSize(0.05);
  axis->Draw();

  MSaveBigPNG("zio", 2.);

  //c_gen_perform->SaveAs("c_gen_perform.eps");

}
