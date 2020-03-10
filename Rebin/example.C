#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "HistoTransform.h"

//int main(int argc, char* argv[]) {
void example() {

  // create example histograms
  TH1* hBkg1 = new TH1F("hBkg1","hBkg1", 1000, -10., 10.);
  TH1* hBkg2 = new TH1F("hBkg2","hBkg2", 1000, -10., 10.);
  TH1* hBkg3 = new TH1F("hBkg3","hBkg3", 1000, -10., 10.);
  TH1* hBkg  = new TH1F("hBkg","hBkg", 1000, -10., 10.);
  TH1* hSig  = new TH1F("hSig","hSig", 1000, -10., 10.);
  hBkg1 -> Sumw2();
  hBkg2 -> Sumw2();
  hBkg3 -> Sumw2();
  hBkg  -> Sumw2();
  hSig  -> Sumw2();

  hBkg1 -> FillRandom("pol0", 10000);
  hBkg2 -> FillRandom("pol0", 10);
  hBkg3 -> FillRandom("pol0", 100);
  hBkg->Add(hBkg1);
  hBkg->Add(hBkg2);
  hBkg->Add(hBkg3);
  hSig  -> FillRandom("gaus", 10000);
  
  hBkg -> Scale(1./100);
  hSig -> Scale(1./1000);

  // create transformation
  HistoTransform histoTrafo;
  histoTrafo.trafoFzSig = 0.3;
  histoTrafo.trafoFzBkg = 4;
  int method = 12; // 12 = "Transformation F"
  float maxUnc = 1;
  vector<int> bins = histoTrafo.getRebinBins(hBkg, hSig, method, maxUnc);

  // draw original
  TCanvas* c1 = new TCanvas();
  hSig -> SetLineColor(kRed);
  hBkg -> DrawCopy("");
  hSig -> DrawCopy("same");

  // transform any histogram
  histoTrafo.rebinHisto(hBkg1, &bins);
  histoTrafo.rebinHisto(hBkg2, &bins);
  histoTrafo.rebinHisto(hBkg3, &bins);
  histoTrafo.rebinHisto(hSig, &bins);

  // draw transformed
  TCanvas* c2 = new TCanvas();
  hBkg1 -> GetYaxis() -> SetRangeUser(0, 40);
  hSig -> Scale(10.); // scale for visibility

  hBkg1 -> SetLineColor(kBlack);
  hBkg1 -> Scale(1/100.);
  hBkg1 -> Draw("");

  hBkg2 -> SetLineColor(kMagenta);
  hBkg2 -> Draw("same");

  hBkg3 -> SetLineColor(kGreen);
  hBkg3 -> Draw("same");

  hSig  -> Draw("same");
}
