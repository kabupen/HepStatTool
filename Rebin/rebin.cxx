#include "HistoTransform.h"

void rebin() 
{
    TFile* infile = TFile::Open("/Users/takeda/cernbox/merged.squirtle.21Feb2020.mc16e_Syst.root", "OPEN");

    for ( auto btag : {"1tag", "2tag"} ) {
        for ( auto mass : {"300", "1000","1200"} ) {

            /* Background lists */
            std::map<std::string, TH1D*> bkg_map;
            for ( auto bkg : {"DY", "DYtt", "VH", "W", "WZ", "WW", "Wtt", "Zbc", "Zl", "Zcc", "Zcl", "Zbl", "Zbb", "ZZ", "Zttcl", "Zttbb", "Zttl", "Zttbc", "Zttcc", "Zttbl", "data", "stopWt", "stops", "stopt", "ttH", "ttbar"} ) {

                if ( infile->Get( Form("BasicKinematics_FullRun2/%s_%s2pjet_0ptv_TauPT", bkg, btag) ) == NULL ) {
                    std::cout << "@@@@ Not found : BasicKinematics_FullRun2/" << bkg << "_" << btag << "2pjet_0ptv_TauPT" << std::endl;
                    continue;
                }
                bkg_map[bkg] = (TH1D*)infile->Get( Form("BasicKinematics_FullRun2/%s_%s2pjet_0ptv_TauPT", bkg, btag) );
            }

            /* Signal */
            TH1D* hSig = (TH1D*)infile->Get( Form("BasicKinematics_FullRun2/LQ3Up%s_%s2pjet_0ptv_TauPT", mass, btag) );
            TH1D* hBkg = nullptr;
            for ( auto itr = bkg_map.begin(); itr != bkg_map.end(); itr++ ){
                if ( itr == bkg_map.begin() ) hBkg = itr->second;
                else                          hBkg->Add(itr->second);
            }

            // create transformation
            HistoTransform histoTrafo;
            histoTrafo.trafoFzSig = 0.3;
            histoTrafo.trafoFzBkg = 4;
            int method = 12; // 12 = "Transformation F"
            float maxUnc = 1;
            std::vector<int> bins = histoTrafo.getRebinBins(hBkg, hSig, method, maxUnc);

            // Output as a root file
            TFile* outfile = new TFile( Form("rebin_%s_%sGeV.root", btag, mass), "RECREATE");
            outfile->mkdir("Before");
            outfile->cd("Before");

            hSig->SetName(Form("LQ3Up%s", mass));
            hSig->Write();

            for ( auto itr = bkg_map.begin(); itr != bkg_map.end(); itr++ ){
                itr->second->SetName(Form("%s", itr->first.c_str()));
                itr->second->Write();
            }

            // transform any histogram
            outfile->cd();
            for ( auto itr = bkg_map.begin(); itr != bkg_map.end(); itr++ ){
                histoTrafo.rebinHisto(itr->second, &bins);
                itr->second-> Write();
            }
            histoTrafo.rebinHisto(hSig, &bins);
            hSig -> Write();
            
            outfile->Write();
            outfile->Close();
            std::cout << outfile->GetName() << std::endl;
        }
    }
}
