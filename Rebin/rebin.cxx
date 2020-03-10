#include "HistoTransform.h"

std::vector<std::string> masses = 
{
    "1000",
    "1200"
};

std::vector<std::string> backgrounds = 
{
    "DY", "DYtt", "VH", "W", "WZ", "WW", "Wtt", "Zbc", "Zl", "Zcc", "Zcl", "Zbl", "Zbb", "ZZ", "Zttcl", "Zttbb", "Zttl", "Zttbc", "Zttcc", "Zttbl", "data", "stopWt", "stops", "stopt", "ttH", "ttbar"
}; 

std::vector<std::string> systematics = 
{
    "SysFT_EFF_extrapolation_from_charm__1up",
    "SysTTbarST__1up",
    "SysFT_EFF_Eigen_C_0__1up",
    "SysFT_EFF_Eigen_Light_2__1down",
    "SysTTbarST__1down",
    "SysFT_EFF_Eigen_Light_4__1down",
    "SysTAUS_TRUEHADTAU_EFF_RECO_HIGHPT__1down",
    "SysFT_EFF_Eigen_Light_0__1down",
    "SysFT_EFF_Eigen_Light_3__1up",
    "SysTAUS_TRUEHADTAU_EFF_JETID_HIGHPT__1down",
    "SysFT_EFF_Eigen_C_1__1up",
    "SysFT_EFF_Eigen_B_0__1down",
    "SysTTbarMBB__1up",
    "SysFT_EFF_Eigen_Light_1__1down",
    "SysTAUS_TRUEHADTAU_EFF_RECO_TOTAL__1down",
    "SysFT_EFF_Eigen_C_2__1up",
    "SysFT_EFF_Eigen_B_2__1up",
    "SysFT_EFF_extrapolation__1up",
    "SysTAUS_TRUEHADTAU_EFF_JETID_HIGHPT__1up",
    "SysFT_EFF_Eigen_B_0__1up",
    "SysFT_EFF_Eigen_B_1__1up",
    "SysFT_EFF_Eigen_B_2__1down",
    "SysFT_EFF_Eigen_C_0__1down",
    "SysFT_EFF_Eigen_C_1__1down",
    "SysFT_EFF_Eigen_Light_0__1up",
    "SysFT_EFF_Eigen_Light_1__1up",
    "SysFT_EFF_Eigen_Light_2__1up",
    "SysFT_EFF_Eigen_Light_3__1down",
    "SysFT_EFF_Eigen_Light_4__1up",
    "SysFT_EFF_extrapolation__1down",
    "SysPRW_DATASF__1up",
    "SysPRW_DATASF__1down",
    "SysTAUS_TRUEHADTAU_EFF_ELEOLR_TOTAL__1up",
    "SysTAUS_TRUEHADTAU_EFF_ELEOLR_TOTAL__1down",
    "SysTAUS_TRUEHADTAU_EFF_JETID_TOTAL__1up",
    "SysTAUS_TRUEHADTAU_EFF_JETID_TOTAL__1down",
    "SysTAUS_TRUEHADTAU_EFF_RECO_TOTAL__1up",
    "SysTAUS_TRUEHADTAU_EFF_RECO_HIGHPT__1up",
    "SysFT_EFF_Eigen_B_1__1down",
    "SysTTbarMBB__1down",
    "SysFT_EFF_Eigen_C_2__1down",
    "SysFT_EFF_extrapolation_from_charm__1down"
};


void rebin() 
{
    TFile* infile = TFile::Open("/Users/takeda/cernbox/merged.squirtle.21Feb2020.mc16e_Syst.root", "OPEN");

    for ( auto btag : {"1tag", "2tag"} ) {
        for ( auto mass : masses ) {

            /* Background lists */
            std::map<std::string, TH1D*> bkg_map;
            for ( auto bkg : backgrounds ) {
                if ( infile->Get( Form("BasicKinematics_FullRun2/%s_%s2pjet_0ptv_TauPT", bkg.c_str(), btag) ) == NULL ) {
                    std::cout << "@@@@ Not found : BasicKinematics_FullRun2/" << bkg << "_" << btag << "2pjet_0ptv_TauPT" << std::endl;
                    continue;
                }
                bkg_map[bkg] = (TH1D*)infile->Get( Form("BasicKinematics_FullRun2/%s_%s2pjet_0ptv_TauPT", bkg.c_str(), btag) );
            }

            /* systematics */
            std::map<std::string, std::vector<TH1D*>> syst_map;
            for ( auto syst : systematics ){
                for ( auto bkg : backgrounds ) { 
                    if ( infile->Get( Form("BasicKinematics_FullRun2/Systematics/%s_%s2pjet_0ptv_TauPT_%s", bkg.c_str(), btag, syst.c_str()) ) == NULL ) {
                        std::cout << "@@@@ Not found : BasicKinematics_FullRun2/" << bkg << "_" << btag << "2pjet_0ptv_TauPT_" << syst << std::endl;
                        continue;
                    }
                    TH1D* htmp = (TH1D*)infile->Get( Form("BasicKinematics_FullRun2/Systematics/%s_%s2pjet_0ptv_TauPT_%s", bkg.c_str(), btag, syst.c_str()) );
                    htmp->SetName( Form("%s_%s", bkg.c_str(), syst.c_str()) );
                    syst_map[syst].push_back(htmp);
                }
                if ( infile->Get( Form("BasicKinematics_FullRun2/Systematics/LQ3Up%s_%s2pjet_0ptv_TauPT_%s", mass.c_str(), btag, syst.c_str()) ) == NULL ) {
                    continue;
                }
                TH1D* htmp = (TH1D*)infile->Get( Form("BasicKinematics_FullRun2/Systematics/LQ3Up%s_%s2pjet_0ptv_TauPT_%s", mass.c_str(), btag, syst.c_str()));
                htmp->SetName( Form("LQ3Up%s_%s", mass.c_str(), syst.c_str()) );
                syst_map[syst].push_back(htmp);
            }


            /* Signal */
            TH1D* hSig = (TH1D*)infile->Get( Form("BasicKinematics_FullRun2/LQ3Up%s_%s2pjet_0ptv_TauPT", mass.c_str(), btag) );
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
            TFile* outfile = new TFile( Form("rebin_%s_%sGeV.root", btag, mass.c_str()), "RECREATE");
            outfile->mkdir("Before");
            outfile->cd("Before");

            // Write the original hists
            hSig->SetName(Form("LQ3Up%s", mass.c_str()));
            hSig->Write();
            // bkg
            for ( auto itr = bkg_map.begin(); itr != bkg_map.end(); itr++ ){
                itr->second->SetName(Form("%s", itr->first.c_str()));
                itr->second->Write();
            }

            // transform systematics
            for ( auto syst : systematics ) {
                outfile->mkdir(syst.c_str());
                outfile->cd(syst.c_str());
                
                // Write the original hists
                outfile->mkdir(Form("%s/Before",syst.c_str()));
                outfile->cd(Form("%s/Before",syst.c_str()));
                for ( auto sys_hist : syst_map[syst] ) {
                    sys_hist->Write();
                }
                outfile->cd();

                // Write the re-binned hists
                outfile->cd(syst.c_str());
                for ( auto sys_hist : syst_map[syst] ) {
                    histoTrafo.rebinHisto( sys_hist, &bins);
                    sys_hist->Write();
                }
                outfile->cd();
            }

            // transform sig & bkg
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
