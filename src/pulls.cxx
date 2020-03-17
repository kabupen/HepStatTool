
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "Math/MinimizerOptions.h"
#include "TPRegexp.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"

#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooNLLVar.h"
#include "RooConstVar.h"
#include "RooFitResult.h"
#include "RooPoisson.h"
#include "RooGaussian.h"

#include "../HepStat/HepMinSvc.h"

using namespace RooFit;
using namespace RooStats;

enum NPType { GammaPoisson, GammaGaussian, GammaProtection, Others, Unknown};
std::pair<double,double> getPrefitErrorForGamma(const RooWorkspace *ws, const RooRealVar* par, NPType type);

void Draw(RooWorkspace* ws, RooAbsReal* nll, TString suffix = "")
{
    auto poi = ws->var("SigXsecOverSM");
    // Plot likelihood scan in parameter frac
    RooPlot* frame = poi->frame(Bins(10), Range(-10.0,10.0));
    nll->plotOn(frame,ShiftToZero());

    // Plot the profile likelihood in frac
//    RooAbsReal* pll_frac = nll->createProfile(*poi);
//    pll_frac->plotOn(frame, LineColor(kRed));

    TCanvas c("c","c");
    frame->Draw();
    c.SaveAs("SigXsecOverSM" + suffix + ".pdf");
}

void pulls()
{
//    const char* inFileName      = "1200GeV_WSMaker_workspace.root";
    std::string inFileName      = "./data/output_combined_VH_model_500GeV.root";
    const char* poiName         = "SigXsecOverSM";
    const char* wsName          = "combined";
    const char* modelConfigName = "ModelConfig";
    const char* dataName        = "obsData";
    std::string folder          = "LQ3Up500GeV";
    double precision = 0.005;

    /* Minimize interface */
    HepMinSvc hepmin;

    // some settings
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);

    // loading the workspace etc.
    std::cout << "Running over workspace: " << inFileName << std::endl;
    system(("mkdir -vp output/" + folder + "/root-files/pulls").c_str());

    TFile file(inFileName.c_str(), "OPEN");
    RooWorkspace* ws = (RooWorkspace*)file.Get(wsName);
    if (!ws) {
        std::cout << "Workspace: " << wsName << " doesn't exist!" << std::endl;
        exit(1);
    }

    ModelConfig* mc = (ModelConfig*)ws->obj(modelConfigName);
    if (!mc) {
        std::cout << "ModelConfig: " << modelConfigName << " doesn't exist!" << std::endl;
        exit(1);
    }

    TString datastring(dataName);
    int commaPos = datastring.Index(",");
    if(commaPos != TString::kNPOS) {
        ws->loadSnapshot(TString(datastring(commaPos+1, datastring.Length())));
        datastring = datastring(0, commaPos);
    }

    RooDataSet* data = (RooDataSet*)ws->data(datastring);
    if (!data) {
        std::cout << "Dataset: " << dataName << " doesn't exist!" << std::endl;
        exit(1);
    }

    std::vector<RooRealVar*> pois;
    TIterator* POI_itr = mc->GetParametersOfInterest()->createIterator();
    while ( RooRealVar*  poi = (RooRealVar*) POI_itr->Next() ) {
        std::cout << "Getting POI " << poi->GetName() << std::endl;
        poi->setVal(1);
        poi->setRange(-5.,5.);
        poi->setConstant(0);
        pois.push_back(poi);
    }

    /* Nuisance params */
    RooArgSet* nuis = (RooArgSet*)mc->GetNuisanceParameters();
    if (!nuis) {
        std::cout << "Nuisance parameter set doesn't exist!" << std::endl;
        exit(1);
    }
    std::vector<TString> nuisance_params;
    TIterator* itr = nuis->createIterator();
    std::cout << "Nuisance parameters" << std::endl;
    while ( RooRealVar* var = (RooRealVar*)itr->Next() ) {
        std::string varName = var->GetName();

        if ( varName.find("ATLAS_norm_All") != string::npos ) {
            std::cout << "Skipping " << varName << std::endl;
            continue;
        }

        // all remaining nuisance parameters
        std::cout << "\t" << varName << std::endl;
        nuisance_params.push_back(var->GetName());
    }
    itr->Reset();

    /* Global obs */
    std::cout << "Global observables" << std::endl;
    RooArgSet* globs = (RooArgSet*)mc->GetGlobalObservables();
    if (!globs) {
        std::cout <<  "GetGlobal observables don't exist!" << std::endl;
        exit(1);
    }
    TIterator* itr_globs = globs->createIterator();
    while ( RooRealVar* var = static_cast<RooRealVar*>(itr_globs->Next()) ){
        std::cout << "\t" << var->GetName() << std::endl;
    }

    /* create nll and do unconditional fit */ 
    // For unconditional fit, the POI should be floated
    // Unconditional maximum likelihood estimation (minimum log-likelihood estimation)
    for (unsigned int i = 0; i < pois.size(); i++) {
        pois[i]->setConstant(0);
        pois[i]->setRange(-50., 50.);
        pois[i]->setVal(1.1); // kick !
    }
    
    RooNLLVar* nll = (RooNLLVar*)mc->GetPdf()->createNLL(*data, Constrain(*nuis), GlobalObservables(*globs), Offset(1), Optimize(2));
    if ( hepmin.minuit(nll, "GlobalFit") ) {
        std::cout << "Finished to minimize the NLL." << std::endl; 
        std::cout << "\tNLL = " << nll->getVal() << std::endl;
        std::cout << "\tmu  = " << pois[0]->getVal() << std::endl;
    } else {
        return;
    }

    // it is good that first fit has decent estimate of the errors, so strategy 1 is good
    // for subsequent fits, we don't care, so go faster
    double nll_hat = nll->getVal();
    std::vector<double> pois_hat;
    for (unsigned int i = 0; i < pois.size(); i++) {
        std::cout << "Unconditional fit result : " << pois[i]->GetName() << " " << pois[i]->getVal() << std::endl;
        pois_hat.push_back(pois[i]->getVal());
    }

    // set all nuisance parameters floating 
    while ( RooRealVar* var = (RooRealVar*)itr->Next() ) {
        std::string varName = var->GetName();
        if ( varName.find("ATLAS_norm_All") != string::npos ) {
            var->setConstant(1);
        }
        var->setConstant(0);
    }
    // set POIs floating
    for (unsigned int i = 0; i < pois.size(); i++) {
        pois[i]->setConstant(0);
    }
    ws->saveSnapshot("tmp_snapshot", *mc->GetPdf()->getParameters(data));

    std::cout << "Made unconditional snapshot" << std::endl;
    std::cout << "----------------------------------BEGINNING FOR POI RANKING--------------------------------------------"<< std::endl;
    return;

    for ( unsigned int iPoi = 0; iPoi < pois.size(); iPoi++ ) {

        double poi_errup   = pois.at(iPoi)->getErrorHi();
        double poi_errdown = pois.at(iPoi)->getErrorLo();
        std::cout << __FILE__ << " " << __LINE__ << " " << pois.at(iPoi)->GetName() << " = " << pois_hat[iPoi] << " +" << fabs(poi_errup) << " /  -" << fabs(poi_errdown) << std::endl;

        /* fix theta at the MLE value +/- postfit uncertainty and minimize again to estimate the change in the POI */
        ws->loadSnapshot("tmp_snapshot");
        pois.at(iPoi)->setVal( pois_hat[iPoi] + poi_errup );
        pois.at(iPoi)->setConstant(1);
        hepmin.minuit(nll, "PoiUpFit"); 
        std::vector<double> pois_up;
        std::vector<double> pois_nom_up;
        for (unsigned int i = 0; i < pois.size(); i++ ) {
            if (iPoi==i)  pois_up.push_back(pois_hat[iPoi]);
            else          pois_up.push_back(pois[i]->getVal());
            pois_nom_up  .push_back(pois[i]->getVal());
        }

        ws->loadSnapshot("tmp_snapshot");
        pois.at(iPoi)->setVal( pois_hat[iPoi] + poi_errdown);
        pois.at(iPoi)->setConstant(1);
        hepmin.minuit(nll, "PoiDownFit"); 
        std::vector<double> pois_down;
        std::vector<double> pois_nom_down;
        for (unsigned int i = 0; i < pois.size(); i++) {
            if (iPoi==i)  pois_down.push_back(pois_hat[iPoi]);
            else          pois_down.push_back(pois[i]->getVal());
            pois_nom_down.push_back(pois[i]->getVal());
        }

        /* Output file in ROOT */
        TString fileName = "data/output/" + folder + "/root-files/pulls/" + pois.at(iPoi)->GetName() + ".root";
        TFile fout(fileName, "RECREATE");

        TH1D h_out(pois.at(iPoi)->GetName(), pois.at(iPoi)->GetName(), 3 + 5 * pois.size(), 0, 3 + 5 * pois.size());

        std::cout << pois_hat[iPoi] << " " << fabs(poi_errup) << " " << fabs(poi_errdown) << std::endl;

        h_out.SetBinContent(1, pois_hat[iPoi]);
        h_out.SetBinContent(2, fabs(poi_errup));
        h_out.SetBinContent(3, fabs(poi_errdown));

        h_out.GetXaxis()->SetBinLabel(1, "poi_hat");
        h_out.GetXaxis()->SetBinLabel(2, "poi_up");
        h_out.GetXaxis()->SetBinLabel(3, "poi_down");

        int bin = 4;
        for (unsigned int i = 0; i < pois.size(); i++) {
            h_out.SetBinContent(bin,   pois_hat[i]);
            h_out.SetBinContent(bin+1, pois_up[i]);
            h_out.SetBinContent(bin+2, pois_down[i]);
            h_out.SetBinContent(bin+3, pois_nom_up[i]);
            h_out.SetBinContent(bin+4, pois_nom_down[i]);

            h_out.GetXaxis()->SetBinLabel(bin,    pois[i]->GetName());
            h_out.GetXaxis()->SetBinLabel(bin+1, "poi_up");
            h_out.GetXaxis()->SetBinLabel(bin+2, "poi_down");
            h_out.GetXaxis()->SetBinLabel(bin+3, "poi_nom_up");
            h_out.GetXaxis()->SetBinLabel(bin+4, "poi_nom_down");

            bin += 5;
        }

        fout.Write();
        fout.Close();
    }
    std::cout << "---------------------------------------------------------END FOR POI RANKING-----------------------------------------------" << std::endl;
    return;

    std::cout << "Nuisance parameter loop : " << nuisance_params.size() << std::endl;
    for (int in = 0; in < nuisance_params.size(); in++) {

        ws->loadSnapshot("tmp_snapshot");

        RooRealVar* nuisance_parameter = (RooRealVar*)nuis->find(nuisance_params.at(in));
        std::cout<<" doing nuip name "<< nuisance_parameter->GetName() << std::endl;

        /* find all unconstrained NFs etc. */
        bool isNorm = false;
        if ( ( static_cast<std::string>(nuisance_parameter->GetName()).find("ATLAS_norm") != std::string::npos ) || 
             ( static_cast<std::string>(nuisance_parameter->GetName()).find("gamma")      != std::string::npos ) || 
             ( static_cast<std::string>(nuisance_parameter->GetName()).find("scale_") && static_cast<std::string>(nuisance_parameter->GetName()).find("QCDscale_") == std::string::npos )) isNorm = true;

        double nuip_hat = nuisance_parameter->getVal();
        nuisance_parameter->setConstant(0);

        ws->saveSnapshot("tmp_snapshot2", *mc->GetPdf()->getParameters(data));

        std::cout << "Computing error for var " << nuisance_parameter->GetName() << " at " << nuisance_parameter->getVal() << std::endl;
        std::cout << __FILE__ << " " << __LINE__ << " " << nuip_hat << std::endl;

        /* should be careful for Poisson PDFs that are not defined for negative NP */
        NPType np_type = NPType::Unknown;
        std::string par_name = nuisance_parameter->GetName();
        if( par_name.find("gamma") != std::string::npos ){
            RooAbsPdf* constraint_term = (RooAbsPdf*)ws->pdf((par_name+"_constraint").c_str());
            if( constraint_term ){
                if( typeid(*constraint_term)==typeid(RooPoisson) ){
                    RooConstVar *var_tau = (RooConstVar*)ws->obj((par_name+"_tau").c_str());
                    if( var_tau->getVal() >0.9 ){
                        np_type = NPType::GammaPoisson;
                    } else {
                        np_type = NPType::GammaProtection;
                    }
                }
                else if(typeid(*constraint_term)==typeid(RooGaussian)){
                    np_type = NPType::GammaGaussian;
                }
            }
        } else {
            np_type = NPType::Others;
        }
        
        double nuip_high = nuip_hat+fabs(nuisance_parameter->getErrorHi());
        double nuip_low  = nuip_hat-fabs(nuisance_parameter->getErrorLo());
        std::cout << __FILE__ << " " << __LINE__ << " " << nuip_high << std::endl;
        std::cout << __FILE__ << " " << __LINE__ << " " << nuip_low  << std::endl;
        if     ( np_type == NPType::GammaPoisson && nuip_low<=0 ){ nuip_low  = getPrefitErrorForGamma(ws, nuisance_parameter, np_type).first;  }
        else if( np_type == NPType::GammaProtection )            { nuip_high = getPrefitErrorForGamma(ws, nuisance_parameter, np_type).second; } // This is better estimation and faster convergence 

        double nuip_errup;
        double nuip_errdown;
        if( np_type != NPType::GammaProtection ){
            ws->loadSnapshot("tmp_snapshot2");
            // nuip_errup = hepmin.findSigma(nll, nll_hat, nuisance_parameter, nuip_hat, +1);
            nuip_errup = nuisance_parameter->getErrorHi();
            ws->loadSnapshot("tmp_snapshot2");
            //nuip_errdown = hepmin.findSigma(nll, nll_hat, nuisance_parameter, nuip_hat, -1);
            nuip_errdown = nuisance_parameter->getErrorLo();
        } else {
            ws->loadSnapshot("tmp_snapshot2");
            // nuip_errup = hepmin.findSigma(nll, nll_hat, nuisance_parameter, nuip_hat, +1);
             nuip_errup = nuisance_parameter->getErrorHi();
            std::cout << "This is a gamma parameter for protection. No need to get -1 sigma." << std::endl;
            nuip_errdown = 0.;
        }
        std::cout << __FILE__ << " " << __LINE__ << " " << nuip_errup << std::endl;
        std::cout << __FILE__ << " " << __LINE__ << " " << nuip_errdown << std::endl;

        std::cout << nuisance_parameter->GetName() << " = " << nuip_hat << " +" << fabs(nuip_errup) << " /  -" << fabs(nuip_errdown) << std::endl;

        // fix theta at the MLE value +/- postfit uncertainty and minimize again to estimate the change in the POI
        ws->loadSnapshot("tmp_snapshot2");
        nuisance_parameter->setVal(nuip_hat+fabs(nuip_errup));
        nuisance_parameter->setConstant(1);
        hepmin.minuit(nll); // RooMinuit(*nll).migrad(); // Fits::minimize(nll);
        std::vector<double> pois_up;
        for (unsigned int i = 0; i < pois.size(); i++) {
            pois_up.push_back(pois[i]->getVal());
        }

        ws->loadSnapshot("tmp_snapshot2");
        nuisance_parameter->setVal(nuip_hat-fabs(nuip_errdown));
        nuisance_parameter->setConstant(1);
        hepmin.minuit(nll); // RooMinuit(*nll).migrad(); // Fits::minimize(nll);
        std::vector<double> pois_down;
        for (unsigned int i = 0; i < pois.size(); i++) {
            pois_down.push_back(pois[i]->getVal());
        }

        // fix theta at the MLE value +/- postfit uncertainty and minimize again to estimate the change in the POI
        ws->loadSnapshot("tmp_snapshot2");
        std::vector<double> pois_nom_up;
        if(!isNorm) {
            nuisance_parameter->setVal(nuip_hat+1.0);
            nuisance_parameter->setConstant(1);
            hepmin.minuit(nll); // RooMinuit(*nll).migrad(); // Fits::minimize(nll);
        }
        for (unsigned int i = 0; i < pois.size(); i++) {
            pois_nom_up.push_back(pois[i]->getVal());
        }

        ws->loadSnapshot("tmp_snapshot2");
        if(!isNorm) {
            nuisance_parameter->setVal(nuip_hat-1.0);
            nuisance_parameter->setConstant(1);
            hepmin.minuit(nll); // RooMinuit(*nll).migrad(); // Fits::minimize(nll);
        }
        std::vector<double> pois_nom_down;
        for (unsigned int i = 0; i < pois.size(); i++) {
            pois_nom_down.push_back(pois[i]->getVal());
        }

        for (unsigned int i = 0; i < pois.size(); i++) {
            std::cout << "Variation of " << pois[i]->GetName() << " = " << pois_up[i] << " (" << pois_nom_up[i] << ") / " << pois_down[i] << " (" << pois_nom_down[i] << ")" << std::endl;
        }

        // store result in root file
        TString fileName =  "output/" + folder + "/root-files/pulls/" + nuisance_parameter->GetName() + ".root";
        TFile fout(fileName, "RECREATE");

        TH1D h_out( nuisance_parameter->GetName(), nuisance_parameter->GetName(), 3 + 5 * pois.size(), 0, 3 + 5 * pois.size());

        h_out.SetBinContent(1, nuip_hat);
        h_out.SetBinContent(2, fabs(nuip_errup));
        h_out.SetBinContent(3, fabs(nuip_errdown));

        h_out.GetXaxis()->SetBinLabel(1, "nuip_hat");
        h_out.GetXaxis()->SetBinLabel(2, "nuip_up");
        h_out.GetXaxis()->SetBinLabel(3, "nuip_down");

        int bin = 4;
        for (unsigned int i = 0; i < pois.size(); i++) {
            h_out.SetBinContent(bin, pois_hat[i]);
            h_out.SetBinContent(bin+1, pois_up[i]);
            h_out.SetBinContent(bin+2, pois_down[i]);
            h_out.SetBinContent(bin+3, pois_nom_up[i]);
            h_out.SetBinContent(bin+4, pois_nom_down[i]);

            h_out.GetXaxis()->SetBinLabel(bin, pois[i]->GetName());
            h_out.GetXaxis()->SetBinLabel(bin+1, "poi_up");
            h_out.GetXaxis()->SetBinLabel(bin+2, "poi_down");
            h_out.GetXaxis()->SetBinLabel(bin+3, "poi_nom_up");
            h_out.GetXaxis()->SetBinLabel(bin+4, "poi_nom_down");

            bin += 5;
        }

        fout.Write();
        fout.Close();
    }
}

std::pair<double,double> getPrefitErrorForGamma(const RooWorkspace *ws, const RooRealVar* par, NPType type)
{
    double xlo = -999.;
    double xhi = -999.;
    std::string par_name = par->GetName();
    if(par_name.find("gamma") != std::string::npos){
        RooAbsPdf* constraint_term = (RooAbsPdf*)ws->pdf((par_name+"_constraint").c_str());
        if(constraint_term){
            if (type==NPType::GammaPoisson || type==NPType::GammaProtection){
                if(typeid(*constraint_term)==typeid(RooPoisson)){
                    RooConstVar *var_tau = (RooConstVar*)ws->obj((par_name+"_tau").c_str());
                    double ylim = 0.5;
                    double tau = var_tau->getVal();
                    double xlim_pois = TMath::Exp(-1-ylim/tau); // <= -ln{(gamma*tau)^tau*exp(-tau)} = -ln{tau^tau}*ylim
                    double xlim_gaus = 1-sqrt(2*ylim/tau);     // <= -ln{exp(-(tau/2)*(gamma-1)^2)} = ylim
                    if(type==NPType::GammaPoisson){
                        if(xlim_gaus>0.){
                            xlo = xlim_gaus;
                        } else {
                            xlo = xlim_pois;
                        }
                        xhi = 2-xlim_gaus; //=1+(1-xlim_gaus)
                        std::cout << "getPrefitErrorForGamma    parameter:" << par_name
                            << ". This is RooPoisson with tau:" << tau
                            << " xlim_pois:" << xlim_pois
                            << " xlim_gaus:" << xlim_gaus
                            << ".";
                    } else {
                        xlo = xlim_pois;
                        //xhi = -(1./tau)*TMath::Log(1-0.6827); // 68.27% qunatile for exp(-tau*gamma)
                        xhi = 1./2./tau;  // <= -ln(exp(-tau*gamma)) = 0.5
                        std::cout << "getPrefitErrorForGamma    parameter:" << par_name
                            << ". This is RooPoisson for protecion with tau:" << tau
                            << ".";
                    }
                }
            }
            else if(type==NPType::GammaGaussian){
                if(typeid(*constraint_term)==typeid(RooGaussian)){
                    type = NPType::GammaGaussian;
                    RooConstVar *var_sigma = (RooConstVar*)ws->obj((par_name+"_sigma").c_str());
                    double sigma = var_sigma->getVal();
                    xlo = 1.-sigma;
                    xhi = 1.+sigma;
                    std::cout << "getPrefitErrorForGamma    parameter:" << par_name
                        << ". This is RooGaussian with sigma:" << sigma
                        << ".";
                }
            }
            std::cout << " xlo:" << xlo << " xhi:" << xhi << std::endl;
        }
    }
    return std::make_pair(xlo, xhi);
}
