// Author      : Stefan Gadatsch
// Email       : gadatsch@nikhef.nl
// Date        : 2013-04-24
// Description : Compute pulls and impact on the POI

#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "Math/MinimizerOptions.h"
#include "TStopwatch.h"
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

#include "findSigma.cxx"

using namespace std;
using namespace RooFit;
using namespace RooStats;

enum NPType { GammaPoisson, GammaGaussian, GammaProtection, Others, Unknown};
NPType getNPType(const RooWorkspace *ws, const RooRealVar *par);
std::pair<double,double> getPrefitErrorForGamma(const RooWorkspace *ws, const RooRealVar* par, NPType type);

// compute pulls of the nuisance parameters and store them in text files. 
// norm and syst parameters will be split among different files
// ROOT.runPulls("output/"+ws+"/workspaces/combined/"+mass+".root", poi, "combined",modelConfig, dataName, ws)
void pulls(
        const char* inFileName = "LQ3LH_v10.output_LQ3_13TeV_output_Systs_lephad_BasicKinematics_FullRun2_TauPT_300/workspaces/combined/300.root",
        const char* poiName = "SigXsecOverSM",
        const char* wsName = "combined",
        const char* modelConfigName = "ModelConfig",
        const char* dataName = "obsData",
        const char* folder = "test",

        const char* variable = NULL
        )
{

    double precision = 0.005;
    int nJobs = 1;
    int iJob = 0;

    // some settings
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);

    // loading the workspace etc.
    std::cout << "Running over workspace: " << inFileName << std::endl;
    system(("mkdir -vp output/" + string(folder) + "/root-files/pulls").c_str());

    TFile* file = new TFile(inFileName);
    RooWorkspace* ws = (RooWorkspace*)file->Get(wsName);
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

    RooArgSet* nuis = (RooArgSet*)mc->GetNuisanceParameters();
    if (!nuis) {
        std::cout << "Nuisance parameter set doesn't exist!" << std::endl;
        exit(1);
    }

    RooArgSet* globs = (RooArgSet*)mc->GetGlobalObservables();
    if (!globs) {
        std::cout <<  "GetGlobal observables don't exist!" << std::endl;
        exit(1);
    }

    /* collect nuisance parameters */
    std::vector<string> vec_nuis;
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
        vec_nuis.push_back(string(var->GetName()));
    }

    itr->Reset();

    /* create nll and do unconditional fit */ 
    // For unconditional fit, the POI should be floated
    // Unconditional maximum likelihood estimation (minimum log-likelihood estimation)
    for (unsigned int i = 0; i < pois.size(); i++) {
        pois[i]->setConstant(0);
        pois[i]->setRange(-5., 5.);
        pois[i]->setVal(1.1); // kick !
    }
    RooNLLVar* nll = (RooNLLVar*)mc->GetPdf()->createNLL(*data, Constrain(*nuis), GlobalObservables(*globs), Offset(1), Optimize(2));
    RooMinuit(*nll).migrad(); // Fits::minimize(nll);
    std::cout << "Finished to minimize the NLL." << std::endl;

    // it is good that first fit has decent estimate of the errors, so strategy 1 is good
    // for subsequent fits, we don't care, so go faster
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
    double nll_hat = nll->getVal();
    std::cout << "NLL getVal() : " << nll->getVal() << std::endl;
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

    for ( unsigned int in = 0; in < pois.size();in++ ) {

        if ( pois.size()<(unsigned int)nJobs && in != (unsigned int)iJob ) continue; 

        double poi_errup;
        double poi_errdown;

        ws->loadSnapshot("tmp_snapshot");
        //poi_errup   = findSigma(nll, nll_hat, pois[in], pois_hat[in]+fabs(pois[in]->getErrorHi()), pois_hat[in], +1, precision); 
        poi_errup   = findSigma(nll, nll_hat, pois[in], pois_hat[in], +1); 
        ws->loadSnapshot("tmp_snapshot");
        //poi_errdown = findSigma(nll, nll_hat, pois[in], pois_hat[in]-fabs(pois[in]->getErrorLo()), pois_hat[in], -1, precision);
        poi_errdown = findSigma(nll, nll_hat, pois[in], pois_hat[in], -1);
        std::cout << __FILE__ << " " << __LINE__ << " " << pois[in]->GetName() << " = " << pois_hat[in] << " +" << fabs(poi_errup) << " /  -" << fabs(poi_errdown) << std::endl;

        // fix theta at the MLE value +/- postfit uncertainty and minimize again to estimate the change in the POI
        ws->loadSnapshot("tmp_snapshot");
        pois[in]->setVal(pois_hat[in]+fabs(poi_errup));
        pois[in]->setConstant(1);
        Fits::minimize(nll);
        std::vector<double> pois_up;
        for (unsigned int i = 0; i < (pois.size()); i++) {
            if (in==i)  pois_up.push_back(pois_hat[in]);
            else pois_up.push_back(pois[i]->getVal());
        }

        ws->loadSnapshot("tmp_snapshot");
        pois[in]->setVal(pois_hat[in]-fabs(poi_errdown));
        pois[in]->setConstant(1);
        Fits::minimize(nll);
        vector<double> pois_down;
        for (unsigned int i = 0; i < pois.size(); i++) {
            if (in==i)  pois_down.push_back(pois_hat[in]);
            else pois_down.push_back(pois[i]->getVal());
        }

        // fix theta at the MLE value +/- postfit uncertainty and minimize again to estimate the change in the POI
        ws->loadSnapshot("tmp_snapshot");
        std::vector<double> pois_nom_up, pois_nom_down;
        for (unsigned int i = 0; i < pois.size(); i++) {
            pois_nom_up.push_back(pois[i]->getVal());
            pois_nom_down.push_back(pois[i]->getVal());
        }

        // store result in root file
        std::stringstream fileName;
        fileName << "output/" << folder << "/root-files/pulls/" << pois[in]->GetName() << ".root";
        TFile fout(fileName.str().c_str(), "recreate");

        TH1D* h_out = new TH1D((pois[in]->GetName()), (pois[in]->GetName()), 3 + 5 * pois.size(), 0, 3 + 5 * pois.size());

        h_out->SetBinContent(1, pois_hat[in]);
        h_out->SetBinContent(2, fabs(poi_errup));
        h_out->SetBinContent(3, fabs(poi_errdown));

        h_out->GetXaxis()->SetBinLabel(1, "poi_hat");
        h_out->GetXaxis()->SetBinLabel(2, "poi_up");
        h_out->GetXaxis()->SetBinLabel(3, "poi_down");

        int bin = 4;
        for (unsigned int i = 0; i < pois.size(); i++) {
            h_out->SetBinContent(bin, pois_hat[i]);
            h_out->SetBinContent(bin+1, pois_up[i]);
            h_out->SetBinContent(bin+2, pois_down[i]);
            h_out->SetBinContent(bin+3, pois_nom_up[i]);
            h_out->SetBinContent(bin+4, pois_nom_down[i]);

            h_out->GetXaxis()->SetBinLabel(bin, pois[i]->GetName());
            h_out->GetXaxis()->SetBinLabel(bin+1, "poi_up");
            h_out->GetXaxis()->SetBinLabel(bin+2, "poi_down");
            h_out->GetXaxis()->SetBinLabel(bin+3, "poi_nom_up");
            h_out->GetXaxis()->SetBinLabel(bin+4, "poi_nom_down");

            bin += 5;

        }

        fout.Write();
        fout.Close();
    }
    //---------------------------------------------------------END FOR POI RANKING-----------------------------------------------//

    std::cout << "Nuisance parameter loop : " << vec_nuis.size() << std::endl;
    //for (int in = 0; in < vec_nuis.size(); in++) {
    for (int in = 0; in < 3; in++) {

        ws->loadSnapshot("tmp_snapshot");

        RooRealVar* nuip = (RooRealVar*)nuis->find(vec_nuis[in].c_str());
        std::cout<<" doing nuip name "<< nuip->GetName() << std::endl;
        string nuipName(nuip->GetName());

        if (variable != NULL && nuipName != string(variable)) continue;

        // find all unconstrained NFs etc.
        bool isNorm = 0;
        if (nuipName.find("ATLAS_norm") != string::npos) isNorm = 1;
        if (nuipName.find("gamma") != string::npos) isNorm = 1;
        if (nuipName.find("scale_") != string::npos && nuipName.find("QCDscale_") == string::npos) isNorm = true;

        double nuip_hat = nuip->getVal();
        nuip->setConstant(0);

        ws->saveSnapshot("tmp_snapshot2", *mc->GetPdf()->getParameters(data));

        std::cout << "Computing error for var " << nuip->GetName() << " at " << nuip->getVal() << std::endl;

        double nuip_errup;
        double nuip_errdown;

        // should be careful for Poisson PDFs that are not defined for negative NP
        double nuip_high = nuip_hat+fabs(nuip->getErrorHi());
        double nuip_low  = nuip_hat-fabs(nuip->getErrorLo());
        NPType np_type = getNPType(ws, nuip);

        if( np_type==NPType::GammaPoisson && nuip_low<=0 ){
            nuip_low = getPrefitErrorForGamma(ws, nuip, np_type).first;
        } else if(np_type==NPType::GammaProtection){
            // This is better estimation and faster convergence
            nuip_high =  getPrefitErrorForGamma(ws, nuip, np_type).second;
        }

        if(np_type!=NPType::GammaProtection){
            ws->loadSnapshot("tmp_snapshot2");
            //nuip_errup = findSigma(nll, nll_hat, nuip, nuip_high, nuip_hat, +1, precision);
            nuip_errup = findSigma(nll, nll_hat, nuip, nuip_hat, +1);
            ws->loadSnapshot("tmp_snapshot2");
//            nuip_errdown = findSigma(nll, nll_hat, nuip, nuip_low, nuip_hat, -1, precision);
            nuip_errdown = findSigma(nll, nll_hat, nuip, nuip_hat, -1);
        } else {
            ws->loadSnapshot("tmp_snapshot2");
//            nuip_errup = findSigma(nll, nll_hat, nuip, nuip_high, nuip_hat, +1, precision);
            nuip_errup = findSigma(nll, nll_hat, nuip, nuip_hat, +1);
            std::cout << "This is a gamma parameter for protection. No need to get -1 sigma." << std::endl;
            nuip_errdown = 0.;
        }

        std::cout << nuip->GetName() << " = " << nuip_hat << " +" << fabs(nuip_errup) << " /  -" << fabs(nuip_errdown) << std::endl;

        // fix theta at the MLE value +/- postfit uncertainty and minimize again to estimate the change in the POI
        ws->loadSnapshot("tmp_snapshot2");
        nuip->setVal(nuip_hat+fabs(nuip_errup));
        nuip->setConstant(1);
        Fits::minimize(nll);
        vector<double> pois_up;
        for (unsigned int i = 0; i < pois.size(); i++) {
            pois_up.push_back(pois[i]->getVal());
        }

        ws->loadSnapshot("tmp_snapshot2");
        nuip->setVal(nuip_hat-fabs(nuip_errdown));
        nuip->setConstant(1);
        Fits::minimize(nll);
        vector<double> pois_down;
        for (unsigned int i = 0; i < pois.size(); i++) {
            pois_down.push_back(pois[i]->getVal());
        }

        // fix theta at the MLE value +/- postfit uncertainty and minimize again to estimate the change in the POI
        ws->loadSnapshot("tmp_snapshot2");
        vector<double> pois_nom_up;
        if(!isNorm) {
            nuip->setVal(nuip_hat+1.0);
            nuip->setConstant(1);
            Fits::minimize(nll);
        }
        for (unsigned int i = 0; i < pois.size(); i++) {
            pois_nom_up.push_back(pois[i]->getVal());
        }

        ws->loadSnapshot("tmp_snapshot2");
        if(!isNorm) {
            nuip->setVal(nuip_hat-1.0);
            nuip->setConstant(1);
            Fits::minimize(nll);
        }
        vector<double> pois_nom_down;
        for (unsigned int i = 0; i < pois.size(); i++) {
            pois_nom_down.push_back(pois[i]->getVal());
        }

        for (unsigned int i = 0; i < pois.size(); i++) {
            std::cout << "Variation of " << pois[i]->GetName() << " = " << pois_up[i] << " (" << pois_nom_up[i] << ") / " << pois_down[i] << " (" << pois_nom_down[i] << ")" << std::endl;
        }

        // store result in root file
        stringstream fileName;
        fileName << "output/" << folder << "/root-files/pulls/" << nuipName << ".root";
        TFile fout(fileName.str().c_str(), "recreate");

        TH1D* h_out = new TH1D(nuipName.c_str(), nuipName.c_str(), 3 + 5 * pois.size(), 0, 3 + 5 * pois.size());

        h_out->SetBinContent(1, nuip_hat);
        h_out->SetBinContent(2, fabs(nuip_errup));
        h_out->SetBinContent(3, fabs(nuip_errdown));

        h_out->GetXaxis()->SetBinLabel(1, "nuip_hat");
        h_out->GetXaxis()->SetBinLabel(2, "nuip_up");
        h_out->GetXaxis()->SetBinLabel(3, "nuip_down");

        int bin = 4;
        for (unsigned int i = 0; i < pois.size(); i++) {
            h_out->SetBinContent(bin, pois_hat[i]);
            h_out->SetBinContent(bin+1, pois_up[i]);
            h_out->SetBinContent(bin+2, pois_down[i]);
            h_out->SetBinContent(bin+3, pois_nom_up[i]);
            h_out->SetBinContent(bin+4, pois_nom_down[i]);

            h_out->GetXaxis()->SetBinLabel(bin, pois[i]->GetName());
            h_out->GetXaxis()->SetBinLabel(bin+1, "poi_up");
            h_out->GetXaxis()->SetBinLabel(bin+2, "poi_down");
            h_out->GetXaxis()->SetBinLabel(bin+3, "poi_nom_up");
            h_out->GetXaxis()->SetBinLabel(bin+4, "poi_nom_down");

            bin += 5;
        }

        fout.Write();
        fout.Close();
    }
}


NPType getNPType(const RooWorkspace *ws, const RooRealVar *par)
{
    NPType type = NPType::Unknown;
    std::string par_name = par->GetName();
    if( par_name.find("gamma") != std::string::npos ){
        RooAbsPdf* constraint_term = (RooAbsPdf*)ws->pdf((par_name+"_constraint").c_str());
        if( constraint_term ){
            if(typeid(*constraint_term)==typeid(RooPoisson)){
                RooConstVar *var_tau = (RooConstVar*)ws->obj((par_name+"_tau").c_str());
                double tau = var_tau->getVal();
                if(tau>0.9){
                    type = NPType::GammaPoisson;
                } else {
                    type = NPType::GammaProtection;
                }
            }
            else if(typeid(*constraint_term)==typeid(RooGaussian)){
                type = NPType::GammaGaussian;
            }
        }
    } else {
        type = NPType::Others;
    }
    return type;
}


std::pair<double,double> getPrefitErrorForGamma(const RooWorkspace *ws, const RooRealVar* par, NPType type){

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
