// Author      : Stefan Gadatsch
// Email       : gadatsch@nikhef.nl
// Date        : 2013-04-26
// Description : Compute uncertainty due to different groups of parameters specified in a XML

#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "Math/MinimizerOptions.h"
#include "TStopwatch.h"
#include "TXMLEngine.h"

#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooNLLVar.h"
#include "RooFitResult.h"

#include "../HepStatTool/HepMinSvc.h"

using namespace RooFit;
using namespace RooStats;

struct settings {
    string inFileName;
    string wsName;
    string modelConfigName;
    string dataName;
    string poiName;
    string xmlName;
    string technique;
    string catecory2eval;
    double precision;
    double corrCutoff;
    bool useMinos;
    string folder;
    string loglevel;
};

// ____________________________________________________________________________|__________
// Compute ranking of systematics specified in xml
void breakdown()
{
    std::string inFileName      = "./src/tutorial/histfact/output_ws/histfact_combined_BSM_model.root";
    std::string wsName          = "combined";
    std::string modelConfigName = "ModelConfig";
    TString dataName        = "obsData";
    std::string poiName         = "SigXsecOverSM";
    std::string xmlName         = "config/breakdown.xml";
    std::string technique       = "add";
    std::string catecory2eval   = "total";
    double precision            = 0.005;
    double corrCutoff           = 0.0;
    std::string folder          = "LQ3Up500GeV";
    std::string loglevel        = "DEBUG";

    // store all settings for passing to other functions
    settings* config = new settings();
    config->inFileName = inFileName;
    config->wsName = wsName;
    config->modelConfigName = modelConfigName;
    config->dataName = dataName;
    config->poiName = poiName;
    config->xmlName = xmlName;
    config->technique = technique;
    config->catecory2eval = catecory2eval;
    config->precision = precision;
    config->corrCutoff = corrCutoff;
    config->folder = folder;
    config->loglevel = loglevel;

    // some settings
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    if (config->loglevel == "DEBUG") {
        ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);
    } else {
        ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    }

    //vector<string> parsed = parseString(config->poiName, ",");

    TFile* file = new TFile(config->inFileName.c_str());

    RooWorkspace* ws = (RooWorkspace*)file->Get(config->wsName.c_str());
    if (!ws) {
        exit(1);
    }

    ModelConfig* mc = (ModelConfig*)ws->obj(config->modelConfigName.c_str());
    if (!mc) {
        exit(1);
    }

    RooDataSet* data = (RooDataSet*)ws->data(dataName);
    if (!data) {
        exit(1);
    }

    RooRealVar* poi_SigXsecOverSM = (RooRealVar*)mc->GetParametersOfInterest()->first();
    poi_SigXsecOverSM->setVal(1);
    poi_SigXsecOverSM->setRange(-10.,10.);
    poi_SigXsecOverSM->setConstant(1);

    RooArgSet* nuis = (RooArgSet*)mc->GetNuisanceParameters();
    if (!nuis) {
        exit(1);
    }
    TIterator* nitr = nuis->createIterator();
    RooRealVar* var;

    RooArgSet* globs = (RooArgSet*)mc->GetGlobalObservables();
    if (!globs) {
        exit(1);
    }

    ws->loadSnapshot("nominalNuis");
    poi_SigXsecOverSM->setRange(-10., 10.);
    poi_SigXsecOverSM->setConstant(0);
    poi_SigXsecOverSM->setVal(1.1); // Kick !

    RooNLLVar* nll = (RooNLLVar*)mc->GetPdf()->createNLL(*data, Constrain(*nuis), GlobalObservables(*globs), Offset(1), NumCPU(4, RooFit::Hybrid), Optimize(2));

    RooFitResult* fitresult = nullptr;
    HepMinSvc hepmin;
    hepmin.minuit(nll, "GlobalFit");

    RooArgSet nuisAndPOI(*mc->GetNuisanceParameters(), *mc->GetParametersOfInterest());
    ws->saveSnapshot("tmp_shot", nuisAndPOI);

    double nll_val_true = nll->getVal();
    double poi_hat = poi_SigXsecOverSM->getVal();

    double poi_up   = poi_SigXsecOverSM->getErrorHi();
    double poi_down = poi_SigXsecOverSM->getErrorLo();

    system(("mkdir -vp ./data/output/" +  string(config->folder) + "/root-files/breakdown_" + string(technique)).c_str());
    stringstream fileName;
    fileName << "./data/output/" << config->folder << "/root-files/breakdown_" << technique << "/" << config->catecory2eval << ".root";
    TFile fout(fileName.str().c_str(), "recreate");

    TH1D* h_out = new TH1D(config->catecory2eval.c_str(), config->catecory2eval.c_str(), 3, 0, 3);

    h_out->SetBinContent(1, poi_hat);
    h_out->SetBinContent(2, fabs(poi_up));
    h_out->SetBinContent(3, fabs(poi_down));

    h_out->GetXaxis()->SetBinLabel(1, poi_SigXsecOverSM->GetName());
    h_out->GetXaxis()->SetBinLabel(2, "poi_up");
    h_out->GetXaxis()->SetBinLabel(3, "poi_down");

    fout.Write();
    fout.Close();
}

