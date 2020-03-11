
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

using namespace RooFit;
using namespace RooStats;

class PullCalc
{
    public:
        PullCalc();
        ~PullCalc();

        void MinimizeNLL();
        void Draw(TString suffix = "");
    private:
        TFile*        m_file;
        RooWorkspace* m_workspace;
        RooStats::ModelConfig*  m_modelConfig;
        RooDataSet*   m_data;
        RooArgSet*    m_nuis;
        RooArgSet*    m_globs;
        std::vector<RooRealVar*> m_pois;
        std::vector<std::string> m_nuis_name;
        RooNLLVar*  m_nll;
};

PullCalc::PullCalc()
{
    const char* inFileName      = "hist2workspace/output_combined_VH_model.root";
    const char* poiName         = "SigXsecOverSM";
    const char* wsName          = "combined";
    const char* modelConfigName = "ModelConfig";
    const char* dataName        = "obsData";
    double precision = 0.005;

    // some settings
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);

    // loading the workspace etc.
    std::cout << "Running over workspace: " << inFileName << std::endl;
    std::system("mkdir -vp output/test/root-files/pulls");

    m_file = new TFile(inFileName, "OPEN");
    m_workspace = (RooWorkspace*)m_file->Get(wsName);
    if (!m_workspace) {
        std::cout << "Workspace: " << wsName << " doesn't exist!" << std::endl;
        exit(1);
    } else {
        std::cout << "Workspace = " << m_workspace->GetName() << "\n" << std::endl;
    }

    m_modelConfig = (RooStats::ModelConfig*)m_workspace->obj(modelConfigName);
    if (!m_modelConfig) {
        std::cout << "ModelConfig: " << modelConfigName << " doesn't exist!" << std::endl;
        exit(1);
    } else {
        std::cout << "ModelConfig = " << m_modelConfig->GetName() << "\n" << std::endl;
    }

    m_data = (RooDataSet*)m_workspace->data(dataName);
    if (!m_data) {
        std::cout << "Dataset: " << dataName << " doesn't exist!" << std::endl;
        exit(1);
    } else {
        std::cout << "Data = " << m_data->GetName() << "\n" << std::endl;
    }

    TIterator* itr_poi = m_modelConfig->GetParametersOfInterest()->createIterator();
    while ( RooRealVar*  poi = (RooRealVar*) itr_poi->Next() ) {
        std::cout << "Getting POI " << poi->GetName() << std::endl;
        poi->setVal(1);
        poi->setRange(-5.,5.);
        poi->setConstant(0);
        m_pois.push_back(poi);
    }

    /* Nuisance params */
    m_nuis = (RooArgSet*)m_modelConfig->GetNuisanceParameters();
    if (!m_nuis) {
        std::cout << "Nuisance parameter set doesn't exist!" << std::endl;
        exit(1);
    }
    TIterator* itr = m_nuis->createIterator();
    std::cout << "Nuisance parameters" << std::endl;
    while ( RooRealVar* var = (RooRealVar*)itr->Next() ) {
        std::string varName = var->GetName();

        if ( varName.find("ATLAS_norm_All") != string::npos ) {
            std::cout << "Skipping " << varName << std::endl;
            continue;
        }

        // all remaining nuisance parameters
        std::cout << "\t" << varName << std::endl;
        m_nuis_name.push_back(std::string(var->GetName()));
    }

    /* Global obs */
    std::cout << "Global observables" << std::endl;
    m_globs = (RooArgSet*)m_modelConfig->GetGlobalObservables();
    if (!m_globs) {
        std::cout <<  "GetGlobal observables don't exist!" << std::endl;
        exit(1);
    }
    TIterator* itr_globs = m_globs->createIterator();
    while ( RooRealVar* var = static_cast<RooRealVar*>(itr_globs->Next()) ){
        std::cout << "\t" << var->GetName() << std::endl;
    }
}

PullCalc::~PullCalc()
{
    delete m_file;
    delete m_workspace;
    delete m_modelConfig;
    delete m_data;
    delete m_nuis;
    delete m_globs;
    //        delete std::vector<RooRealVar*> m_pois;
    //        delete std::vector<std::string> m_nuis_name;
    delete m_nll;
}

void PullCalc::MinimizeNLL()
{
    /* For unconditional fit, the POI should be floated 
     * Unconditional maximum likelihood estimation (minimum log-likelihood estimation)
     */
    std::cout << "Create nll" << std::endl;
    for (unsigned int i = 0; i < m_pois.size(); i++) {
        m_pois[i]->setConstant(0);
        m_pois[i]->setRange(-50., 50.);
        m_pois[i]->setVal(1.1); // kick !
    }

    m_nll = (RooNLLVar*)m_modelConfig->GetPdf()->createNLL(*m_data, Constrain(*m_nuis), GlobalObservables(*m_globs));

    std::cout << "Minimize the created NLL" << std::endl;
    RooMinuit(*m_nll).migrad(); 
    std::cout << "Finished to minimize the NLL." << std::endl;
    std::cout << m_nll->getVal() << " " << m_pois[0]->getVal() << std::endl;
}

void PullCalc::Draw(TString suffix)
{
    auto poi = m_workspace->var("SigXsecOverSM");

    RooPlot* frame = poi->frame(Bins(10), Range(-50.0,50.0));
    m_nll->plotOn(frame, ShiftToZero());

    TCanvas c;
    frame->Draw();
    c.SaveAs("SigXsecOverSM" + suffix + ".pdf");
}


