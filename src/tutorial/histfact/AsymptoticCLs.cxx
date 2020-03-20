
using namespace RooFit;
using namespace RooStats;

void AsymptoticCLs()
{
    TString filename = "output_ws/histfact_combined_BSM_model.root";
    bool fileExist = !gSystem->AccessPathName(filename); // note opposite return code
    // if file does not exists generate with histfactory
    if (!fileExist) {
        std::cout << "@@@@ File not found : " << filename << std::endl;
        return;
    }
    // Try to open the file
    TFile *file = TFile::Open(filename);

    RooWorkspace *w = (RooWorkspace*)file->Get("combined");
    w->Print();

    RooAbsData *data = w->data("obsData");

    // (mUseVectorStore = true
    RooAbsData::setDefaultStorageType(RooAbsData::Vector);
    data->convertToVectorStore();

    // get models from WS
    ModelConfig *sbModel = (ModelConfig *)w->obj("ModelConfig");
    sbModel->SetName("S+B_Model");
    sbModel->SetSnapshot(*sbModel->GetParametersOfInterest()); 
    
    // Create B only model by using S+B model with poi=0
    // Set the POI as 0
    ModelConfig* bModel = (ModelConfig *)sbModel->Clone();
    bModel->SetName("B_Model");
    RooRealVar *var_poi = dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());
    if (!var_poi) return 0;
    double oldval = var_poi->getVal();
    var_poi->setVal(0);
    bModel->SetSnapshot(RooArgSet(*var_poi));

    // check model  has global observables when there are nuisance pdf
    // for the hybrid case the globals are not needed
    bool hasNuisParam = (sbModel->GetNuisanceParameters() && sbModel->GetNuisanceParameters()->getSize() > 0);
    bool hasGlobalObs = (sbModel->GetGlobalObservables()  && sbModel->GetGlobalObservables()->getSize()  > 0);
    if (hasNuisParam && !hasGlobalObs) {
        // try to see if model has nuisance parameters first
        RooAbsPdf *constrPdf = RooStats::MakeNuisancePdf(*sbModel, "nuisanceConstraintPdf_sbmodel");
        if (constrPdf) {
            std::cout << "Model" << sbModel->GetName() << " has nuisance parameters but no global observables associated" << std::endl;
            std::cout << "The effect of the nuisance parameters will not be treated correctly " << std::endl;
        }
    }

    // save all initial parameters of the model including the global observables
    RooArgSet initialParameters;
    RooArgSet*allParams = sbModel->GetPdf()->getParameters(*data);
    allParams->snapshot(initialParameters);
    delete allParams;

    // run first a data fit
    std::string minimizer_type = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
    std::cout << "Using as minimizer for computing the test statistic"<< minimizer_type << std::endl;

    // create the HypoTest calculator class (Asimov data generated with nominal values)
    AsymptoticCalculator::SetPrintLevel(0);
    HypoTestCalculatorGeneric *hc = new AsymptoticCalculator(*data, *bModel, *sbModel, true); 

    ((AsymptoticCalculator *)hc)->SetOneSided(true);

    // Get the result
    RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);

    HypoTestInverter hypo(*hc);
    hypo.SetConfidenceLevel(0.95);
    hypo.UseCLs(true);
    hypo.SetVerbose(true);
    const double poimin = 0; 
    const double poimax = 5;
    const int npoints = 20; 
    hypo.SetFixedScan(npoints, poimin, poimax);

    HypoTestInverterResult* hypo_inverter_result = hypo.GetInterval();
    if (!hypo_inverter_result) { std::cerr << "Error running the HypoTestInverter - Exit " << std::endl; return; }

    // analyze result produced by the inverter, optionally save it in a file

    double lowerLimit = hypo_inverter_result->LowerLimit();
    double llError    = hypo_inverter_result->LowerLimitEstimatedError();
    double upperLimit = hypo_inverter_result->UpperLimit();
    double ulError    = hypo_inverter_result->UpperLimitEstimatedError();

    if (lowerLimit < upperLimit * (1. - 1.E-4) && lowerLimit != 0)
        std::cout << "The computed lower limit is: " << lowerLimit << " +/- " << llError << std::endl;
    std::cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError << std::endl;

    // compute expected limit
    std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
    std::cout << " expected limit (median) " << hypo_inverter_result->GetExpectedUpperLimit(0) << std::endl;
    std::cout << " expected limit (-1 sig) " << hypo_inverter_result->GetExpectedUpperLimit(-1) << std::endl;
    std::cout << " expected limit (+1 sig) " << hypo_inverter_result->GetExpectedUpperLimit(1) << std::endl;
    std::cout << " expected limit (-2 sig) " << hypo_inverter_result->GetExpectedUpperLimit(-2) << std::endl;
    std::cout << " expected limit (+2 sig) " << hypo_inverter_result->GetExpectedUpperLimit(2) << std::endl;

    // write result in a file
    TString outputname = TString::Format("Asym_CLs_%s.root", (npoints < 0) ? "auto" : "grid");
    TFile *fileOut = new TFile( outputname, "RECREATE");
    hypo_inverter_result->Write();
    fileOut->Close();

    // plot the result ( p values vs scan points)
    TString plotTitle = TString::Format("Asymptotic CL Scan for workspace %s", hypo_inverter_result->GetName() );
    HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot", plotTitle, hypo_inverter_result);

    // plot in a new canvas with style
    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->SetLogy(false);
    plot->Draw("CLb 2CL"); // plot all and Clb

    // plot test statistics distributions for the two hypothesis
    TCanvas *c2 = new TCanvas("c2", "c2");
    const int nEntries = hypo_inverter_result->ArraySize();
    if (nEntries > 1) {
        int ny = TMath::CeilNint(TMath::Sqrt(nEntries));
        int nx = TMath::CeilNint(double(nEntries) / ny);
        c2->Divide(nx, ny);
    }
    for (int i = 0; i < nEntries; i++) {
        c2->cd(i + 1);
        SamplingDistPlot *pl = plot->MakeTestStatPlot(i);
//        pl->SetXRange(0,100);
//        pl->SetYRange(0,100);
        //pl->SetLogYaxis(true);
        pl->Draw();
    }
}
