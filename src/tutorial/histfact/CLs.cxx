
using namespace RooFit;
using namespace RooStats;
using namespace std;

// structure defining the options
struct HypoTestInvOptions {

    bool plotHypoTestResult    = true;  // plot test statistic result at each point
    bool writeResult           = true;  // write HypoTestInverterResult in a file
    TString resultFileName;             // file with results (by default is built automatically using the workspace input file name)
    bool optimize              = true;  // optimize evaluation of test statistic
    bool useVectorStore        = true;  // convert data to use new roofit data store
    bool generateBinned        = false; // generate binned data sets
    bool noSystematics         = false; // force all systematics to be off (i.e. set all nuisance parameters as constat to their nominal values)
    double nToysRatio          = 2;     // ratio Ntoys S+b/ntoysB
    double maxPOI              = -1;    // max value used of POI (in case of auto scan)
    bool useProof              = false; // use Proof Lite when using toys (for freq or hybrid)
    int nworkers               = 0;     // number of worker for ProofLite (default use all available cores)
    bool enableDetailedOutput  = false; // enable detailed output with all fit information for each toys (output will be written in result file)
    bool rebuild               = false; // re-do extra toys for computing expected limits and rebuild test stat distributions (N.B this requires much more CPU (factor is equivalent to nToyToRebuild)
    int nToyToRebuild          = 100;   // number of toys used to rebuild
    int rebuildParamValues     = 0;     // = 0   do a profile of all the parameters on the B (alt snapshot) before performing a
                                        // rebuild operation (default)
                                        // = 1   use initial workspace parameters with B snapshot values
                                        // = 2   use all initial workspace parameters with B
                                        // Otherwise the rebuild will be performed using
    int initialFit = -1;                // do a first  fit to the model (-1 : default, 0 skip fit, 1 do always fit)
    int randomSeed = -1;                // random seed (if = -1: use default value, if = 0 always random )
    // NOTE: Proof uses automatically a random seed

    int nAsimovBins = 0; // number of bins in observables used for Asimov data sets (0 is the default and it is given by
    // workspace, typically is 100)

    bool reuseAltToys         = false; // reuse same toys for alternate hypothesis (if set one gets more stable bands)
    double confLevel          = 0.95;   // confidence level value
    std::string minimizerType = ""; // minimizer type (default is what is in ROOT::Math::MinimizerOptions::DefaultMinimizerType()
    std::string massValue     = "";     // extra string to tag output file of result
    int printLevel            = 0;             // print level for debugging PL test statistics and calculators
    bool useNLLOffset         = false;      // use NLL offset when fitting (this increase stability of fits)
};

void CLs()
{
    const char *infile = 0; 
    const char *wsName = "combined";
    const char *modelSBName = "ModelConfig";
    const char *modelBName = "";
    const char *dataName = "obsData"; 
    int calculatorType = 3; 
    int testStatType = 3;
    bool useCLs = true; 
    int npoints = 6; 
    double poimin = 0; 
    double poimax = 5;
    int ntoys = 1000; 
    bool useNumberCounting = false; 
    const char *nuisPriorName = 0;
    /*
       Other Parameter to pass in tutorial
       apart from standard for filename, ws, modelconfig and data

       type = 0 Freq calculator
       type = 1 Hybrid calculator
       type = 2 Asymptotic calculator
       type = 3 Asymptotic calculator using nominal Asimov data sets (not using fitted parameter values but nominal ones)

       testStatType = 0 LEP
       = 1 Tevatron
       = 2 Profile Likelihood
       = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)
       = 4 Profiel Likelihood signed ( pll = -pll if mu < mu_hat)
       = 5 Max Likelihood Estimate as test statistic
       = 6 Number of observed event as test statistic

       useCLs          scan for CLs (otherwise for CLs+b)

       npoints:        number of points to scan , for autoscan set npoints = -1
       
       poimin,poimax:  min/max value to scan in case of fixed scans
       (if min >  max, try to find automatically)
       
       ntoys:         number of toys to use
       
       useNumberCounting:  set to true when using number counting events
       
       nuisPriorName:   name of prior for the nuisance. This is often expressed as constraint term in the global model
       It is needed only when using the HybridCalculator (type=1)
       If not given by default the prior pdf from ModelConfig is used.
       
       extra options are available as global parameters of the macro. They major ones are:
       
       plotHypoTestResult   plot result of tests at each point (TS distributions) (default is true)
       useProof             use Proof   (default is true)
       writeResult          write result of scan (default is true)
       rebuild              rebuild scan for expected limits (require extra toys) (default is false)
       generateBinned       generate binned data sets for toys (default is false) - be careful not to activate with
       a too large (>=3) number of observables
       nToyRatio            ratio of S+B/B toys (default is 2)
    */

    TString filename(infile);
    if (filename.IsNull()) {
        filename = "results/example_combined_GaussExample_model.root";
        bool fileExist = !gSystem->AccessPathName(filename); // note opposite return code

        // if file does not exists generate with histfactory
        if (!fileExist) {
            // Normally this would be run on the command line
            cout << "will run standard hist2workspace example" << endl;
            gROOT->ProcessLine(".! prepareHistFactory .");
            gROOT->ProcessLine(".! hist2workspace config/example.xml");
            cout << "\n\n---------------------" << endl;
            cout << "Done creating example input" << endl;
            cout << "---------------------\n\n" << endl;
        }
    }
    // Try to open the file
    TFile *file = TFile::Open(filename);

    // if input file was specified byt not found, quit
    if (!file) {
        cout << "StandardRooStatsDemoMacro: Input file " << filename << " is not found" << endl;
        return;
    }

    RooWorkspace *w = dynamic_cast<RooWorkspace *>(file->Get(wsName));
    HypoTestInverterResult *r = 0;
    std::cout << w << "\t" << filename << std::endl;

    std::cout << "Running HypoTestInverter on the workspace " << w->GetName() << std::endl;

    w->Print();

    RooAbsData *data = w->data(dataName);
    std::cout << "Using data set " << dataName << std::endl;

    // (mUseVectorStore = true
    RooAbsData::setDefaultStorageType(RooAbsData::Vector);
    data->convertToVectorStore();

    // get models from WS
    // get the modelConfig out of the file
    std::cout << "Get models modelBName=" << modelBName << " modelSBName=" << modelSBName << std::endl;
    ModelConfig *bModel  = (ModelConfig *)w->obj(modelBName);
    ModelConfig *sbModel = (ModelConfig *)w->obj(modelSBName);

    if (!sbModel) { return 0; }
    // check the model
    if (!sbModel->GetPdf())                  { Error("StandardHypoTestDemo", "Model %s has no pdf ", modelSBName); return 0; }
    if (!sbModel->GetParametersOfInterest()) { Error("StandardHypoTestDemo", "Model %s has no poi ", modelSBName); return 0; }
    if (!sbModel->GetObservables())          { Error("StandardHypoTestInvDemo", "Model %s has no observables ", modelSBName); return 0; }
    if (!sbModel->GetSnapshot())             { Info("StandardHypoTestInvDemo", "Model %s has no snapshot  - make one using model poi", modelSBName); sbModel->SetSnapshot(*sbModel->GetParametersOfInterest()); }

    // if (!bModel || bModel == sbModel) 
    std::cout << "The background model " <<  modelBName << " does not exist" << std::endl;
    std::cout << "Copy it from ModelConfig " << modelSBName << " and set POI to zero" << std::endl;
    bModel = (ModelConfig *)sbModel->Clone();
    bModel->SetName(TString(modelSBName) + TString("_with_poi_0"));
    RooRealVar *var = dynamic_cast<RooRealVar *>(bModel->GetParametersOfInterest()->first());
    if (!var) return 0;
    double oldval = var->getVal();
    var->setVal(0);
    bModel->SetSnapshot(RooArgSet(*var));
    var->setVal(oldval);

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
    const RooArgSet *poiSet = sbModel->GetParametersOfInterest();
    RooRealVar *poi = (RooRealVar *)poiSet->first();
    std::cout << "StandardHypoTestInvDemo : POI initial value:   " << poi->GetName() << " = " << poi->getVal() << std::endl;

    double poihat = 0;

    std::string minimizer_type = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
    std::cout << "Using as minimizer for computing the test statistic"<< minimizer_type << std::endl;

    // build test statistics and hypotest calculators for running the inverter
    SimpleLikelihoodRatioTestStat slrts(*sbModel->GetPdf(), *bModel->GetPdf());

    // null parameters must includes snapshot of poi plus the nuisance values
    RooArgSet nullParams(*sbModel->GetSnapshot());
    if (sbModel->GetNuisanceParameters()) nullParams.add(*sbModel->GetNuisanceParameters());
    if (sbModel->GetSnapshot())           slrts.SetNullParameters(nullParams);

    RooArgSet altParams(*bModel->GetSnapshot());
    if (bModel->GetNuisanceParameters()) altParams.add(*bModel->GetNuisanceParameters());
    if (bModel->GetSnapshot())           slrts.SetAltParameters(altParams);

    AsymptoticCalculator::SetPrintLevel(0);

    // create the HypoTest calculator class
    // Asimov data generated with nominal values
    HypoTestCalculatorGeneric *hc = new AsymptoticCalculator(*data, *bModel, *sbModel, true); // for using Asimov data generated with nominal values

    ((AsymptoticCalculator *)hc)->SetOneSided(true);

    // Get the result
    RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);

    HypoTestInverter hypo(*hc);
    hypo.SetConfidenceLevel(0.95);
    hypo.UseCLs(true);
    hypo.SetVerbose(true);
    hypo.SetFixedScan(npoints, poimin, poimax);

    r = hypo.GetInterval();
    if (!r) { std::cerr << "Error running the HypoTestInverter - Exit " << std::endl; return; }

    // analyze result produced by the inverter, optionally save it in a file

    double lowerLimit = r->LowerLimit();
    double llError    = r->LowerLimitEstimatedError();
    double upperLimit = r->UpperLimit();
    double ulError    = r->UpperLimitEstimatedError();

    if (lowerLimit < upperLimit * (1. - 1.E-4) && lowerLimit != 0)
        std::cout << "The computed lower limit is: " << lowerLimit << " +/- " << llError << std::endl;
    std::cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError << std::endl;

    // compute expected limit
    std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
    std::cout << " expected limit (median) " << r->GetExpectedUpperLimit(0) << std::endl;
    std::cout << " expected limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << std::endl;
    std::cout << " expected limit (+1 sig) " << r->GetExpectedUpperLimit(1) << std::endl;
    std::cout << " expected limit (-2 sig) " << r->GetExpectedUpperLimit(-2) << std::endl;
    std::cout << " expected limit (+2 sig) " << r->GetExpectedUpperLimit(2) << std::endl;

    // write result in a file
    if ( r != NULL ) {

        // write to a file the results
        const char *calcType  = "Asym";
        const char *limitType = "CLs";
        const char *scanType  = (npoints < 0) ? "auto" : "grid";

        TString outputname = TString::Format("%s_%s_%s_ts%d.root", calcType, limitType, scanType, testStatType);
        TFile *fileOut = new TFile( outputname, "RECREATE");
        r->Write();
        fileOut->Close();
    }

    // plot the result ( p values vs scan points)
    std::string typeName = "Asymptotic";

    const char *resultName = r->GetName();
    TString plotTitle = TString::Format("%s CL Scan for workspace %s", typeName.c_str(), resultName);
    HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot", plotTitle, r);

    // plot in a new canvas with style
    TString c1Name = TString::Format("%s_Scan", typeName.c_str());
    TCanvas *c1 = new TCanvas(c1Name);
    c1->SetLogy(false);

    plot->Draw("CLb 2CL"); // plot all and Clb

    const int nEntries = r->ArraySize();

    // plot test statistics distributions for the two hypothesis
    TCanvas *c2 = new TCanvas("c2");
    if (nEntries > 1) {
        int ny = TMath::CeilNint(TMath::Sqrt(nEntries));
        int nx = TMath::CeilNint(double(nEntries) / ny);
        c2->Divide(nx, ny);
    }
    for (int i = 0; i < nEntries; i++) {
        if (nEntries > 1)
            c2->cd(i + 1);
        SamplingDistPlot *pl = plot->MakeTestStatPlot(i);
        pl->SetLogYaxis(true);
        pl->Draw();
    }
    gPad = c1;
}
