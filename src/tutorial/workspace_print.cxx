
using namespace RooFit; 
using namespace RooStats; 

void workspace_print()
{
    //construct the model
    RooWorkspace w("w");

    w.factory("Gaussian::constraint_b(nuisance_b[41.7,0,100],41.7,4.6)"); //constrained b to be positive - "truncated gaussian"
    w.factory("Gaussian::constraint_acc(nuisance_acc[0.71,0,1],0.71,0.09)"); //constrained acc in range 0-1
    w.factory("Gaussian::constraint_lumi(nuisance_lumi[5.0,0.0,10.0],5.0,0.195)"); //constrained lumi from 0 to 10.0
    w.factory("prod::s(sigma[0,100],nuisance_lumi,nuisance_acc)");
    w.factory("sum::mean(s,nuisance_b)");
    w.factory("Poisson::pois(n[61,0,100],mean)");
    w.factory("PROD::model(pois,constraint_b,constraint_lumi,constraint_acc)");
    
    //define RooArgSets for convenience
    w.defineSet("obs","n"); //observables
    w.defineSet("poi","sigma"); //parameters of interest
    w.defineSet("np","nuisance_b,nuisance_lumi,nuisance_acc"); //nuisance parameters

    RooDataSet data("data", "data", *w.set("obs"));
    data.add(*w.set("obs")); //actually add the data

    w.Print();
    w.var("sigma")->Print();
    
    return;

    RooAbsReal* nll = w.pdf("model")->createNLL(data);
    w.pdf("model")->fitTo(data,Minos(*w.set("poi")),Save(),Hesse(false));
    std::cout << "@@@@ " << w.var("sigma")->getErrorHi() << " " << w.var("sigma")->getErrorLo() << std::endl;

    RooAbsReal* pll = nll->createProfile(*w.set("poi"));
    w.pdf("model")->fitTo(data,Minos(*w.set("poi")),Save(),Hesse(false));
    std::cout << "@@@@ " << w.var("sigma")->getErrorHi() << " " << w.var("sigma")->getErrorLo() << std::endl;
    
    RooPlot* frame = w.var("sigma")->frame();
    nll->plotOn(frame, ShiftToZero()); //the ShiftToZero option puts the minimum at 0 on the y-axis
    frame->Draw();
    pll->plotOn(frame,LineColor(kRed));
    frame->Draw();
    
    RooFitResult* res = w.pdf("model")->fitTo(data,Minos(*w.set("poi")),Save(),Hesse(false));
    

    if(res->status()==0) {
        w.var("sigma")->Print();
    } else {
        cout << "Likelihood maximization failed" << endl;
    }

}
