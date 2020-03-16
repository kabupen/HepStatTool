
using namespace RooFit;
using namespace RooStats;

class HepMinSvc
{
    public:
        void minuit(RooAbsReal* nll);
        double findSigma( RooAbsReal* nll, double nll_unconditional_fit_value, RooRealVar* poi, const double poi_best_fit, double N);
};

void HepMinSvc::minuit(RooAbsReal* nll)
{
    RooMinimizer m(*nll);
    m.setPrintLevel(-1);
    m.minimize("Minuit2", "Migrad");
}

//get sigma assuming nll -2logLR is poiabolic in poi
double HepMinSvc::findSigma( RooAbsReal* nll, double nll_unconditional_fit_value, RooRealVar* poi, const double poi_best_fit, double N)
{
    double precision = 0.005;

    int direction = int(N/fabs(N));
    double poi_guess = 0;
    if      ( direction > 0 ) poi_guess = poi->getErrorHi();
    else if ( direction < 0 ) poi_guess = poi->getErrorLo();

    bool isConst = poi->isConstant();
    poi->setConstant(true);

    std::map<double, double> guess_to_corr;
    double damping_factor = 1.0;
    double damping_factor_pre = damping_factor;
    double val_pre = poi_guess - 10 * precision;

    double tmu    = 0;
    int nrDamping = 1;
    int nrItr     = 0;

    while ( fabs(val_pre-poi_guess) > precision ) {
        std::cout << "----------------------" << std::endl;
        std::cout << "Starting iteration " << nrItr << " of " << nll->GetName() << " and poiameter " << poi->GetName() << std::endl;
        std::cout << "poi_best_fit = " << poi_best_fit  << std::endl;
        std::cout << "poi_guess    = " << poi_guess << std::endl;
        val_pre = poi_guess;
        damping_factor_pre = damping_factor;

        /* sigma guess */
        if ( poi_guess > 0 && poi->getMax() < poi_guess ) poi->setMax(2*poi_guess);
        if ( poi_guess < 0 && poi->getMin() > poi_guess ) poi->setMin(2*poi_guess);
        poi->setVal(poi_guess);

        // Fits::minimize(nll); /* original */
        HepMinSvc::minuit(nll);

        tmu = 2 * ( nll->getVal() - nll_unconditional_fit_value );

        double sigma_guess = fabs( poi_guess - poi_best_fit )/sqrt(tmu);

        /* Correction factors */
        double corr = damping_factor*(val_pre - poi_best_fit - N*sigma_guess);

        for ( std::map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++ ) {
            if (fabs(itr->first - val_pre) < direction*val_pre*0.02) {
                damping_factor *= 0.8;
                cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
                if (nrDamping++ > 10) {
                    nrDamping = 1;
                    damping_factor = 1.0;
                }
                corr *= damping_factor;
                break;
            }
        }

        //subtract off the difference in the new and damped correction
        poi_guess -= corr;
        guess_to_corr[val_pre] = corr;

        std::cout << "NLL:            " << nll->GetName() << " = " << nll->getVal()  << std::endl;
        std::cout << "delta(NLL):     " << nll->getVal()-nll_unconditional_fit_value << std::endl;
        std::cout << "N*sigma(pre):   " << fabs(val_pre-poi_best_fit)                     << std::endl;
        std::cout << "sigma(guess):   " << sigma_guess                               << std::endl;
        std::cout << "poi(guess):     " << poi_guess+corr                            << std::endl;
        std::cout << "true val:       " << poi_best_fit                                   << std::endl;
        std::cout << "tmu:            " << tmu                                       << std::endl;
        std::cout << "Precision:      " << direction*poi_guess*precision             << std::endl;
        std::cout << "Correction:     " << (-corr<0?" ":"")  << -corr                << std::endl;
        std::cout << "N*sigma(guess): " << fabs(poi_guess-poi_best_fit)                   << std::endl;
        std::cout << std::endl;

        nrItr++;
        if (nrItr > 25) {
            cout << "Infinite loop detected in getSigma(). Please intervene." << endl;
            break;
        }
    }

    std::cout << "Found sigma for nll " << nll->GetName() << ": " << (poi_guess-poi_best_fit)/N << std::endl;
    std::cout << "Finished in " << nrItr << " iterations." << std::endl;
    poi->setConstant(isConst);
    return ( poi_guess - poi_best_fit )/N;
}
