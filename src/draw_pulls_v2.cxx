// Author      : Stefan Gadatsch
// Email       : gadatsch@nikhef.nl
// Date        : 2013-04-13
// Description : Draw pulls of nuisance parameters, rank by importance

#include "TCanvas.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TSystemFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TROOT.h"
#include <TPRegexp.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>

#include "../HepStatTool/PullPlotHelper.h"
#include "../HepStatTool/drawPlot.C"
#include "root2ascii.cxx"

// global style options
TString energy1 = "8";
TString energy2 = "13";
TString lumi1 = "20.3";
TString lumi2 = "36.1";
TString lumi3 = "43.8";
TString lumi4 = "79.8";
TString lumi5 = "139.0";
bool doHorizontal          = false; // produce a horizontal plot
bool drawErrorBars         = false; // draw bars visualising the total, stat and syst uncertainty
bool drawHatchedBands      = false; // draw hatched bands around delta muhat = 0 to visualise the total, stat and syst uncertainty
bool drawParamNames        = true;  // show the nuisance parameter labels
bool drawPostfitImpactBand = true;  // && (mode != error); // draw a band with the variation of muhat due to each theta (postfit uncertainty)
bool drawPrefitImpactBand  = true;  // && (mode != error); // draw a band with the variation of muhat due to each theta (prefit uncertainty)
bool drawStandardBand      = false;  // draw 1 sigma standard deviation bar
bool drawRedNorm           = true;  //draw the Normalization parameters with Red dots
bool useRelativeImpact     = false;  // switch to delta muhat / delta muhat tot for the top axis
bool scaleGammas           = false; // switch to scale gamma nuisance parameters to 1/MCStat.error -1
int useBreakdown           = 0;     // 0 = add, 1 = sub
//double scale_poi           = 1.25;  // zoom the impact axis
bool isVH                  = true; // switch only when using the flag removeSigUnc. True is VH, false is VZ 
bool removeSigUnc          = false; // remove signal uncertainties from the plot
int showTopParameters      = 15; // -1 to show all parameters
double showHighImpact      = 0.0;   // sigma_comp / sigma_tot threshold
double UserlabelSize       = 0.03; // label size
double hatch_width         = 2.5; // postfit hatch width
Color_t color_standardband    = kRed-4;
Color_t color_totalerror   = kBlue-4;
Color_t color_staterror    = kGreen+1;
Color_t color_systerror    = kMagenta-4;
Color_t color_pulls        = kGray+2;
Color_t color_normalization = kRed-4;
Color_t color_prefit       = kYellow-7;
Color_t color_postfit      = kBlue-4;
bool m_postFitOrder          = true;
bool rankPOI_top           = false;

vector<string> getLabel(const char* fileName, int nrPars);
TString translateNPname(TString internalName, bool isMVA);
TString translateGammaStatName(TString internalName);
bool isSignalNP(string NPname);

void draw_pulls_core(
        std::string mass = "125", 
        std::string cardName = "", 
        float scale_factor = 1.7,  
        std::string overlayCard="", 
        bool postFitOrder = true, 
        std::string POIname = "SigXsecOverSM", 
        int POIpos = 0, 
        int POItotal = 1) 
{
    gStyle->SetHatchesLineWidth(hatch_width);

    m_postFitOrder = postFitOrder;
    std::vector<std::string> parsed = parseString(cardName, ":");
    string cardOpts;
    if (parsed.size() > 1) {
        cardOpts = parsed[1];
    }
    cardName = parsed[0];
    computeFlags(cardName);

    applyOverlay(overlayCard , overlay , "");

    TCanvas* c1 = new TCanvas("c1","c1",1024,1448);
    TPad *pad1  = new TPad("pad1", "pad1",  0.0, 0.0,  1.0,  1.0, 0);

    /* pad1 */
    if (drawParamNames) pad1->SetLeftMargin(0.30);//Graph size
    else                pad1->SetLeftMargin(0.05);
    if (drawErrorBars)  pad1->SetTopMargin(0.10);
    else                pad1->SetTopMargin(0.09);
    pad1->SetRightMargin(0.05);
    pad1->SetBottomMargin(0.09);
    pad1->Draw();

    labelPosX   = 0.06;
    channelPosX = 0.33;
    channelPosY = 0.19;
    markerSize  = 0.8;

    // Start (Old draw pulls2)
    std::cout << "INFO::Drawing pulls: " << cardName << " for mH = " << mass << " GeV" << std::endl;

    gStyle->SetHatchesLineWidth(hatch_width);

    // load and initialize ascii files
    std::string file_pulls             = "./data/output/" + cardName + "/ascii/pulls.txt";
    std::string file_breakdown_add     = "./data/output/" + cardName + "/ascii/breakdown_add.txt";
    std::string file_pulls_id          = "./data/output/" + cardName + "/ascii/pulls_id.txt";
    std::string file_breakdown_add_id  = "./data/output/" + cardName + "/ascii/breakdown_add_id.txt";

    FileHolder* pulls = new FileHolder( file_pulls         , 3+POItotal*5);
    FileHolder* cats  = new FileHolder( file_breakdown_add , POItotal*3  );

    // get the values from the ascii files
    int nrNuis = pulls->getMassPointsSize();
    int nrCats = cats ->getMassPointsSize();
    std::cout << "# nuisance parameters = " << nrNuis << std::endl;
    std::cout << "# nrCats = " << nrCats << std::endl;

    std::vector<double> points_nuis = pulls->getMassPoints();
    std::vector<double> points_cats = cats ->getMassPoints();

    for (int i = 0; i < nrNuis; i++) points_nuis[i] = i + 0.5;

    std::vector<double> val               = pulls->getCol(0);
    std::vector<double> up                = pulls->getCol(1);
    std::vector<double> down              = pulls->getCol(2);
    std::vector<double> poi_hat           = pulls->getCol(3+5*POIpos);
    std::vector<double> poi_up            = pulls->getCol(4+5*POIpos);
    std::vector<double> poi_down          = pulls->getCol(5+5*POIpos);
    std::vector<double> poi_nom_up        = pulls->getCol(6+5*POIpos);
    std::vector<double> poi_nom_down      = pulls->getCol(7+5*POIpos);

    std::vector<double> poi_up_sign       = pulls->getCol(4+5*POIpos);
    std::vector<double> poi_down_sign     = pulls->getCol(5+5*POIpos);
    std::vector<double> poi_nom_up_sign   = pulls->getCol(6+5*POIpos);
    std::vector<double> poi_nom_down_sign = pulls->getCol(7+5*POIpos);

    std::vector<double> cats_val          = cats->getCol(0+3*POIpos);
    std::vector<double> cats_up           = cats->getCol(1+3*POIpos);
    std::vector<double> cats_down         = cats->getCol(2+3*POIpos);

    std::cout << "Set correct values for the poi" << std::endl;
    float scale_theta = scale_factor;
    float scale_poi   = scale_factor;

    for (int i = 0; i < nrNuis; i++) {
        val[i] *= scale_theta;

        poi_up[i]   = poi_up[i]   - poi_hat[i];
        poi_down[i] = poi_down[i] - poi_hat[i];

        poi_nom_up[i]   = poi_nom_up[i]   - poi_hat[i];
        poi_nom_down[i] = poi_nom_down[i] - poi_hat[i];
        //
        if (poi_up[i] < 0) {
            poi_up_sign[i]       = 0.;
            poi_down_sign[i]     = fabs(poi_up[i]);
            poi_nom_up_sign[i]   = 0.;
            poi_nom_down_sign[i] = fabs(poi_nom_up[i]);
        } else {
            poi_up_sign[i]       = fabs(poi_up[i]);
            poi_down_sign[i]     = 0.;
            poi_nom_up_sign[i]   = fabs(poi_nom_up[i]);
            poi_nom_down_sign[i] = 0.;
        }

        if ( poi_up[i]     < 0 ) std::swap(poi_up[i]    , poi_down[i]);
        if ( poi_nom_up[i] < 0 ) std::swap(poi_nom_up[i], poi_nom_down[i]);

        poi_up[i]       = fabs(poi_up[i]);
        poi_down[i]     = fabs(poi_down[i]);
        poi_nom_up[i]   = fabs(poi_nom_up[i]);
        poi_nom_down[i] = fabs(poi_nom_down[i]);
        poi_hat[i]      = 0;
    }

    std::cout << "Do the sum for printout at the end" << std::endl;
    double sum_poi2 = 0.;
    for (int i = 0; i < nrNuis; ++i) {  
        double up   = fabs(poi_up[i] - poi_hat[i]);
        double down = fabs(poi_up[i] - poi_hat[i]);
        sum_poi2 += pow((up+down)/2, 2);
    }

    std::cout << "Find maximal error due to a single nuisance parameter" << std::endl;
    double max_poi = 0.;
    for (int i = 0; i < nrNuis; ++i) {
        if (poi_up[i]   > max_poi) max_poi = poi_up[i];
        if (poi_down[i] > max_poi) max_poi = poi_down[i];
    }

    // get labels
    std::vector<std::string> labels;
    std::vector<std::string> cats_labels;
    int nlines      = 0;
    int cats_nlines = 0;

    std::ifstream idFile ( file_pulls_id );
    std::ifstream idFile3( file_breakdown_add_id ); 

    while (1) {
        if (!idFile.good() || nlines > nrNuis-1) break;
        std::string label;
        idFile >> label;
        labels.push_back(label);
        cout << "added: " << label << endl;
        nlines++;
    }
    while (1) {
        if (!idFile3.good() || cats_nlines > nrCats-1) break;
        string cat_label;
        idFile3 >> cat_label;
        cats_labels.push_back(cat_label);
        cats_nlines++;
    }

    std::cout << "Map of category uncertainties" << std::endl;
    std::map<string, vector<double> > cat_uncerts;
    for (int i = 0; i < nrCats; i++) {
        string index = cats_labels[i];
        cout << i << " " << index << " " << cats_val[i] << " " << cats_up[i] << " " << cats_down[i] << endl;
        cat_uncerts[index].push_back(cats_val[i]);
        cat_uncerts[index].push_back(cats_up[i]);
        cat_uncerts[index].push_back(cats_down[i]);
    }

    double sigma_tot_hi  = cat_uncerts["total"][1];
    double sigma_tot_lo  = cat_uncerts["total"][2];

    double sigma_stat_hi = 0.;
    double sigma_stat_lo = 0.;
    double sigma_syst_hi = 0.;
    double sigma_syst_lo = 0.;
    // TODO: can probably drop this

    std::cout << "Dump everything in maps" << std::endl;
    std::map< std::string, std::vector<double> > nuis_map;
    for ( int i = 0; i < nrNuis; i++) {
        std::string index = labels[i];      
        nuis_map[index].push_back(val[i]);
        nuis_map[index].push_back(up[i]);
        nuis_map[index].push_back(down[i]);
        nuis_map[index].push_back(poi_hat[i]);
        nuis_map[index].push_back(poi_up[i]);
        nuis_map[index].push_back(poi_down[i]);
        nuis_map[index].push_back(poi_nom_up[i]);
        nuis_map[index].push_back(poi_nom_down[i]);
    }

    std::cout << "Dump everything in maps" << std::endl;
    std::map<std::string, std::vector<double> > nuis_map_sign;
    for (int i = 0; i < nrNuis; i++) {
        string index = labels[i];
        nuis_map_sign[index].push_back(val[i]);
        nuis_map_sign[index].push_back(up[i]);
        nuis_map_sign[index].push_back(down[i]);
        nuis_map_sign[index].push_back(poi_hat[i]);
        nuis_map_sign[index].push_back(poi_up_sign[i]);
        nuis_map_sign[index].push_back(poi_down_sign[i]);
        nuis_map_sign[index].push_back(poi_nom_up_sign[i]);
        nuis_map_sign[index].push_back(poi_nom_down_sign[i]);
    }

    // Getting the vectors back
    nrNuis    = labels.size();

    for (int i = 0; i < nrNuis-1; i++) {
        for (int j = 0; j < nrNuis-1-i; j++) {
            if (strcmp(labels[i].c_str(),labels[i+1].c_str())) {
                std::swap(labels[j], labels[j+1]);
            }
        }
    }

    for (int i = 0; i < nrNuis; i++) {
        if(isSignalNP(labels[i])){
            val[i]          = 0;
            up[i]           = 0;
            down[i]         = 0;
            poi_hat[i]      = 0;
            poi_up[i]       = 0;
            poi_down[i]     = 0;
            poi_nom_up[i]   = 0;
            poi_nom_down[i] = 0;
            poi_up_sign[i]       = 0;
            poi_down_sign[i]     = 0;
            poi_nom_up_sign[i]   = 0;
            poi_nom_down_sign[i] = 0;
        }else{
            val[i]          = nuis_map[labels[i]][0];
            up[i]           = nuis_map[labels[i]][1];
            down[i]         = nuis_map[labels[i]][2];
            poi_hat[i]      = nuis_map[labels[i]][3];
            poi_up[i]       = nuis_map[labels[i]][4];
            poi_down[i]     = nuis_map[labels[i]][5];
            poi_nom_up[i]   = nuis_map[labels[i]][6];
            poi_nom_down[i] = nuis_map[labels[i]][7];
            poi_up_sign[i]       = nuis_map_sign[labels[i]][4];
            poi_down_sign[i]     = nuis_map_sign[labels[i]][5];
            poi_nom_up_sign[i]   = nuis_map_sign[labels[i]][6];
            poi_nom_down_sign[i] = nuis_map_sign[labels[i]][7];
        }
    }

    // sort poi values by variation size
    int test_idx=-1;
    int NPoi=POItotal;
    for (int i = 0; i < nrNuis-NPoi; i++) {
        for (int j = 0; j < nrNuis-NPoi-i; j++) {
            bool doSwap = false;
            if (m_postFitOrder) { doSwap = poi_up[j]+poi_down[j] > poi_up[j+1]+poi_down[j+1]; } 
            else                { doSwap = poi_nom_up[j]+poi_nom_down[j] > poi_nom_up[j+1]+poi_nom_down[j+1]; }
            if (doSwap) {
                // swap postfit poi
                swap(poi_up[j], poi_up[j+1]);
                swap(poi_down[j], poi_down[j+1]);
                swap(poi_up_sign[j], poi_up_sign[j+1]);
                swap(poi_down_sign[j], poi_down_sign[j+1]);

                // swap prefit poi
                swap(poi_nom_up[j], poi_nom_up[j+1]);
                swap(poi_nom_down[j], poi_nom_down[j+1]);

                // swap pulls
                swap(up[j], up[j+1]);
                swap(down[j], down[j+1]);
                swap(val[j], val[j+1]);

                // swap names
                swap(labels[j], labels[j+1]);
            }
        }
    }

    // make the 1 sigma boxes
    std::vector<double> boxup;
    std::vector<double> boxdown;
    std::vector<double> cenup;
    std::vector<double> cendown;

    for (int i = 0; i < nrNuis; i++) {
        boxup.push_back(1.*scale_theta);
        boxdown.push_back(1.*scale_theta);
        double height = 0.5;
        cenup.push_back(height);
        cendown.push_back(height);
    }

    // make the 1 sigma boxes
    double* statboxup   = new double[nrNuis];
    double* statboxdown = new double[nrNuis];
    double* systboxup   = new double[nrNuis];
    double* systboxdown = new double[nrNuis];

    for (int i = 0; i < nrNuis; i++) {
        statboxup[i]   = sigma_stat_hi * scale_poi / max_poi;
        statboxdown[i] = sigma_stat_lo * scale_poi / max_poi;

        systboxup[i]   = sigma_syst_hi * scale_poi / max_poi;
        systboxdown[i] = sigma_syst_lo * scale_poi / max_poi;
    }

    // make the final arrays for plotting, in particular remove parameters
    int nrNuis2remove = 0;
    for (int i = 0; i < nrNuis; i++) {
        // print up and down effect
        cout << "Rank " << nrNuis - i << ":  \t" << fabs(poi_down[i]-poi_hat[i]) << "  \t" << fabs(poi_up[i]-poi_hat[i]) << "  \t" << labels[i] << endl;
        if ( fabs(poi_down[i]-poi_hat[i]) > 150 ) cout << " --  " << poi_down[i] << " " << poi_hat[i] << " " << poi_up[i] <<endl;

        if ((fabs(poi_down[i]) + fabs(poi_up[i])) / (sigma_tot_lo + sigma_tot_hi) < showHighImpact ){
            cout << "WARNING::Removing " << labels[i] << ". Below threshold." << endl;
            nrNuis2remove++;
        }
    }

    if (showTopParameters != -1) nrNuis2remove = std::max(0, nrNuis - showTopParameters);

    labels.erase(labels.begin(), labels.begin() + nrNuis2remove);
    points_nuis.erase(points_nuis.end() - nrNuis2remove, points_nuis.end());

    val.erase(val.begin(), val.begin() + nrNuis2remove);
    down.erase(down.begin(), down.begin() + nrNuis2remove);
    up.erase(up.begin(), up.begin() + nrNuis2remove);

    poi_hat      .erase(poi_hat.begin()      , poi_hat      .begin() + nrNuis2remove );
    poi_down     .erase(poi_down.begin()     , poi_down     .begin() + nrNuis2remove );
    poi_up       .erase(poi_up.begin()       , poi_up       .begin() + nrNuis2remove );
    poi_down_sign.erase(poi_down_sign.begin(), poi_down_sign.begin() + nrNuis2remove );
    poi_up_sign  .erase(poi_up_sign.begin()  , poi_up_sign  .begin() + nrNuis2remove );
    poi_nom_down .erase(poi_nom_down.begin() , poi_nom_down .begin() + nrNuis2remove );
    poi_nom_up   .erase(poi_nom_up.begin()   , poi_nom_up   .begin() + nrNuis2remove );
    boxdown      .erase(boxdown.begin()      , boxdown      .begin() + nrNuis2remove );
    boxup        .erase(boxup.begin()        , boxup        .begin() + nrNuis2remove );
    cendown      .erase(cendown.begin()      , cendown      .begin() + nrNuis2remove );
    cenup        .erase(cenup.begin()        , cenup        .begin() + nrNuis2remove );

    nrNuis -= nrNuis2remove;

    int offset = ceil(2 * nrNuis / 10); // used for space to plot the labels and legend

    for (int i = 0; i < nrNuis; i++) {
        poi_up[i] = fabs(poi_up[i]) * scale_poi / max_poi;
        poi_down[i] = fabs(poi_down[i]) * scale_poi / max_poi;
        poi_up_sign[i] = fabs(poi_up_sign[i]) * scale_poi / max_poi;
        poi_down_sign[i] = fabs(poi_down_sign[i]) * scale_poi / max_poi;

        poi_nom_up[i] = fabs(poi_nom_up[i]) * scale_poi / max_poi;
        poi_nom_down[i] = fabs(poi_nom_down[i]) * scale_poi / max_poi;

        if( labels[i].find("gamma") != std::string::npos ){
            poi_nom_up[i] = 0;
            poi_nom_down[i] = 0;
        }

        if (useRelativeImpact) {
            poi_up[i] /= sigma_tot_hi;
            poi_down[i] /= sigma_tot_lo;
            poi_up_sign[i] /= sigma_tot_hi;
            poi_down_sign[i] /= sigma_tot_lo;
            poi_nom_up[i] /= sigma_tot_hi;
            poi_nom_down[i] /= sigma_tot_lo;
        }

        up[i] = fabs(up[i]) * scale_theta;
        down[i] = fabs(down[i]) * scale_theta;
    }

    // change to the right pad
    pad1->cd();
    // make plot of pulls for nuisance parameters
    markerSize = 2;

    // make plot of 1 sigma boxes
    PullPlotHelper pull_helper;
    TGraphAsymmErrors* gr1s = pull_helper.makeGraphErr("", nrNuis, &val[0], &points_nuis[0], &boxdown[0], &boxup[0], NULL, NULL);
    gr1s->SetLineColor(color_standardband);
    gr1s->SetMarkerColor(color_standardband);
    gr1s->SetLineStyle(1);
    gr1s->SetLineWidth(2);
    gr1s->SetMarkerSize(markerSize);
    gr1s->GetXaxis()->SetTitleOffset(1.2);

    // make plot for the POI change for postfit uncertainties
    double hatch_width=4;
    TGraphAsymmErrors* gr_poi = pull_helper.makeGraphErr("", nrNuis, &poi_hat[0], &points_nuis[0], &poi_down[0], &poi_up[0], &cenup[0], &cendown[0]);
    gr_poi->SetLineColor(color_postfit);
    gr_poi->SetFillColor(color_postfit);
    gr_poi->SetFillStyle(0);
    gr_poi->SetLineWidth(2);
    gr_poi->SetMarkerSize(0);
    
    // make plot for the POI change for postfit uncertainties positive side
    vector<double> poi_zero;
    for(int ii=0;ii<nrNuis;ii++)
        poi_zero.push_back(0.);
    TGraphAsymmErrors* gr_poi_pos = pull_helper.makeGraphErr("", nrNuis, &poi_hat[0], &points_nuis[0], &poi_down_sign[0], &poi_up_sign[0], &cenup[0], &cendown[0]);
    gr_poi_pos->SetLineColor(color_postfit);
    gr_poi_pos->SetFillColor(color_postfit);
    gr_poi_pos->SetFillStyle(3354);//3004);
    gr_poi_pos->SetLineWidth(2);
    gr_poi_pos->SetMarkerSize(0);

    // make plot for the POI change for prefit uncertainties
    TGraphAsymmErrors* gr_poi_nom = pull_helper.makeGraphErr("", nrNuis, &poi_hat[0], &points_nuis[0], &poi_nom_down[0], &poi_nom_up[0], &cenup[0], &cendown[0]);
    gr_poi_nom->SetLineColor(color_prefit);
    gr_poi_nom->SetFillColor(color_prefit);
    gr_poi_nom->SetLineWidth(1);
    gr_poi_nom->SetMarkerSize(0);

    double border_lo = -sigma_tot_lo / max_poi;
    double border_hi = sigma_tot_hi / max_poi;
    // different shades for better readability
    int nrShades = ceil((nrNuis+1)/2);
    std::vector<double> shadeCenter;
    std::vector<double>  shadePoints;
    std::vector<double>  shadeWidth;
    std::vector<double>  shadeHeight;
    for (int ishade = 0; ishade < nrShades; ishade++) {
        shadeCenter.push_back(0.0);
        shadePoints.push_back(2.0*ishade+0.5);
        shadeWidth.push_back(10.); // TODO: should not be hardcoded
        shadeHeight.push_back(0.5);
    }

    TGraphAsymmErrors* gr_shades = pull_helper.makeGraphErr("", nrShades, &shadeCenter[0], &shadePoints[0], &shadeWidth[0], &shadeWidth[0], &shadeHeight[0], &shadeHeight[0]);
    gr_shades->SetLineColor(18);
    gr_shades->SetFillColor(18);
    gr_shades->SetFillStyle(1001);
    gr_shades->SetLineWidth(1);
    gr_shades->SetMarkerSize(0);

    //Added
    TH2F *h = new TH2F("h", "", 1, border_lo, border_hi, nrNuis+offset+1, -offset, nrNuis+1);
    vector<int> isNorm;
    for (int i = offset; i < nrNuis+offset; i++){
        if(labels[i-offset].find("norm")!= string::npos)
            isNorm.push_back(1);
        else
            isNorm.push_back(0);
        bool isBDT = false;
        if(cardName.find("MVA") != string::npos)
            isBDT = true;
        TString tmpName = labels[i-offset];
        bool is8TeV = tmpName.Contains("8TeV");
        tmpName = tmpName.ReplaceAll("_8TeV","");
        TString newLabels = translateNPname(tmpName, isBDT);
        cout<<(i-offset)<<"th NP Renamed "<<labels[i-offset]<<" -> "<<newLabels<<endl;
        h->GetYaxis()->SetBinLabel(i+1, drawParamNames?newLabels.Data():"");
    }

    std::vector<double> points_nuis_norm;
    std::vector<double> val_norm;
    std::vector<double> up_norm;
    std::vector<double> down_norm;
    std::vector<double> points_nuis_pull;
    std::vector<double> val_pull;
    std::vector<double> up_pull;
    std::vector<double> down_pull;
    for(int i = 0;i < nrNuis; i++){
        if(isNorm.at(i) == 1){
            points_nuis_norm.push_back(points_nuis.at(i));
            val_norm.push_back(val.at(i));
            up_norm.push_back(up.at(i));
            down_norm.push_back(down.at(i));
        }
        else{
            points_nuis_pull.push_back(points_nuis.at(i));
            if(labels[i].find("gamma_stat")!= std::string::npos && scaleGammas){
                double fit = ((val.at(i)/scale_theta) - 1) / (up.at(i)/scale_theta);
                val_pull.push_back(fit*scale_theta);
                up_pull.push_back(scale_theta);
                down_pull.push_back(scale_theta);

            }else{
                val_pull.push_back(val.at(i));
                up_pull.push_back(up.at(i));
                down_pull.push_back(down.at(i));
            }
        }
    }

    TGraphAsymmErrors* gr = pull_helper.makeGraphErr("", val_pull.size(), &val_pull[0], &points_nuis_pull[0], &down_pull[0], &up_pull[0], NULL, NULL);
    for ( auto itr : val_pull ) 
        std::cout << "@@@@ val_pull = " << itr << std::endl;
    for ( auto itr : points_nuis_pull ) 
        std::cout << "@@@@ points_nuis_pull = " << itr << std::endl;
    for ( auto itr : down_pull ) 
        std::cout << "@@@@ down_pull = " << itr << std::endl;
    for ( auto itr : up_pull ) 
        std::cout << "@@@@ up_pull = " << itr << std::endl;
    gr->SetLineColor(kBlack);
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(20);
    gr->SetLineStyle(1);
    gr->SetLineWidth(3);
    gr->SetMarkerSize(5);
    gr->GetXaxis()->SetTitleOffset(1.2);

    TGraphAsymmErrors* gr_norm = pull_helper.makeGraphErr("", val_norm.size(), &val_norm[0], &points_nuis_norm[0], &down_norm[0], &up_norm[0], NULL, NULL);
    gr_norm->SetLineColor(kRed);
    gr_norm->SetMarkerColor(kRed);
    gr_norm->SetMarkerStyle(24);
    gr_norm->SetLineStyle(1);
    gr_norm->SetLineWidth(3);
    gr_norm->SetMarkerSize(markerSize);
    gr_norm->GetXaxis()->SetTitleOffset(1.2);
    //
    h->LabelsOption("h");
    double labelSize = 1./nrNuis;
    h->SetLabelSize(labelSize>UserlabelSize?UserlabelSize:labelSize,"Y");//label size
    h->GetXaxis()->SetLabelColor(kWhite);
    h->GetXaxis()->SetAxisColor(kWhite);
    h->GetYaxis()->SetLabelColor(kBlack);
    h->GetYaxis()->SetAxisColor(kBlack);
    h->GetYaxis()->SetTickLength(0.);
    h->SetStats(0);
    h->Draw("h");
    // TODO: order should be the same for overlay, so just do it once

    // axis for the POI correlation
    TGaxis *axis_poi = new TGaxis(border_lo, nrNuis+1, border_hi, nrNuis+1, (-sigma_tot_lo) / scale_poi, (sigma_tot_hi) / scale_poi, 510, "-L");
    std::cout <<" Top range from "<< (-sigma_tot_lo) / scale_poi <<" to "<< (sigma_tot_hi) / scale_poi << std::endl;
    axis_poi->ImportAxisAttributes(h->GetXaxis());
    axis_poi->SetName("axis_poi");
    if (useRelativeImpact) axis_poi->SetTitle("#Delta#mu/#Delta#mu_{tot}");
    else                   axis_poi->SetTitle("#Delta#mu");
    axis_poi->SetTitleOffset(1.1);
    axis_poi->SetLineColor(kBlue);
    axis_poi->SetLabelColor(kBlue);
    axis_poi->SetTitleColor(kBlue);
    axis_poi->SetLabelSize(0.034);
    axis_poi->SetTitleSize(0.034);

    // axis for the nuisance parameter pull
    TGaxis *axis_theta = new TGaxis(border_lo, -offset, border_hi, -offset, (-sigma_tot_lo / max_poi) / scale_theta, (sigma_tot_hi / max_poi) / scale_theta, 510, "+R");
    std::cout << max_poi      << std::endl;
    std::cout << sigma_tot_hi << std::endl;
    std::cout << sigma_tot_lo << std::endl;
    std::cout << scale_theta  << std::endl;
    std::cout <<" Bottom range from "<< (-sigma_tot_lo / max_poi) / scale_theta <<" to "<< (sigma_tot_hi / max_poi) / scale_theta << std::endl;
    axis_theta->ImportAxisAttributes(h->GetXaxis());
    axis_theta->SetName("axis_theta");
    axis_theta->SetTitle("");
    axis_theta->SetTitleOffset(1.1);
    axis_theta->SetLineColor(kBlack);
    axis_theta->SetLabelColor(kBlack);
    axis_theta->SetTitleColor(kBlack);
    axis_theta->SetLabelSize(0.034);
    axis_theta->SetTitleSize(0.034);

    // axis for the nuisance parameter labels
    TGaxis *axis_label = new TGaxis(border_lo, 0, border_lo, nrNuis+1, 0, nrNuis+1, 0, "-R");
    axis_label->SetLineColor(kBlack);
    axis_label->SetTitleColor(kWhite);
    axis_label->SetLabelSize(0);
    axis_label->SetNdivisions(0);

    // some line definitions
    TLine l;
    l.SetLineWidth(2);
    l.SetLineColor(color_pulls);
    l.SetLineStyle(2);

    TLine l_stat;
    l_stat.SetLineWidth(2);
    l_stat.SetLineColor(color_staterror);
    l_stat.SetLineStyle(2);

    TLine l_syst;
    l_syst.SetLineWidth(2);
    l_syst.SetLineColor(color_systerror);
    l_syst.SetLineStyle(2);

    TLine l_tot;
    l_tot.SetLineWidth(2);
    l_tot.SetLineColor(color_totalerror);
    l_tot.SetLineStyle(2);

    // draw the nuisance parameter pulls including error bands and impact on poi
    gr_shades->Draw("p2");
    //gr_shades->GetXaxis()->SetRange(-2,2);

    if (drawPrefitImpactBand) {
        gr_poi_nom->Draw("p2");
    }
    if (drawPostfitImpactBand) {
        gr_poi->Draw("p2");
        gr_poi_pos->Draw("p2 same");
    }
    // draw axes
    if (drawPrefitImpactBand || drawPostfitImpactBand || drawErrorBars) axis_poi->Draw();
    axis_theta->Draw();
    axis_label->Draw();

    // draw +-1 and 0 sigma lines for pulls
    l.DrawLine( 0.              , 0.,  0.              , nrNuis);
    l.DrawLine( 1. * scale_theta, 0.,  1. * scale_theta, nrNuis);
    l.DrawLine(-1. * scale_theta, 0., -1. * scale_theta, nrNuis);
    gStyle->SetEndErrorSize(5.0);
    if(drawStandardBand){
        gr1s->Draw("p");
    }
    // draw syst and stat errors
    if (drawErrorBars) {
        l_stat.SetLineStyle(1);
        l_syst.SetLineStyle(1);
        l_tot.SetLineStyle(1);

        l_stat.DrawLine(-sigma_stat_lo * scale_poi / max_poi, 1.07*nrNuis,  sigma_stat_hi * scale_poi / max_poi, 1.07*nrNuis);
        l_stat.DrawLine( sigma_stat_hi * scale_poi / max_poi, 1.06*nrNuis,  sigma_stat_hi * scale_poi / max_poi, 1.08*nrNuis);
        l_stat.DrawLine(-sigma_stat_lo * scale_poi / max_poi, 1.06*nrNuis, -sigma_stat_lo * scale_poi / max_poi, 1.08*nrNuis);

        l_syst.DrawLine(-sigma_syst_lo * scale_poi / max_poi, 1.10*nrNuis,  sigma_syst_hi * scale_poi / max_poi, 1.10*nrNuis);
        l_syst.DrawLine( sigma_syst_hi * scale_poi / max_poi, 1.09*nrNuis,  sigma_syst_hi * scale_poi / max_poi, 1.11*nrNuis);
        l_syst.DrawLine(-sigma_syst_lo * scale_poi / max_poi, 1.09*nrNuis, -sigma_syst_lo * scale_poi / max_poi, 1.11*nrNuis);

        l_tot.DrawLine(-sigma_tot_lo * scale_poi / max_poi, 1.13*nrNuis,  sigma_tot_hi * scale_poi / max_poi, 1.13*nrNuis);
        l_tot.DrawLine( sigma_tot_hi * scale_poi / max_poi, 1.12*nrNuis,  sigma_tot_hi * scale_poi / max_poi, 1.14*nrNuis);
        l_tot.DrawLine(-sigma_tot_lo * scale_poi / max_poi, 1.12*nrNuis, -sigma_tot_lo * scale_poi / max_poi, 1.14*nrNuis);

        TLatex t_stat;
        TLatex t_syst;
        TLatex t_tot;

        t_stat.SetTextSize(0.03);
        t_stat.SetTextAlign(32);
        t_stat.SetTextColor(color_staterror);
        t_stat.DrawLatex((-sigma_stat_lo-0.025) * scale_poi / max_poi, 1.07*nrNuis, "statistics");

        t_syst.SetTextSize(0.03);
        t_syst.SetTextAlign(32);
        t_syst.SetTextColor(color_systerror);
        t_syst.DrawLatex((-sigma_syst_lo-0.025) * scale_poi / max_poi, 1.10*nrNuis, "systematics");

        t_tot.SetTextSize(0.03);
        t_tot.SetTextAlign(32);
        t_tot.SetTextColor(color_totalerror);
        t_tot.DrawLatex((-sigma_tot_lo-0.025) * scale_poi / max_poi, 1.13*nrNuis, "total");

        t_stat.Draw();
        t_syst.Draw();
        t_tot.Draw();
    }

    gr->Draw("p");
    if(drawRedNorm) 
        gr_norm->Draw("p same");

    pad1->SetTicks(0, 0);

    c1->SetLogy(0);

    TLegend* leg = makeLeg();
    leg->SetX1(channelPosX + 0.27);
    leg->SetY1(channelPosY-0.0775);
    leg->SetX2(channelPosX + 0.77);
    leg->SetY2(channelPosY+0.02);
    leg->SetTextSize(0.0225);

    leg->AddEntry(gr, "Pull: (#hat{#theta} - #theta_{0})/#Delta#theta","lp");
    if(drawRedNorm) leg->AddEntry(gr_norm, "Normalisation","lp");
    if(drawStandardBand){
        leg->AddEntry(gr1s, "1 standard deviation","l");
    }
    if (drawPostfitImpactBand) {
        leg->AddEntry(gr_poi_pos, "+1#sigma Postfit Impact on #mu", "f");//LOST YOUR HAT?
        leg->AddEntry(gr_poi, "-1#sigma Postfit Impact on #mu", "f");
    }

    leg->Draw();

    std::cout << "total unc = " << (fabs(sigma_tot_hi) + fabs(sigma_tot_lo)) / 2 << std::endl;
    std::cout << "sum of sq = " << sqrt(sum_poi2) << std::endl;

    pad1->cd();

    labelPosY = channelPosY-0.02;
    //ATLASLabel(labelPosX,labelPosY,"\nInternal",1);

    TLatex p;
    p.SetNDC();
    p.SetTextFont(42);
    p.DrawLatex(labelPosX,labelPosY-0.04,labelTxt.c_str());

    TLatex t3;
    t3.SetTextSize(0.025);
    t3.SetNDC();
    TString lumi = lumi5;
    if      ( parsed[0].find("mc16a_")  != std::string::npos) lumi = lumi2;
    else if ( parsed[0].find("mc16d_")  != std::string::npos) lumi = lumi3;
    else if ( parsed[0].find("mc16ad_") != std::string::npos) lumi = lumi4;
    else if ( parsed[0].find("mc16ade_")!= std::string::npos) lumi = lumi5;

    TString stex = "#sqrt{s} = " + energy2 + " TeV";//, #int Ldt = " + lumi2 + " fb^{-1}";
    TString lumitex2 = lumi + " fb^{-1}"; // Updated luminosity style
    t3.DrawLatex(channelPosX-0.020, channelPosY, stex.Data());
    t3.DrawLatex(channelPosX-0.020, channelPosY-0.04, lumitex2.Data());

    TLatex t2;
    t2.SetTextSize(0.025);
    t2.SetNDC();

    char const* tmp = getenv("ANALYSISTYPE");
    if ( tmp == NULL ) {
        std::cout << "Couldn't getenv ANALYSISTYPE" << std::endl;
    } else {
        std::string s( tmp );
    }

    t2.DrawLatex(channelPosX-0.020, channelPosY-0.075, ("m_{A}="+mass+" GeV").c_str());

    std::string saveName = "./data/output/" + cardName + "/pulls_" + POIname + "_";
    if (!m_postFitOrder) saveName += "prefit_";
    saveName += mass;
    std::cout << saveName << std::endl;
    
    c1->Draw();
    c1->SaveAs( Form("%s.pdf", saveName.c_str()) ) ;
    delete c1;
}

// before we are using draw_pulls_core, now we want to do the multiple mus, so we add a shell outside draw_pulls_core
void draw_pulls_v2(std::string mass = "125", std::string cardName = "LQ3Up500GeV", float scale_factor = 1.7, std::string overlayCard="", bool postFitOrder = true) 
{
    std::vector<std::string> parsed = parseString(cardName, ":");
    cardName = parsed[0];

    // TODO: check if following two lines are needed?
    computeFlags(cardName);
    applyOverlay(overlayCard , overlay , "");

    // Create ascii tables
    root2ascii( "./data/output/LQ3Up500GeV/root-files/pulls" );
    root2ascii( "./data/output/LQ3Up500GeV/root-files/breakdown_add" );
    std::cout << "Finished to create the ascii files" << std::endl;

    // find out how many POIs 
    TFile* infile = NULL;
    infile = new TFile( "PullPlot/output/test/root-files/breakdown_add/total.root");
    TH1D* hist = (TH1D*)infile->Get("total");
    int nBins = hist->GetNbinsX();
    if ( nBins%3 != 0) {
        cout << "Wired number of bins, it should contain central/up/down for each POI " <<endl;
        exit(1);
    }

    std::cout << "Start the poi loop"<< std::endl;
    int POItotal = nBins/3;
    for ( int POIpos = 0; POIpos < POItotal; POIpos++) {
        std::string POIname = hist->GetXaxis()->GetBinLabel(POIpos*3+1);
        std::cout << "Do the ranking plot for POI: "<< POIname << " as " << POIpos+1 << " (starting from 1) of " << POItotal << std::endl;

        // for each POI do the plot
        draw_pulls_core( mass, cardName, scale_factor, overlayCard, postFitOrder, POIname, POIpos, POItotal);
    }

    infile->Close();
}

// ____________________________________________________________________________|__________
// Return vector of strings from textfile
vector<string> getLabel(const char* fileName, int nrPars) {
    vector<string> tmp_labels;
    Int_t nlines = 0;
    ifstream idFile(fileName);
    while (1) {
        if (!idFile.good() || nlines > nrPars) break;
        string label;
        idFile >> label;
        tmp_labels.push_back(label);
        nlines++;
    }
    return tmp_labels;
}

TString translateNPname(TString internalName, bool isMVA) {    
    if(internalName == "ATLAS_norm_Wbb") return "W+HF normalisation";
    if(internalName == "ATLAS_norm_Wcl") return "W+cl normalisation";
    if(internalName == "ATLAS_norm_Zbb") return "Z+HF normalisation";
if(internalName == "ATLAS_norm_Zcl") return "Z+cl normalisation";
if(internalName == "ATLAS_norm_ttbar") return "t#bar{t} normalisation";
if(internalName == "ATLAS_norm_ttbar_L0") return "0-lepton t#bar{t} normalisation";
if(internalName == "ATLAS_norm_ttbar_L1") return "1-lepton t#bar{t} normalisation";
if(internalName == "ATLAS_norm_ttbar_L2") return "2-lepton t#bar{t} normalisation";
if(internalName == "alpha_ATLAS_LUMI_2012") return "Luminosity";
if(internalName == "alpha_ATLAS_LUMI_2015") return "Luminosity";
if(internalName == "alpha_ATLAS_LUMI_2015_2016") return "Luminosity";
if(internalName == "alpha_ATLAS_LUMI_2015_2017") return "Luminosity";
if(internalName == "alpha_ATLAS_LUMI_2015_2018") return "Luminosity";
//Run 2 b-tagging
if(internalName == "alpha_SysFT_EFF_Eigen_Light_0")return "Light-flavour tagging efficiency 0";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_1")return "Light-flavour tagging efficiency 1";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_2")return "Light-flavour tagging efficiency 2";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_3")return "Light-flavour tagging efficiency 3";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_4")return "Light-flavour tagging efficiency 4";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_5")return "Light-flavour tagging efficiency 5";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_6")return "Light-flavour tagging efficiency 6";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_7")return "Light-flavour tagging efficiency 7";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_8")return "Light-flavour tagging efficiency 8";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_9")return "Light-flavour tagging efficiency 9";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_10")return "Light-flavour tagging efficiency 10";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_11")return "Light-flavour tagging efficiency 11";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_12")return "Light-flavour tagging efficiency 12";
if(internalName == "alpha_SysFT_EFF_Eigen_Light_13")return "Light-flavour tagging efficiency 13";
if(internalName == "alpha_SysFT_EFF_Eigen_C_0")return "c-jet tagging efficiency 0";
if(internalName == "alpha_SysFT_EFF_Eigen_C_1")return "c-jet tagging efficiency 1";
if(internalName == "alpha_SysFT_EFF_Eigen_C_2")return "c-jet tagging efficiency 2";
if(internalName == "alpha_SysFT_EFF_Eigen_C_3")return "c-jet tagging efficiency 3";
if(internalName == "alpha_SysFT_EFF_Eigen_B_0")return "b-jet tagging efficiency 0";
if(internalName == "alpha_SysFT_EFF_Eigen_B_1")return "b-jet tagging efficiency 1";
if(internalName == "alpha_SysFT_EFF_Eigen_B_2")return "b-jet tagging efficiency 2";
if(internalName == "alpha_SysFT_EFF_Eigen_B_3")return "b-jet tagging efficiency 3";
if(internalName == "alpha_SysFT_EFF_Eigen_B_4")return "b-jet tagging efficiency 4";
if(internalName == "alpha_SysFT_EFF_Eigen_B_5")return "b-jet tagging efficiency 5";
if(internalName == "alpha_SysFT_EFF_extrapolation")return "b-jet tagging extrapolation";

// 19NP Scheme
if(internalName == "alpha_SysJET_19NP_JET_EffectiveNP_1")return"JES NP 1";
if(internalName == "alpha_SysJET_19NP_JET_EffectiveNP_2")return"JES NP 2";
if(internalName == "alpha_SysJET_19NP_JET_EffectiveNP_3")return"JES NP 3";
if(internalName == "alpha_SysJET_19NP_JET_EffectiveNP_4")return"JES NP 4";
if(internalName == "alpha_SysJET_19NP_JET_EffectiveNP_5")return"JES NP 5";
if(internalName == "alpha_SysJET_19NP_JET_EffectiveNP_6restTerm")return"JES NP 6";
if(internalName == "alpha_SysJET_19NP_JET_EtaIntercalibration_Modelling")return"JES #eta intercalibration modelling";
if(internalName == "alpha_SysJET_19NP_JET_EtaIntercalibration_TotalStat")return"JES #eta intercalibration stat.";
if(internalName == "alpha_SysJET_19NP_JET_EtaIntercalibration_NonClosure")return"JES #eta intercalibration non-closure";
if(internalName == "alpha_SysJET_19NP_JET_Pileup_OffsetMu")return"JES PU offset(#mu)";
if(internalName == "alpha_SysJET_19NP_JET_Pileup_OffsetNPV")return"JES PU offset(N_{PV})";
if(internalName == "alpha_SysJET_19NP_JET_Pileup_PtTerm")return"JES PU p_{T} term";
if(internalName == "alpha_SysJET_19NP_JET_Pileup_RhoTopology")return"JES PU #rho topology";
if(internalName == "alpha_SysJET_19NP_JET_BJES_Response")return"b-jets response";
if(internalName == "alpha_SysJET_19NP_JET_PunchThrough_MCTYPE")return"JES punch-through MC Type";
if(internalName == "alpha_SysJET_19NP_JET_SingleParticle_HighPt")return"JES single particle hi-p_{T}";
if(internalName == "alpha_SysJET_19NP_JET_Flavor_Response")return"JES flavour response";
if(internalName == "alpha_SysJET_19NP_JET_Flavor_Composition_Top")return"JES top flavour composition";
if(internalName == "alpha_SysJET_19NP_JET_Flavor_Composition_Zjets")return"JES Z+jets flavour composition";
if(internalName == "alpha_SysJET_19NP_JET_Flavor_Composition_Wjets")return"JES W+jets flavour composition";

if(internalName == "alpha_SysJET_JER_SINGLE_NP")return"Jet energy resolution";

// Correlate
if(internalName == "alpha_SysJET_19NP_JET_Flavor_Composition")return"JES flavour composition";

if(internalName == "alpha_SysBJetReso") return "b-jet energy resolution";
if(internalName == "alpha_SysFTBTagB0Effic_Y2012") return "b-jet tagging efficiency 0";
if(internalName == "alpha_SysBTagB0Effic_Y2012") return "b-jet tagging efficiency 0";
if(internalName == "alpha_SysBTagB1Effic_Y2012") return "b-jet tagging efficiency 1";
if(internalName == "alpha_SysBTagB2Effic_Y2012") return "b-jet tagging efficiency 2";
if(internalName == "alpha_SysBTagB3Effic_Y2012") return "b-jet tagging efficiency 3";
if(internalName == "alpha_SysBTagB4Effic_Y2012") return "b-jet tagging efficiency 4";
if(internalName == "alpha_SysBTagB5Effic_Y2012") return "b-jet tagging efficiency 5";
if(internalName == "alpha_SysBTagB6Effic_Y2012") return "b-jet tagging efficiency 6";
if(internalName == "alpha_SysBTagB7Effic_Y2012") return "b-jet tagging efficiency 7";
if(internalName == "alpha_SysBTagB8Effic_Y2012") return "b-jet tagging efficiency 8";
if(internalName == "alpha_SysBTagB9Effic_Y2012") return "b-jet tagging efficiency 9";
if(internalName == "alpha_SysBTagBPythia8_Y2012") return "#splitline{Generator dependence of}{    tagging efficiency}";
if(internalName == "alpha_SysBTagBSherpa_Y2012") return "Generator dependence of tagging eff.(for Sherpa)";
if(internalName == "alpha_SysBTagC0Effic_Y2012") return "c-jet tagging efficiency 0";
if(internalName == "alpha_SysBTagC1Effic_Y2012") return "c-jet tagging efficiency 1";
if(internalName == "alpha_SysBTagC2Effic_Y2012") return "c-jet tagging efficiency 2";
if(internalName == "alpha_SysBTagC3Effic_Y2012") return "c-jet tagging efficiency 3";
if(internalName == "alpha_SysBTagC4Effic_Y2012") return "c-jet tagging efficiency 4";
if(internalName == "alpha_SysBTagC5Effic_Y2012") return "c-jet tagging efficiency 5";
if(internalName == "alpha_SysBTagC6Effic_Y2012") return "c-jet tagging efficiency 6";
if(internalName == "alpha_SysBTagC7Effic_Y2012") return "c-jet tagging efficiency 7";
if(internalName == "alpha_SysBTagC8Effic_Y2012") return "c-jet tagging efficiency 8";
if(internalName == "alpha_SysBTagC9Effic_Y2012") return "c-jet tagging efficiency 9";
if(internalName == "alpha_SysBTagC10Effic_Y2012") return "c-jet tagging efficiency 10";
if(internalName == "alpha_SysBTagC11Effic_Y2012") return "c-jet tagging efficiency 11";
if(internalName == "alpha_SysBTagC12Effic_Y2012") return "c-jet tagging efficiency 12";
if(internalName == "alpha_SysBTagC13Effic_Y2012") return "c-jet tagging efficiency 13";
if(internalName == "alpha_SysBTagC14Effic_Y2012") return "c-jet tagging efficiency 14";
if(internalName == "alpha_SysBTagCPythia8_Y2012") return "Generator dependence of tagging eff. (for charm)";
if(internalName == "alpha_SysBTagCSherpa_Y2012") return "Generator dependence of tagging eff. (for Sherpa, charm)";
if(internalName == "alpha_SysBTagL0Effic_Y2012") return "l-jet tagging efficiency 0";
if(internalName == "alpha_SysBTagL1Effic_Y2012") return "l-jet tagging efficiency 1";
if(internalName == "alpha_SysBTagL2Effic_Y2012") return "l-jet tagging efficiency 2";
if(internalName == "alpha_SysBTagL3Effic_Y2012") return "l-jet tagging efficiency 3";
if(internalName == "alpha_SysBTagL4Effic_Y2012") return "l-jet tagging efficiency 4";
if(internalName == "alpha_SysBTagL5Effic_Y2012") return "l-jet tagging efficiency 5";
if(internalName == "alpha_SysBTagL6Effic_Y2012") return "l-jet tagging efficiency 6";
if(internalName == "alpha_SysBTagL7Effic_Y2012") return "l-jet tagging efficiency 7";
if(internalName == "alpha_SysBTagL8Effic_Y2012") return "l-jet tagging efficiency 8";
if(internalName == "alpha_SysBTagL9Effic_Y2012") return "l-jet tagging efficiency 9";
if(internalName == "alpha_SysElecE") return "Electron energy";
//if(internalName == "alpha_SysElecEffic") return "Electron efficiency";
if(internalName == "alpha_SysJVF_Y2012") return "jet vertex fraction";
if(internalName == "alpha_SysJetBE") return "b-jet energy";//
if(internalName == "alpha_SysJetEResol_Y2012") return "Jet energy resolution";
if(internalName == "alpha_SysJetEtaModel") return "Jet energy scale eta modeling";
//if(internalName == "alpha_SysJetEtaStat_Y2012")
if(internalName == "alpha_SysJetFlavB") return "Modeling of b-jet response";
//if(internalName == "alpha_SysJetFlavComp_Top")
//if(internalName == "alpha_SysJetFlavComp_VHVV")
//if(internalName == "alpha_SysJetFlavComp_Wjets")
//if(internalName == "alpha_SysJetFlavComp_Zjets")
//if(internalName == "alpha_SysJetFlavResp_Top")
//if(internalName == "alpha_SysJetFlavResp_VHVV")
//if(internalName == "alpha_SysJetFlavResp_Wjets")
//if(internalName == "alpha_SysJetFlavResp_Zjets")
//if(internalName == "alpha_SysJetMu")
if(internalName == "alpha_SysJetNP1_Y2012") return "Jet energy scale 1";
if(internalName == "alpha_SysJetNP2_Y2012") return "Jet energy scale 2";
if(internalName == "alpha_SysJetNP3_Y2012") return "Jet energy scale 3";
if(internalName == "alpha_SysJetNP4_Y2012") return "Jet energy scale 4";
if(internalName == "alpha_SysJetNP5_Y2012") return "Jet energy scale 5";
if(internalName == "alpha_SysJetNP6_rest_Y2012") return "Jet energy scale 6 and rest";
//if(internalName == "alpha_SysJetNPV")
if(internalName == "alpha_SysJetNonClos") return "Jet energy scale MC non-closure";
//if(internalName == "alpha_SysJetPilePt_Y2012")
//if(internalName == "alpha_SysJetPileRho_Y2012")
//if(internalName == "alpha_SysLepVeto")
if(internalName == "alpha_SysMETResoSoftTerms_Y2012") return "resolution E_{T}^{miss,SoftTerm}";
if(internalName == "alpha_SysMET_SoftTrk_Scale") return "resolution E_{T}^{miss,SoftTerm}";
if(internalName == "alpha_SysMET_SoftTrk_ResoPara") return "E_{T}^{miss} soft term resolution";// #parallel";
if(internalName == "alpha_SysMET_SoftTrk_ResoPerp") return "E_{T}^{miss} soft term resolution";// #perp";
if(internalName == "alpha_SysMETScaleSoftTerms_Y2012") return "scale E_{T}^{miss,SoftTerm}";
if(internalName == "alpha_SysMETTrigTop") return "E_{T}^{miss} top trigger";
if(internalName == "alpha_SysMJ_El_METstr") return "Multijet template";
// if(internalName == "alpha_SysMJ_El_flavor") return "Multijet flavor (electron)";
// if(internalName == "alpha_SysMUONS_MS") return "Muon MS";
// if(internalName == "alpha_SysEG_SCALE_ALL") return "Electron scale";
//if(internalName == "alpha_SysMJElCaloIso_T1")
//if(internalName == "alpha_SysMJElCaloIso_T2")
//if(internalName == "alpha_SysMJElDR")
//if(internalName == "alpha_SysMJElNorm_J2_T1")
//if(internalName == "alpha_SysMJElNorm_J2_TTypell")
//if(internalName == "alpha_SysMJElNorm_J2_TTypemm")
//if(internalName == "alpha_SysMJElNorm_J2_TTypett")
//if(internalName == "alpha_SysMJElNorm_J3_T1")
//if(internalName == "alpha_SysMJElNorm_J3_T2")
//if(internalName == "alpha_SysMJElPtV")
//if(internalName == "alpha_SysMJElTrkIso_T1_J2")
//if(internalName == "alpha_SysMJElTrkIso_T1_J3")
//if(internalName == "alpha_SysMJElTrkIso_T2_J2")
//if(internalName == "alpha_SysMJElTrkIso_T2_J3")
//if(internalName == "alpha_SysMJMuNorm_J2_T1")
//if(internalName == "alpha_SysMJMuNorm_J2_TTypell")
//if(internalName == "alpha_SysMJMuNorm_J2_TTypemm")
//if(internalName == "alpha_SysMJMuNorm_J2_TTypett")
//if(internalName == "alpha_SysMJMuNorm_J3_T1")
//if(internalName == "alpha_SysMJMuNorm_J3_T2")
//if(internalName == "alpha_SysMJMuTrkIso_T1_J2")
//if(internalName == "alpha_SysMJMuTrkIso_T1_J3")
//if(internalName == "alpha_SysMJMuTrkIso_T2_J2")
//if(internalName == "alpha_SysMJMuTrkIso_T2_J3")
//if(internalName == "alpha_SysMJ_J2_T1_L0_Y2012")
//if(internalName == "alpha_SysMJ_J2_T1_L0_Y2012_B1")
if(internalName == "alpha_SysMJNorm") return "Multijet normalisation)";
if(internalName == "alpha_SysMJ_J2_T2_L0_Y2012") return "Zerolepton Multijet (2-jet, p_{T}^{V} > 120 GeV)";
if(internalName == "alpha_SysMJ_J2_T2_L0_Y2012_B1") return "Zerolepton Multijet (2-jet, p_{T}^{V} < 120 GeV)";
//if(internalName == "alpha_SysMJ_J3_T1_L0_Y2012")
//if(internalName == "alpha_SysMJ_J3_T2_L0_Y2012")
if(internalName == "alpha_SysMJ_L2_Y2012") return "2-lepton multijet";
//if(internalName == "alpha_SysMuonEffic")
//if(internalName == "alpha_SysSChanAcerMC")
//if(internalName == "alpha_SysSChanAcerMCPS")
if(internalName == "alpha_SysTChanPtB2") return "Single top t-channel acceptance";
if(internalName == "alpha_SysTheoryAccPS") return "Signal acceptance (parton shower)";
if(internalName == "alpha_SysTheoryPDFAccPS") return "VH acc. for PS/UE tunes";
if(internalName == "alpha_SysstoptAcc") return "Single top t-ch. acc.";
if(internalName == "alpha_SysstopWtAcc") return "Single top Wt acc.";
//if(internalName == "alpha_SysTheoryAccPDF_ggZH")
//if(internalName == "alpha_SysTheoryAccPDF_qqVH")
//if(internalName == "alpha_SysTheoryAcc_J2_ggZH")
//if(internalName == "alpha_SysTheoryAcc_J2_qqVH")
//if(internalName == "alpha_SysTheoryAcc_J3_ggZH")
//if(internalName == "alpha_SysTheoryAcc_J3_qqVH")
//if(internalName == "alpha_SysTheoryBRbb")
//if(internalName == "alpha_SysTheoryPDF_ggZH")
//if(internalName == "alpha_SysTheoryPDF_qqVH")
if(internalName == "alpha_SysTheoryQCDscale_ggZH") return "QCD scale for ggZH";
if(internalName == "alpha_SysTheoryQCDscale_qqVH") return "QCD scale for qqVH";
//if(internalName == "alpha_SysTheoryVHPt")
if(internalName == "alpha_SysTheoryVPtQCD") return "#splitline{signal shape uncertainty}{    from QCD/PDF}";
if(internalName == "alpha_SysTheoryVPtQCD_ggZH") return "#splitline{ggZH shape uncertainty}{    from QCD/PDF}";
if(internalName == "alpha_SysTheoryVPtQCD_qqVH") return "#splitline{qqVH shape uncertainty}{    from QCD/PDF}";
if(internalName == "alpha_SysTopPt") return "t#bar{t} truth avg p_{T} modeling";
if(internalName == "alpha_SysTruthTagDR_Y2012") return "truth tagging DR modeling";
if(internalName == "alpha_SysTtbarMBBCont") return "t#bar{t} m_{jj} shape";
//if(internalName == "alpha_SysTtbarMetCont")
//if(internalName == "alpha_SysVVJetPDFAlphaPt")
//if(internalName == "alpha_SysVVJetScalePtST1")
//if(internalName == "alpha_SysVVJetScalePtST2")
//if(internalName == "alpha_SysVVMbb_WW")
//if(internalName == "alpha_SysVVMbb_WZ")
//if(internalName == "alpha_SysVVMbb_ZZ")
//if(internalName == "alpha_SysWDPhi_J2_Wcl")
//if(internalName == "alpha_SysWDPhi_J2_Whf")
//if(internalName == "alpha_SysWDPhi_J2_Wl")
//if(internalName == "alpha_SysWDPhi_J3_Wcl")
//if(internalName == "alpha_SysWDPhi_J3_Whf")
//if(internalName == "alpha_SysWDPhi_J3_Wl")
if(internalName == "alpha_SysWMbb_B0_WbbORcc"){
    if(isMVA)
        return "#splitline{W+b#bar{b}, W+c#bar{c} m_{jj} shape}{    (p_{T}^{V} < 120 GeV)}";
    else
        return "#splitline{W+b#bar{b}, W+c#bar{c} m_{jj} shape}{    (p_{T}^{V} < 90 GeV)}";
}
if(internalName == "alpha_SysWMbb_B1_WbbORcc") return "#splitline{W+b#bar{b}, W+c#bar{c} m_{jj} shape}{    (90 < p_{T}^{V} < 120 GeV)}";
if(internalName == "alpha_SysWMbb_WbbORcc")    return "#splitline{W+b#bar{b}, W+c#bar{c} m_{jj} shape}{    (p_{T}^{V} > 120 GeV)}";
if(internalName == "alpha_SysWMbb_WbcORbl") return "W+bc, W+bl m_{jj} shape";
//if(internalName == "alpha_SysWMbb_Wcl")
//if(internalName == "alpha_SysWMbb_Wl")
if(internalName == "alpha_SysWPtV_J2_Whf") return "W+HF p_{T}^{V} shape (2-jet) ";
if(internalName == "alpha_SysWPtV_J3_Whf") return "W+HF p_{T}^{V} shape (3-jet) ";
if(internalName == "alpha_SysWbcWbbRatio") return "W+bc to W+b#bar{b} normalisation";
if(internalName == "alpha_SysWblWbbRatio") return "#splitline{W+bl to W+b#bar{b} normalisation}{    (p_{T}^{V} > 120 GeV)}";
if(internalName == "alpha_SysWblWbbRatio_B0"){
    if(isMVA)
        return "#splitline{W+bl to W+b#bar{b} normalisation}{    (p_{T}^{V} < 120 GeV)}";
    else
        return "#splitline{W+bl to W+b#bar{b} normalisation}{    (p_{T}^{V} < 90 GeV)}";
}
if(internalName == "alpha_SysWblWbbRatio_B1") return "#splitline{W+bl to W+b#bar{b} normalisation}{    (90 GeV < p_{T}^{V} < 120 GeV)}";
//if(internalName == "alpha_SysWccWbbRatio")
//if(internalName == "alpha_SysWclNorm_J3")
//if(internalName == "alpha_SysWhfNorm_J3")
//if(internalName == "alpha_SysWlNorm")
//if(internalName == "alpha_SysWlNorm_J3")
//if(internalName == "alpha_SysWtChanAcerMC")
//if(internalName == "alpha_SysWtChanPythiaHerwig")
if(internalName == "alpha_SysZDPhi_J2_ZbORc") return "Z+b#bar{b}, Z+c#bar{c} d#phi shape (2-jet)";
if(internalName == "alpha_SysZDPhi_J2_Zl") return "Z+l d#phi shape (2-jet)";
if(internalName == "alpha_SysZDPhi_J3_ZbORc") return "Z+b#bar{b}, Z+c#bar{c} d#phi shape (3-jet)";
if(internalName == "alpha_SysZDPhi_J3_Zl") return "Z+l d#phi shape (3-jet)";
if(internalName == "alpha_SysZMbb_ZbORc") return "Z+b#bar{b}, Z+c#bar{c} m_{jj} shape";
if(internalName == "alpha_SysZMbb_Zl") return "Z+light m_{jj} shape";
if(internalName == "alpha_SysZMbb") return "Z+jets m_{bb} shape";
if(internalName == "alpha_SysWMbb") return "W+jets m_{bb} shape";
if(internalName == "alpha_SysVVMbbME") return "Diboson m_{bb} shape";
if(internalName == "alpha_SysTTbarMBB") return "t#bar{t} m_{bb} shape";
if(internalName == "alpha_SysTTbarMBB_L0") return "0-lepton t#bar{t} m_{bb} shape";
if(internalName == "alpha_SysTTbarMBB_L2") return "2-lepton t#bar{t} m_{bb} shape";
if(internalName == "alpha_SysStoptMBB") return "Single top t-ch. m_{bb} shape";
if(internalName == "alpha_SysStopWtMBB") return "Single top Wt m_{bb} shape";
if(internalName == "alpha_SysStopWtothACC") return "Single top acceptance (Wt#rightarrowother)";
if(internalName == "alpha_SysZPtV_ZbORc") return "Z+b#bar{b}, Z+c#bar{c} p_{T}^{V}";
if(internalName == "alpha_SysZPtV_Zl") return "Z+light p_{T}^{V}";
if(internalName == "alpha_SysZPtV") return "Z+jets p_{T}^{V}";
if(internalName == "alpha_SysWPtV") return "W+jets p_{T}^{V}";
if(internalName == "alpha_SysTTbarPTV") return "t#bar{t} p_{T}^{V}";
if(internalName == "alpha_SysTTbarPTV_L2") return "2-lepton t#bar{t} p_{T}^{V}";
if(internalName == "alpha_SysStoptPTV") return "Single top t-ch. p_{T}^{V}";
if(internalName == "alpha_SysStopWtPTV") return "Single top Wt p_{T}^{V}";
if(internalName == "alpha_SysZbbNorm_L0") return "0-lepton Z+HF normalisation";//b#bar{b}
if(internalName == "alpha_SysZbbNorm_L2") return "2-lepton Z+HF normalisation";//b#bar{b}
if(internalName == "alpha_SysZbbNorm") return "Z+HF normalisation";//b#bar{b}
if(internalName == "alpha_SysWbbNorm_L1") return "1-lepton W+HF normalisation";//b#bar{b}
if(internalName == "alpha_SysWbbNorm") return "W+HF normalisation";//b#bar{b}
if(internalName == "alpha_SysZZNorm_L0") return "0-lepton ZZ normalisation";
if(internalName == "alpha_SysZZNorm_L2") return "2-lepton ZZ normalisation";
if(internalName == "alpha_SysZZNorm") return "ZZ normalisation";
if(internalName == "alpha_SysWZNorm_L0") return "0-lepton WZ normalisation";
if(internalName == "alpha_SysWZNorm_L2") return "2-lepton WZ normalisation";
if(internalName == "alpha_SysWZNorm") return "WZ normalisation";
if(internalName == "alpha_SysVZNorm_J2") return "VZ normalisation (2-jet)";
if(internalName == "alpha_SysVZNorm_J3") return "VZ normalisation (3-jet)";
TString anaType = std::getenv("ANALYSISTYPE");
if(anaType == "boostedVHbbRun2"){
    if(internalName == "alpha_SysVZNorm") return "VZ acceptance";
}else{
    if(internalName == "alpha_SysVZNorm") return "VZ normalisation";
}

if(internalName == "alpha_SysZbbNorm_J2") return "Z+HF normalisation (2-jet)";//b#bar{b}
if(internalName == "alpha_SysZbbNorm_J3") return "Z+HF normalisation (3-jet)";//b#bar{b}
if(internalName == "alpha_SysWbbNorm_J2") return "W+HF normalisation (2-jet)";//b#bar{b}
if(internalName == "alpha_SysWbbNorm_J3") return "W+HF normalisation (3-jet)";//b#bar{b}
if(internalName == "alpha_SysWbbNorm") return "W+HF normalisation";//b#bar{b}
if(internalName == "alpha_SysZbcZbbRatio") return "Z+bc to Z+b#bar{b} normalisation";
if(internalName == "alpha_SysZblZbbRatio") return "Z+bl to Z+HF normalisation";
if(internalName == "alpha_SysZblZbbRatio_J2") return "Z+bl to Z+b#bar{b} normalisation (2-jet)";
if(internalName == "alpha_SysZblZbbRatio_J3") return "Z+bl to Z+b#bar{b} normalisation (3-jet)";
if(internalName == "alpha_SysZccZbbRatio") return "Z+bc to Z+b#bar{b} normalisation";
if(internalName == "alpha_SysZclNorm") return "Z+cl normalisation";
//if(internalName == "alpha_SysZlNorm")
//if(internalName == "alpha_SysZlNorm_J3")
if(internalName == "alpha_SysstopWtNorm") return "Single top Wt normalisation";
//if(internalName == "alpha_SysstopsNorm")
//if(internalName == "alpha_SysstoptNorm")
if(internalName == "alpha_SysTheoryUEPSAcc") return "VH acceptance (PS/UE)";
if(internalName == "alpha_SysTheoryUEPSAcc_J3") return "VH acceptance in 3-jet (PS/UE)";
if(internalName == "alpha_SysTheoryAcc_J2_qqVH") return "VH acceptance (QCD scales)";
if(internalName == "ATLAS_norm_Zbb_J3") return "Z + HF 3-jet scale factor";
if(internalName == "alpha_SysTheoryAcc_J3_qqVH") return "VH 2-to-3 jets acc. (QCD)";
if(internalName.Contains("1NP_JET_Flavor_Composition_VV")) return "VV JES flav. comp. uncert.";
if(internalName == "alpha_SysttbarHighPtV") return "t#bar{t} high p_{T}^{V} normalisation";
if(internalName == "alpha_SysttbarNorm_J3") return "t#bar{t} normalisation (3-jet)";
if(internalName == "alpha_SysttbarNorm_J3_L2") return "t#bar{t} normalisation (3-jet, dilepton)";
if(internalName == "alpha_SysMJ_L2_Y2012_Spctopemucr") return "2-lepton Multijet (top C.R.)";
if(internalName == "gamma_stat_Region_BMin150_J2_T2_L1_Y2015_distmva_DWhfSR_bin_13") return "MC stat. 1 lepton 2 jet, bin 13";
if(internalName == "gamma_stat_Region_BMin150_J2_T2_L1_Y2015_distmva_DWhfSR_bin_12") return "MC stat. 1 lepton 2 jet, bin 12";
//boosted VHbb renaming
if( (internalName == "alpha_SysFATJET_JMR_VH") || (internalName == "alpha_SysMJJMR_Higgs")) return "large-R jet JMR VH";
if( (internalName == "alpha_SysFATJET_JMR_VV") || (internalName == "alpha_SysMJJMR_Diboson")) return "large-R jet JMR VV";
if( (internalName == "alpha_SysFATJET_JMR_Vjets") ) return "large-R jet JMR V+jets";
if(internalName == "alpha_SysFATJET_JMR_Rest") return "large-R jet JMR V+jets, single top";
if(internalName == "alpha_SysFATJET_JMR_Vjets") return "large-R jet JMR V+jets";
if(internalName == "alpha_SysFATJET_Medium_JET_Comb_Modelling_Kin_Vjets") return "large-R JMS/JES Modelling V+jets"; 
if(internalName == "alpha_SysFATJET_Medium_JET_Comb_Modelling_Kin_Diboson") return "large-R JMS/JES Modelling VV"; 
if(internalName == "alpha_SysFATJET_Medium_JET_Comb_Modelling_Kin_Rest") return "large-R JMS/JES Modelling V+jets, single top";
if(internalName == "alpha_SysFATJET_Medium_JET_Comb_Modelling_Kin_VV") return "large-R JMS/JES Modelling VH, VV";
if(internalName == "alpha_SysFATJET_Medium_JET_Comb_Modelling_Kin_Diboson") return "large-R JMS/JES Modelling VV";
if(internalName == "alpha_SysFATJET_Medium_JET_Comb_Baseline_Kin_Vjets") return "large-R JMS/JES Baseline V+jets";
if(internalName == "alpha_SysFATJET_Medium_JET_Comb_Baseline_Kin_Diboson") return "large-R JMS/JES Baseline VV"; 
if(internalName == "alpha_SysFATJET_Medium_JET_Comb_Baseline_Kin_VV") return "large-R JMS/JES Baseline VH, VV";
if(internalName == "alpha_SysFATJET_Medium_JET_Comb_Baseline_Kin_Rest") return "large-R JMS/JES Baseline V+jets, single top";
if(internalName == "alpha_SysFATJET_Medium_JET_Comb_Baseline_Kin_Diboson") return "large-R JMS/JES Baseline VV";
if(internalName == "alpha_SysTTbarMJPS") return "t#bar{t} m_{J} shape PS";
if(internalName == "alpha_SysZMbbBoosted") return "Z+jets m_{J} shape #mu_{R}";
if(internalName == "alpha_SysWMbbBoosted") return "W+jets m_{J} shape #mu_{R}";
if(internalName == "alpha_SysStopWtMJ") return "Single top Wt m_{J} shape DS/DR";
if(internalName == "alpha_SysWbbNorm_L0") return "0- to 1-lepton W+HF ratio";
if(internalName == "alpha_SysVVPP8MJBoosted_Diboson") return "VZ m_{J} shape PwPy";
if(internalName == "alpha_SysRatioHPLP_DSRnoaddbjetsr_J0_Whf") return "W+HF HP/LP ratio";
if(internalName == "alpha_SysRatioHPLP_DSRnoaddbjetsr_J0_Stop") return "single top HP/LP ratio";
if(internalName == "alpha_SysRatioHPLP_DSRnoaddbjetsr_L1_J0_Ttbar") return "1-lepton t#bar{t} HP/LP ratio";
if(internalName == "alpha_SysRatioHPLP_VZ_DSRnoaddbjetsr_J0") return "VZ HP/LP ratio";
if(internalName == "alpha_SysRatioMedHighPtv_BMin400_Stop") return "Single top med/high p_{T}^{V} ratio";
if(internalName == "alpha_SysRatioMedHighPtv_VZ_BMin400") return "VZ med/high p_{T}^{V} ratio";
if(internalName == "alpha_SysRatioSRCR_DSRnoaddbjetsr_L1_Ttbar") return "1-lepton t#bar{t} SR/CR ratio";
if(internalName == "SigXsecOverSM_VZwithVH") return "Signal strength #mu_{VZ}";
if(internalName == "SigXsecOverSM") return "Signal strength #mu_{VH}";
if(internalName == "alpha_SysVVPP8MJBoosted") return "VV m_{J} shape PP8";
if(internalName == "alpha_SysVVPP8MJBoosted_Diboson") return "VV m_{J} shape PP8";
if(internalName == "alpha_SysJET_CR_JET_JER_EffectiveNP_1") return "Effective small-R JER NP 1";
if(internalName == "alpha_SysRatiottbarVRjets_ttbar2VR") return "t#bar{t} 2/3+ assoc. VR-jet ratio";
if(internalName == "alpha_QCDScaleDelta75_qqVH") return "QCD Scale #Delta_{75}^{qq}";
if(internalName == "alpha_SysTheoryDelta1_qqVH") return "QCD Scale #Delta_{1}^{qq}";
if(internalName == "gamma_stat_Region_BMax400_BMin250_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L1_distmBB_J0_bin_6") return "MC stat. 1L med. p_{T}^{V} HPSR, bin 6";
if(internalName == "gamma_stat_Region_BMax400_BMin250_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L1_distmBB_J0_bin_7") return "MC stat. 1L med. p_{T}^{V} HPSR, bin 7";
if(internalName == "gamma_stat_Region_BMax400_BMin250_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L1_distmBB_J0_bin_3") return "MC stat. 1L med. p_{T}^{V} HPSR, bin 3";
if(internalName == "gamma_stat_Region_BMax400_BMin250_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L1_distmBB_J0_bin_0") return "MC stat. 1L med. p_{T}^{V} HPSR, bin 0";
if(internalName == "gamma_stat_Region_BMax400_BMin250_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L0_distmBB_J0_bin_6") return "MC stat. 0L med. p_{T}^{V} HPSR, bin 6";
if(internalName == "gamma_stat_Region_BMax400_BMin250_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L0_distmBB_J0_bin_7") return "MC stat. 0L med. p_{T}^{V} HPSR, bin 7";
if(internalName == "gamma_stat_Region_BMax400_BMin250_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L0_distmBB_J0_bin_3") return "MC stat. 0L med. p_{T}^{V} HPSR, bin 3";
if(internalName == "gamma_stat_Region_BMax400_BMin250_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L0_distmBB_J0_bin_0") return "MC stat. 0L med. p_{T}^{V} HPSR, bin 0";

if(internalName == "gamma_stat_Region_BMin400_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L0_distmBB_J0_bin_6") return "MC stat. 0L high p_{T}^{V} HPSR, bin 6";
if(internalName == "gamma_stat_Region_BMin400_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L1_distmBB_J0_bin_6") return "MC stat. 1L high p_{T}^{V} HPSR, bin 6";
if(internalName == "gamma_stat_Region_BMin400_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L0_distmBB_J0_bin_7") return "MC stat. 0L high p_{T}^{V} HPSR, bin 7";
if(internalName == "gamma_stat_Region_BMin400_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L1_distmBB_J0_bin_7") return "MC stat. 1L high p_{T}^{V} HPSR, bin 7";
if(internalName == "gamma_stat_Region_BMin400_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L0_distmBB_J0_bin_0") return "MC stat. 0L high p_{T}^{V} HPSR, bin 0";
if(internalName == "gamma_stat_Region_BMin400_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L0_distmBB_J0_bin_3") return "MC stat. 0L high p_{T}^{V} HPSR, bin 3";
if(internalName == "gamma_stat_Region_distmBB_J1_L0_T2_DSRnoaddbjetsr_Y6051_incJet1_Fat1_incFat1_BMin250_BMax400_bin_0") return "MC stat. 0L med. p_{T}^{V} LPSR, bin 0";
if(internalName == "gamma_stat_Region_distmBB_J1_L0_T2_DSRnoaddbjetsr_Y6051_incJet1_Fat1_incFat1_BMin250_BMax400_bin_3") return "MC stat. 0L med. p_{T}^{V} LPSR, bin 3";
if(internalName == "gamma_stat_Region_BMax400_BMin250_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L1_distmBB_J0_bin_6") return "MC stat. 1L med. p_{T}^{V} HPSR, bin 6";
if(internalName == "gamma_stat_Region_BMin400_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L1_distmBB_J0_bin_3") return "MC stat. 1L high p_{T}^{V} HPSR, bin 3";

if(internalName == "gamma_stat_Region_BMin400_incFat1_Fat1_Y6051_DSRnoaddbjetsr_T2_L1_distmBB_J0_bin_1") return "MC stat. 1L high p_{T}^{V} HPSR, bin 0";
if(internalName == "gamma_stat_Region_distmBB_J1_L1_T2_DSRnoaddbjetsr_Y6051_incJet1_Fat1_incFat1_BMin250_BMax400_bin_3") return "MC stat. 1L med. p_{T}^{V} LPSR, bin 3";



if(string(internalName.Data()).find("distmBB_J0")!=string::npos) return translateGammaStatName(internalName);

return internalName;
}

TString translateGammaStatName(TString internalName)
{   

    TPRegexp pat1("gamma_stat_Region_distmBB_J0_L([0-9])_T2_DSR([a-zA-Z]*)_Y[0-9]+_incJet1_Fat1_incFat1_BMin250_BMax400_bin_([0-9]+)");
    TPRegexp pat2("gamma_stat_Region_BMin400_incFat1_Fat1_incJet1_Y[0-9]+_DSR([a-zA-Z]*)_T2_L([0-9])_distmBB_J0_bin_([0-9]+)");

    TString lep;
    TString rg;    
    TString bin;
    TString ptv;

    TObjArray * res = pat1.MatchS(internalName);
    if(res->GetEntries()!=0)
    {
        lep = ((TObjString*)res->At(1))->GetString();
        rg = ((TObjString*)res->At(2))->GetString();
        bin = ((TObjString*)res->At(3))->GetString();
        ptv = "med. p_{T}^{V}";
    }
    else
    {
        res = pat2.MatchS(internalName);
        if(res->GetEntries()!=0)
        {
            rg = ((TObjString*)res->At(1))->GetString();
            lep = ((TObjString*)res->At(2))->GetString();
            bin = ((TObjString*)res->At(3))->GetString();
            ptv = "high p_{T}^{V}";
        }   
        else
        {
            return internalName;
        }
    }    
    TString reg = (string(rg.Data()).find("topaddbjetcr")!=string::npos) ? "topCR" : "SR";

    TString name = TString("MC stat. ") + lep + "L " + ptv + " " + reg + ", bin " + bin;    

    return name;
}



bool isSignalNP (string NPname) {

    if(removeSigUnc){
        if(isVH){ 
            if( (NPname.find("Theory")  != std::string::npos)    || 
                    (NPname.find("SysVH")   != std::string::npos)    || 
                    (NPname.find("ggVH")    != std::string::npos)    || 
                    (NPname.find("ggZH")    != std::string::npos)    || 
                    (NPname.find("qqVH")    != std::string::npos))
                return true;
        }else{
            if( (NPname.find("VZQCD")   != std::string::npos)    || 
                    (NPname.find("ZZQCD")   != std::string::npos)    || 
                    (NPname.find("ZZUEPS")  != std::string::npos)    || 
                    (NPname.find("SysVV")   != std::string::npos)    || 
                    (NPname.find("SysVZ")   != std::string::npos))
                return true;
        }
    }

    return false;

}
