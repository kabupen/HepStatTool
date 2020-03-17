
void root2ascii( TString dir_root )
{
    // output/test/root-files/pulls
    // output/test/root-files/breakdown_add
    TString dir_ascii = dir_root;
    dir_ascii.ReplaceAll("root-files", "ascii");
    std::system(("mkdir -vp "+dir_ascii));

    enum ConvertMode {pulls, breakdown};

    ConvertMode mode;
    if      ( dir_ascii.Contains("pulls") )    { mode = pulls; } 
    else if ( dir_ascii.Contains("breakdown")) { mode = pulls; } 
    else {
        std::cout << "ERROR::Something went wrong." << std::endl;
        exit(1);
    }

    std::vector<TString> list;
    TSystemDirectory  dire(dir_root, dir_root);
    TList* files = dire.GetListOfFiles();
    TIter next(files);
    while( TSystemFile* file = (TSystemFile*)next() ) {
        TString fname = file->GetName();
        if(file->IsDirectory()) continue;
        list.push_back(fname.Data());
    }

    std::stringstream outFileName;
    outFileName << dir_ascii;

    std::ofstream outFile      ((outFileName.str() + ".txt")      .c_str());
    std::ofstream outFile_id   ((outFileName.str() + "_id.txt")   .c_str());
    std::ofstream outFile_nf   ((outFileName.str() + "_nf.txt")   .c_str());
    std::ofstream outFile_nf_id((outFileName.str() + "_nf_id.txt").c_str());

    int nrNPs = 0;
    int nrNFs = 0;

    for( auto itr = list.begin(); itr != list.end(); itr++ ) {
        std::stringstream fileName;
        fileName << dir_root << "/" << *itr;

        bool isNorm   = false;
        bool isInfNan = false;

        TFile* infile = NULL;
        infile = new TFile(fileName.str().c_str());
        if ( infile && !infile->IsOpen() ) {
            cout << "ERROR: Could not open file " << fileName.str() << endl;
            continue;
        }

        TString histoName = itr->ReplaceAll(".root", "");

        TH1D* hist = (TH1D*)infile->Get(histoName);

        int nrBins = hist->GetNbinsX();

        /* Check the input file validity */
        for (int bin = 1; bin <= nrBins; bin++) {
            double number = hist->GetBinContent(bin);

            if ( number > 10e9 )   isInfNan = true; // check inf
            if ( number != number) isInfNan = true; // check nan
        }
        if ( isInfNan ) {
            cout << "WARNING::Skipping " << *itr << " because of inf/nan" << endl;
            continue;
        }

        ( isNorm ? outFile_nf : outFile ) << ( isNorm ? nrNPs++ : nrNFs++ ) << " ";

        for (int bin = 1; bin <= nrBins; bin++) {
            double number = hist->GetBinContent(bin);
            ( isNorm ? outFile_nf : outFile ) << number << ((bin < nrBins)?" ":"\n");
        }

        if   ( isNorm ) outFile_nf_id << * itr << "\n";
        else            outFile_id    << * itr << "\n";
        if( mode==breakdown ) outFile_id << "~~*~*~*~*~*~*~*~*~*~ " << * itr << "\n";

        infile->Close();
    }

    std::cout << "Writing to file: " << outFileName.str() << "*.txt" << std::endl;

    outFile      .close();
    outFile_id   .close();
    outFile_nf   .close();
    outFile_nf_id.close();
}
