#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <set>

class FileHolder
{

    public:
        FileHolder();
        FileHolder(std::string, int);
        ~FileHolder();

    public:
        void skipFirst() {skipFirstLine=true;}
        bool readInFile(std::string name);
        void writeOutFile(std::string name);

        void setUseStrings(bool flag=true) {useStrings=flag;}
        void setNrCols(int n) {nrCols = n;}
        void setRate(double mass, int col, double val);

        double getRate(int row, int col) { return rates[row][col];}
        double getRateByMass(double mass, int col);
        std::vector<double> getCol(int col);
        std::vector<double> getMassPoints() { return massPoints;};

        void copy(FileHolder& alt);
        void removeMass(double mass);
        void addMass(double mass, std::vector<double>& rate);
        void addFile(FileHolder& other);
        void addCol(std::vector<double>& addRates, int col);
        void useMedian(int col);

        int  nrCols;
        bool skipFirstLine;
        bool useStrings;
        std::vector<std::vector<double> > rates; // row: mass, col: sample
        std::vector<double>               massPoints;
        std::vector<std::string>          massPointsS;
        std::string                       _name;
        int m_min_mass = 90;
        int m_max_mass = 600;
};

FileHolder::FileHolder() :
    nrCols(8),
    skipFirstLine(false),
    useStrings(false)
{}

FileHolder::FileHolder(std::string file, int nrcols) :
        skipFirstLine(false),
        useStrings(false)
{
    FileHolder::setNrCols(nrcols);
    FileHolder::readInFile(file);
    std::vector<double> vec_massPoints_tmp = FileHolder::getMassPoints();
    int nrNumbers = vec_massPoints_tmp.size();
    for (int i=0; i < nrNumbers; i++) {
        if ( vec_massPoints_tmp[i] < m_min_mass || vec_massPoints_tmp[i] > m_max_mass) {
            FileHolder::removeMass(vec_massPoints_tmp[i]);
        }
    }
}

bool FileHolder::readInFile(std::string name )
{
    std::cout << "Reading in file: " << name << endl;
    ifstream inFile(name.c_str());
    _name=name;
    if (inFile.fail()) {
        std::cout << "ERROR::Couldn't open file: " << name << endl;
        return false;
    }

    if (skipFirstLine) {
        std::string junk;
        getline(inFile, junk);
    }

    int nrItr = 0;
    while (!inFile.eof()) {
        double mass;
        std::string massS;

        if (useStrings) {
            inFile >> massS;
        }
        else {
            inFile >> mass;
        }

        if (inFile.eof()) break;

        bool fill=false;
        if (!massPoints.size()/*bug in inputs*/) fill=true;
        else if (mass != massPoints.back()) fill=true;

        std::vector<double> numbers;
        double number;
        for (int i=0;i<nrCols;i++) {
            inFile >> number;
            numbers.push_back(number);
        }

        if (fill) {
            massPoints.push_back(mass);
            rates.push_back(numbers);
        }

        nrItr++;
        //std::cout << "mass = " << mass << endl;
        if (nrItr > 500) {
            std::cout << "ERROR::Line overflow detected. Exiting." << endl;
            return false;
        }
    }
    std::cout << "Done" << endl;
    inFile.close();
    return true;
}

void FileHolder::writeOutFile(std::string name)
{
    std::cout << "Writing file: " << name << endl;
    ofstream outFile(name.c_str());
    if (outFile.fail())
    {
        std::cout << "Error writing to file: " << name << endl;
        return;
    }

    int nrPoints = massPoints.size();
    for (int i=0;i<nrPoints;i++)
    {
        //std::cout << "adding point: " << massPoints[i] << endl;
        outFile << massPoints[i];
        for (int j=0;j<nrCols;j++)
        {
            //std::cout << "->rate=" << rates[i][j] << endl;
            outFile << " " << rates[i][j];
        }
        outFile << "\n";
    }
    outFile.close();
}

std::vector<double> FileHolder::getCol(int col)
{
    std::vector<double> vec_col;
    int nrPoints = massPoints.size();
    for (int i=0;i<nrPoints;i++)
    {
        vec_col.push_back(rates[i][col]);
    }
    return vec_col;
}

double FileHolder::getRateByMass(double mass, int col)
{
    int nrPoints = massPoints.size();
    for (int im=0;im<nrPoints;im++)
    {
        if (mass == massPoints[im])
        {
            return rates[im][col];
        }
    }
    return 0.;
}

void FileHolder::setRate(double mass, int col, double val)
{
    for (int imass=0;imass<(int)massPoints.size();imass++)
    {
        if (mass == massPoints[imass])
        {
            rates[imass][col] = val;
            break;
        }
    }
}

void FileHolder::copy(FileHolder& alt)
{
    alt.nrCols=nrCols;
    alt.rates=rates;
    alt.massPoints=massPoints;
}

void FileHolder::removeMass(double mass)
{
    std::vector<double> newPoints;
    std::vector<std::vector<double> > newRates;
    int nrPoints = massPoints.size();
    for (int i=0;i<nrPoints;i++)
    {
        if (massPoints[i] == mass) continue;
        newPoints.push_back(massPoints[i]);
        newRates.push_back(rates[i]);
    }
    massPoints = newPoints;
    rates = newRates;
}

void FileHolder::addMass(double mass, std::vector<double>& rate)
{
    //std::cout << "adding mass: " << mass << endl;
    std::vector<double> newPoints;
    std::vector<std::vector<double> > newRates;

    bool found = false;
    int nrPoints = massPoints.size();
    for (int i=0;i<nrPoints;i++)
    {
        if ((massPoints[i] < mass && !found) || (massPoints[i] > mass && found))
        {
            newPoints.push_back(massPoints[i]);
            newRates.push_back(rates[i]);
        }
        else if (massPoints[i] > mass && !found) 
        {
            //std::cout << "found" << endl;

            newPoints.push_back(mass);
            newRates.push_back(rate);

            newPoints.push_back(massPoints[i]);
            newRates.push_back(rates[i]);

            found = true;
        }
    }
    massPoints = newPoints;
    rates = newRates;
}

void FileHolder::addFile(FileHolder& other)
{
    std::vector<double> otherPoints = other.massPoints;
    std::vector<std::vector<double> > otherRates = other.rates;
    for (int i=0;i<(int)otherPoints.size();i++)
    {
        addMass(otherPoints[i], otherRates[i]);
    }
}

void FileHolder::addCol(std::vector<double>& addRates, int col)
{
    int nrPoints = massPoints.size();
    for (int i=0;i<nrPoints;i++)
    {
        std::vector<double> newRates;
        for (int j=0;j<col;j++)
        {
            newRates.push_back(rates[i][j]);
        }
        newRates.push_back(addRates[i]);
        for (int j=col;j<(int)rates[i].size();j++)
        {
            newRates.push_back(rates[i][j]);
        }
        rates[i]=newRates;
    }
    nrCols++;
}

void FileHolder::useMedian(int col)
{
    //std::cout << "in use median" << endl;
    set<double> vals;
    int nrRates = rates.size();
    for (int i=0;i<nrRates;i++)
    {
        vals.insert(rates[i][col]);
    }
    set<double>::iterator itr=vals.begin();
    int nrVals = vals.size();
    for (int i=0;i<nrVals/2;i++)
    {
        itr++;
    }

    for (int i=0;i<nrRates;i++)
    {
        rates[i][col] = *itr;
    }
    //std::cout << "out use median" << endl;
}
