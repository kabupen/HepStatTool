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
        bool readInFile(std::string name);
        void writeOutFile(std::string name);

        void setNrCols(int n) {m_nrCols = n;}
        void setRate(double mass, int col, double val);

        double getRate(int row, int col) { return m_rates[row][col];}
        double getRateByMass(double mass, int col);
        std::vector<double> getCol(int col);
        std::vector<double> getMassPoints() { return m_massPoints; }
        double getMassPointsSize() { return m_massPoints.size(); }

        void copy(FileHolder& alt);
        void removeMass(double mass);
        void addMass(double mass, std::vector<double>& rate);
        void addFile(FileHolder& other);
        void addCol(std::vector<double>& addm_rates, int col);
        void useMedian(int col);

        int  m_nrCols;
        std::vector<std::vector<double> > m_rates; // row: mass, col: sample
        std::vector<double>               m_massPoints;
        int m_min_mass = 90;
        int m_max_mass = 600;
};

FileHolder::FileHolder() :
    m_nrCols(8)
{}

FileHolder::FileHolder(std::string file, int columns_num) 
{
    if ( !FileHolder::readInFile(file) ) 
        exit(1);

    FileHolder::setNrCols(columns_num);
    std::vector<double> tmpvec = FileHolder::getMassPoints();
    int nrNumbers = tmpvec.size();
    for (int i=0; i < nrNumbers; i++) {
        //        if ( tmpvec.at(i) < m_min_mass || tmpvec.at(i) > m_max_mass) {
        //            FileHolder::removeMass(tmpvec[i]);
        //        }
    }
}

FileHolder::~FileHolder()
{
}

bool FileHolder::readInFile( std::string name )
{
    std::cout << "Reading in file: " << name << std::endl;
    std::ifstream input_file(name.c_str());
    if ( input_file.fail() ) {
        std::cout << "ERROR::Couldn't open file: " << name << std::endl;
        return false;
    }

    if ( input_file.eof() ) return true;

    std::string line;
    std::vector<double> tmp;
    for ( int iRow = 0; getline(input_file, line); iRow++ ) {
        m_massPoints.push_back(iRow);

        std::istringstream stream(line);
        double data;
        for ( int col = 0; stream >> data; col++ ) { 
            if ( col == 0 ) continue;
            std::cout << "\t" << data << std::endl;
            tmp.push_back(data);
        }
        m_rates.push_back(tmp);
    }
    input_file.close();
    return true;
}

void FileHolder::writeOutFile(std::string name)
{
    std::cout << "Writing file: " << name << std::endl;
    ofstream outFile(name.c_str());
    if (outFile.fail())
    {
        std::cout << "Error writing to file: " << name << std::endl;
        return;
    }

    int nrPoints = m_massPoints.size();
    for (int i=0;i<nrPoints;i++)
    {
        //std::cout << "adding point: " << m_massPoints[i] << std::endl;
        outFile << m_massPoints[i];
        for (int j=0;j<m_nrCols;j++)
        {
            //std::cout << "->rate=" << m_rates[i][j] << std::endl;
            outFile << " " << m_rates[i][j];
        }
        outFile << "\n";
    }
    outFile.close();
}

std::vector<double> FileHolder::getCol(int column)
{
    std::vector<double> vec_column;
    for (int i=0;i< m_massPoints.size();i++) {
        vec_column.push_back(m_rates[i][column]);
    }
    return vec_column;
}

double FileHolder::getRateByMass(double mass, int col)
{
    int nrPoints = m_massPoints.size();
    for (int im=0;im<nrPoints;im++)
    {
        if (mass == m_massPoints[im])
        {
            return m_rates[im][col];
        }
    }
    return 0.;
}

void FileHolder::setRate(double mass, int col, double val)
{
    for (int imass=0;imass<(int)m_massPoints.size();imass++)
    {
        if (mass == m_massPoints[imass])
        {
            m_rates[imass][col] = val;
            break;
        }
    }
}

void FileHolder::copy(FileHolder& alt)
{
    alt.m_nrCols=m_nrCols;
    alt.m_rates=m_rates;
    alt.m_massPoints=m_massPoints;
}

void FileHolder::removeMass(double mass)
{
    std::vector<double> newPoints;
    std::vector<std::vector<double> > new_rates;
    int nrPoints = m_massPoints.size();
    for (int i=0;i<nrPoints;i++) {
        if ( m_massPoints.at(i) == mass ) continue;
        newPoints.push_back(m_massPoints[i]);
        new_rates.push_back(m_rates[i]);
    }
    m_massPoints = newPoints;
    m_rates      = new_rates;
}

void FileHolder::addMass(double mass, std::vector<double>& rate)
{
    //std::cout << "adding mass: " << mass << std::endl;
    std::vector<double> newPoints;
    std::vector<std::vector<double> > newm_rates;

    bool found = false;
    int nrPoints = m_massPoints.size();
    for (int i=0;i<nrPoints;i++)
    {
        if ((m_massPoints[i] < mass && !found) || (m_massPoints[i] > mass && found))
        {
            newPoints.push_back(m_massPoints[i]);
            newm_rates.push_back(m_rates[i]);
        }
        else if (m_massPoints[i] > mass && !found) 
        {
            //std::cout << "found" << std::endl;

            newPoints.push_back(mass);
            newm_rates.push_back(rate);

            newPoints.push_back(m_massPoints[i]);
            newm_rates.push_back(m_rates[i]);

            found = true;
        }
    }
    m_massPoints = newPoints;
    m_rates = newm_rates;
}

void FileHolder::addFile(FileHolder& other)
{
    std::vector<double> otherPoints = other.m_massPoints;
    std::vector<std::vector<double> > otherm_rates = other.m_rates;
    for (int i=0;i<(int)otherPoints.size();i++)
    {
        addMass(otherPoints[i], otherm_rates[i]);
    }
}

void FileHolder::addCol(std::vector<double>& addm_rates, int col)
{
    int nrPoints = m_massPoints.size();
    for (int i=0;i<nrPoints;i++)
    {
        std::vector<double> newm_rates;
        for (int j=0;j<col;j++)
        {
            newm_rates.push_back(m_rates[i][j]);
        }
        newm_rates.push_back(addm_rates[i]);
        for (int j=col;j<(int)m_rates[i].size();j++)
        {
            newm_rates.push_back(m_rates[i][j]);
        }
        m_rates[i]=newm_rates;
    }
    m_nrCols++;
}

void FileHolder::useMedian(int col)
{
    //std::cout << "in use median" << std::endl;
    set<double> vals;
    int nrm_rates = m_rates.size();
    for (int i=0;i<nrm_rates;i++)
    {
        vals.insert(m_rates[i][col]);
    }
    set<double>::iterator itr=vals.begin();
    int nrVals = vals.size();
    for (int i=0;i<nrVals/2;i++)
    {
        itr++;
    }

    for (int i=0;i<nrm_rates;i++)
    {
        m_rates[i][col] = *itr;
    }
    //std::cout << "out use median" << std::endl;
}
