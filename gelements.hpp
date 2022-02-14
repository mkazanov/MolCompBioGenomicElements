//
//  gelements.hpp
//  MolCompBioGenomicElements
//
//  Created by Marat Kazanov on 07/02/2022.
//  Copyright © 2022 Marat Kazanov. All rights reserved.
//

#ifndef gelements_hpp
#define gelements_hpp

#define NONDDNA_G4 0

#include <string>
#include <bitset>
#include <map>
#include <vector>

using namespace std;

class CGenomicElementFileFormat{
public:
    char delimiter;
    int chrNo;
    int startposNo;
    int endposNo;
    int strandNo;
    int nameNo;
    string infoNos;
    int isHeader;
    CGenomicElementFileFormat(){};
    CGenomicElementFileFormat(char delimiter_,
                   int chrNo_,
                   int startposNo_,
                   int endposNo_,
                   int strandNo_,
                   int nameNo_,
                   string infoNos_,
                   int isHeader_);
};

class CGenomicElement {
public:
    unsigned long startpos;
    unsigned long endpos;
    char strand;
    string elementName;
    string info;
    CGenomicElement(unsigned long startpos_, unsigned long endpos_);
    CGenomicElement(string starpos_, string endpos_, string strand_, string elementName_, string info_);
    bool operator< (const CGenomicElement &right) const
    {
        if (startpos < right.startpos)
            return true;
        else if (startpos == right.startpos)
            return endpos > right.endpos;
        else
            return false;
    }
};

class CGenomicElements {
public:
    map<int,CGenomicElementFileFormat> fileFormat;
    CGenomicElements();
    ~CGenomicElements();
    vector<CGenomicElement>* gelements;
    vector<unsigned long>* gestarts;
    vector<unsigned long>* geends;
    bitset<250000000>* gebits;
    void LoadGenomicElements(string path, int fileFormatType);
    void PrepareForSearch();
    unsigned long GetNearestElementByPos(int chrNum, unsigned long pos);
    void SaveToFile(string path);
    //void MergeIntervals(vector<CHumanGene>& ret);
};

#endif /* gelements_hpp */
