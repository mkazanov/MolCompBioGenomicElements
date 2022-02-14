//
//  gelements.cpp
//  MolCompBioGenomicElements
//
//  Created by Marat Kazanov on 07/02/2022.
//  Copyright © 2022 Marat Kazanov. All rights reserved.
//

#include "gelements.hpp"
#include <iostream>
#include <fstream>
#include "service.h"
#include "ghuman.hpp"

CGenomicElementFileFormat::CGenomicElementFileFormat(char delimiter_,
                                                     int chrNo_,
                                                     int startposNo_,
                                                     int endposNo_,
                                                     int strandNo_,
                                                     int nameNo_,
                                                     string infoNos_,
                                                     int isHeader_)
{
    delimiter = delimiter_;
    chrNo = chrNo_;
    startposNo = startposNo_;
    endposNo = endposNo_;
    strandNo = strandNo_;
    nameNo = nameNo_;
    infoNos = infoNos_;
    isHeader = isHeader_;
}

CGenomicElement::CGenomicElement(unsigned long startpos_, unsigned long endpos_)
{
    startpos = startpos_;
    endpos = endpos_;
}

CGenomicElement::CGenomicElement(string startpos_, string endpos_, string strand_, string elementName_, string info_)
{
    startpos = str2ul(startpos_);
    endpos = str2ul(endpos_);
    strand = strand_[0];
    elementName = elementName_;
    info = info_;
}

CGenomicElements::CGenomicElements()
{
    // Non-B DNA G4 file format
    fileFormat[NONDDNA_G4] = CGenomicElementFileFormat('\t', //separator
                                        0, // chr field num
                                        3, // startpos field num
                                        4, // endpos field num
                                        7, // strand field num
                                        -1, // name field num
                                        "5,6,8,9,10,11,12", // info field num
                                        1 // is header
                                        );
    
    gelements = new vector<CGenomicElement>[24];
    gestarts = new vector<unsigned long>[24];
    geends = new vector<unsigned long>[24];
    gebits = new bitset<250000000>[24];
}

CGenomicElements::~CGenomicElements()
{
    delete[] gelements;
    delete[] gestarts;
    delete[] geends;
    delete[] gebits;
}

void CGenomicElements::LoadGenomicElements(string path, int fileFormatType)
{
    cout << "Loading genomic elements from file: " << path << '\n';
    
    int i;
    string line;
    int chrNum;
    string startpos;
    string endpos;
    string strand;
    string elementName;
    int infoNo;
    vector<string> infoNos;
    string info;
    vector<string> flds,hflds;
    string colname;
    
    ifstream f(path.c_str());
    if (!f.is_open())
    {
        printf("File not exists\n");
        return;
    }
    
    if(fileFormat[fileFormatType].isHeader)
        getline(f, line);
        hflds = splitd(line,fileFormat[fileFormatType].delimiter);
    
    while(getline(f, line))
    {
        if (line.length() != 0)
        {
            flds = splitd(line,fileFormat[fileFormatType].delimiter);
            chrNum = CHumanGenome::GetChrNum(flds[fileFormat[fileFormatType].chrNo]);
            if(chrNum < 0 || chrNum > 23)
                continue;
            startpos = flds[fileFormat[fileFormatType].startposNo];
            endpos = flds[fileFormat[fileFormatType].endposNo];
            
            if(fileFormat[fileFormatType].strandNo != -1)
                strand = flds[fileFormat[fileFormatType].strandNo][0];
            else
                strand = '\0';
            
            if(fileFormat[fileFormatType].nameNo != -1)
                elementName = flds[fileFormat[fileFormatType].nameNo];
            else
                elementName = "";
            
            info = "";
            if(fileFormat[fileFormatType].infoNos != "-1")
            {
               infoNos = splitd(fileFormat[fileFormatType].infoNos,',');
               for(i=0;i<infoNos.size();i++)
               {
                   infoNo = str2i(infoNos[i]);
                   if(hflds.size()==0)
                       colname = "col"+i2str(infoNo);
                   else
                       colname = hflds[infoNo];
                   info = info + colname + ":" + flds[infoNo] + ";";
               }
            }
            gelements[chrNum].push_back(CGenomicElement(startpos,
                                                        endpos,
                                                        strand,
                                                        elementName,
                                                        info));
            gestarts[chrNum].push_back(str2ul(startpos));
            geends[chrNum].push_back(str2ul(endpos));
        }
    }
    
    f.close();
}

void CGenomicElements::PrepareForSearch()
{
    int chrNum;
    vector<CGenomicElement>::iterator it;
    
    
    // set bits
    unsigned long i;
    for(chrNum=0;chrNum<24;chrNum++)
    {
        for(it=gelements[chrNum].begin();it!=gelements[chrNum].end();++it)
        {
            for(i=it->startpos;i<=it->endpos;i++)
                gebits[chrNum][i-1] = 1;
        }
    }
    
    for(chrNum=0;chrNum<24;chrNum++)
        sort(gestarts[chrNum].begin(),gestarts[chrNum].end());
    
    for(chrNum=0;chrNum<24;chrNum++)
        sort(geends[chrNum].begin(),geends[chrNum].end());
}

unsigned long CGenomicElements::GetNearestElementByPos(int chrNum, unsigned long pos)
{
    unsigned long dist1,dist2;
    vector<unsigned long>::iterator it;
    
    if(gebits[chrNum][pos-1] == 1)
        return(0);
    
    it = upper_bound(gestarts[chrNum].begin(),gestarts[chrNum].end(),pos);
    if(it != gestarts[chrNum].end())
        dist1 = *it - pos;
    else
        dist1 = 0;
    
    it = lower_bound(geends[chrNum].begin(),geends[chrNum].end(),pos);
    --it;
    if(it != geends[chrNum].begin())
        dist2 = pos - *it;
    else
        dist2 = 0;
    
    if(dist1 == 0 & dist2 == 0)
        return(-1);
    else if(dist1 !=0 & dist2 == 0)
        return(dist1);
    else if(dist1 == 0 & dist2 != 0)
        return(dist2);
    else if(dist1 < dist2)
        return(dist1);
    else
        return(dist2);
}
