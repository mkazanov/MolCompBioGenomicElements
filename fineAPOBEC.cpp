//
//  fineAPOBEC.cpp
//  MolCompBioGenomicElements
//
//  Created by Marat Kazanov on 13/02/2022.
//  Copyright © 2022 Marat Kazanov. All rights reserved.
//

#include "fineAPOBEC.hpp"
#include "options.h"
#include <fstream>
#include "service.h"

vector<string> GetAPOBECsamples(int threshold)
{
    string line;
    vector<string> resSamples;
    vector<string> flds;
    
    ifstream f(APOBEC_SAMPLES_LIST);
    if (!f.is_open())
    {
        printf("File with list of samples does not exists\n");
        exit(1);
    }
    
    getline(f, line);
    while(1)
    {
        getline(f, line);
        if(f.eof())
            break;
        flds = splitd(line,'\t');
        if(flds[threshold] == "1")
            resSamples.push_back(flds[SAMPLE_COLUMN_NO]);
    }
    
    return(resSamples);
}
