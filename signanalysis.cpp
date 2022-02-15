  //
//  apobec.cpp
//  APOBECXP
//
//  Created by Marat Kazanov on 28/11/2018.
//  Copyright © 2018 Marat Kazanov. All rights reserved.
//

#include <iostream>

#include "signanalysis.hpp"
#include <map>
#include <fstream>
#include "options.h"
#include <ctime>
#include <cstring>
#include "dna.hpp"

CResultsKey::CResultsKey(string motif_, string cancer_, string sample_, int rtbin_, int expbin_)
{
    motif = motif_;
    cancer = cancer_;
    sample = sample_;
    rtbin = rtbin_;
    expbin = expbin_;
}

CResultsKey::CResultsKey(string cancer_, string sample_, int rtbin_, int expbin_)
{
    motif = "X";
    cancer = cancer_;
    sample = sample_;
    rtbin = rtbin_;
    expbin = expbin_;
}

CResultsValue::CResultsValue(unsigned long mutCnt_, unsigned long leadingCnt_, unsigned long laggingCnt_)
{
    mutCnt = mutCnt_;
    leadingCnt = leadingCnt_;
    laggingCnt = laggingCnt_;
}


CResultsValue::CResultsValue(unsigned long mutCnt_, unsigned long plusStrandConsistent_, unsigned long minusStrandConsistent_, unsigned long plusStrandAll_, unsigned long minusStrandAll_)
{
    mutCnt = mutCnt_;
    plusStrandConsistent = plusStrandConsistent_;
    minusStrandConsistent = minusStrandConsistent_;
    plusStrandAll = plusStrandAll_;
    minusStrandAll = minusStrandAll_;
}

CResultsValue::CResultsValue(unsigned long mutCnt_, unsigned long leadingCnt_, unsigned long laggingCnt_,unsigned long plusStrandConsistent_, unsigned long minusStrandConsistent_, unsigned long plusStrandAll_, unsigned long minusStrandAll_)
{
    mutCnt = mutCnt_;
    leadingCnt = leadingCnt_;
    laggingCnt = laggingCnt_;
    plusStrandConsistent = plusStrandConsistent_;
    minusStrandConsistent = minusStrandConsistent_;
    plusStrandAll = plusStrandAll_;
    minusStrandAll = minusStrandAll_;
}

void CSignatureAnalysis:: ClassifyMutations(CMutations& muts, vector<string> cancers, vector<string> samples, CHumanGenome* phuman_)
{
    CHumanGenome* phuman;
    
    // Load genome
    if(phuman_)
        phuman = phuman_;
    else
    {
        phuman = new CHumanGenome();
        phuman->InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");
    }
    
    // Filter mutations
    muts.FilterBySignature(signatureMuts,signatures,(*phuman),cancers,samples,&otherMuts);
    
    signatureMuts.GetUniqueCancersSamples();
    
    // Free human genome
    if(!phuman_)
    {
        for(int i=0;i<phuman->chrCnt;i++)
            delete phuman->dna[i];
        delete phuman->dna;
    }
}

void CSignatureAnalysis::ClearMutations()
{
    signatureMuts.ClearMutations();
    otherMuts.ClearMutations();
}
















