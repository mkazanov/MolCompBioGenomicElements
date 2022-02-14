//
//  main.cpp
//  MolCompBioGenomicElements
//
//  Created by Marat Kazanov on 07/02/2022.
//  Copyright © 2022 Marat Kazanov. All rights reserved.
//

#include <iostream>
#include "gelements.hpp"
#include "mutation.hpp"
#include "options.h"
#include "fineAPOBEC.hpp"

int main(int argc, const char * argv[]) {
    
    int i;
    string path;
    
    CHumanGenome human;
    human.InitializeHuman("37", HUMAN_PATH, ".fa", "FASTA");

    vector<string> cancers,curcancers;
    vector<string> samples;
    vector<string>::iterator cit;
    
    cancers.push_back("Bladder-TCC");
    cancers.push_back("Breast-AdenoCA");
    cancers.push_back("Breast-DCIS");
    cancers.push_back("Breast-LobularCA");
    cancers.push_back("Cervix-AdenoCA");
    cancers.push_back("Cervix-SCC");
    cancers.push_back("Head-SCC");
    cancers.push_back("Lung-AdenoCA");
    cancers.push_back("Lung-SCC");

    ofstream f;
    path = string(RESULTS_FOLDER) + "/samples.txt";
    f.open(path.c_str());
    
    CMutations m;
    set<CCancerSample>::iterator it;
    for(cit=cancers.begin();cit!=cancers.end();cit++)
    {
        curcancers.push_back((*cit));
        path = CANCER_MUTATIONS_PATH + (*cit) + "_MAF_Oct12_2016_sorted_aTn_anz4.txt";
        m.LoadMutations(CANCER_MUTATIONS_FILEFORMAT, path, curcancers, samples, 1, &human);
        m.GetUniqueCancersSamples();
        
        for(it=m.cancerSample.begin();it!=m.cancerSample.end();it++)
        {
            f << *cit << '\t' << it->sample << '\n';
        }
        m.ClearMutations();
        curcancers.clear();
    }

    f.close();
    return 0;
    
    CGenomicElements g4;
    
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr1_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr2_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr3_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr4_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr5_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr6_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr7_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr8_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr9_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr10_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr11_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr12_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr13_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr14_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr15_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr16_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr17_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr18_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr19_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr20_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr21_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chr22_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chrX_GQ.tsv", NONDDNA_G4);
    g4.LoadGenomicElements("/Users/mar/BIO/BIODATA/nonBDNADB/human_hg19.tsv/gq/tsv/chrY_GQ.tsv", NONDDNA_G4);

    for(i=0;i<24;i++)
    {
        cout << g4.gelements[i].size() << '\n';
    }
    
    int res;
    g4.PrepareForSearch();
    res = g4.GetNearestElementByPos(0, 11029);
    cout << res << '\n';
    return 0;
}
