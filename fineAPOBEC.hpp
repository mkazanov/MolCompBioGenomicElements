//
//  fineAPOBEC.hpp
//  MolCompBioGenomicElements
//
//  Created by Marat Kazanov on 13/02/2022.
//  Copyright © 2022 Marat Kazanov. All rights reserved.
//

#ifndef fineAPOBEC_hpp
#define fineAPOBEC_hpp

#define SAMPLE_COLUMN_NO 1
#define THRESHOLD_HIGH_COLUMN_NO  2   //Column number (4.0 threshold)
#define THRESHOLD_MIDDLE_COLUMN_NO    3   //Column number (4.3 threshold)
#define THRESHOLD_LOW_COLUMN_NO   4   //Column number (4.6 threshold)

#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

vector<string> GetAPOBECsamples(int threshold);

#endif /* fineAPOBEC_hpp */
