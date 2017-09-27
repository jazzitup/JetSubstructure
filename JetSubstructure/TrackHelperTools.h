#ifndef __TRACKHELPERTOOLS_H__
#define __TRACKHELPERTOOLS_H__

#include <string>
#include <vector>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <iostream>
#include "BaseClass.h"

using namespace std;

namespace TrackHelperTools
{
	int getTypeTruth(int barcode, int pdg, int status);
}
#endif
