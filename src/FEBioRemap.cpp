/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#include <FECore/sdk.h>
#include "FERemapPlot.h"
#include "FEMMGResetMesh.h"

#include "FEFixedTimeCriterion.h"
#include "FEQualityCriterion.h"
#include "FEAnyCriterion.h"

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

FECORE_PLUGIN void GetPluginVersion(int& major, int& minor, int& patch)
{
	major = 0;
	minor = 0;
	patch = 0;
}

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);

	// create the module
	const char* info = \
		"{ "
		"   \"title\" : \"FEBio Solution Re-mapping\","
		"   \"info\"  : \"Functionality for remapping solution fields onto remeshed mesh to prevent mesh distortion during large strain simulation.\","
		"   \"author\": \"Henrik Spietz\","
		"   \"version\": \"0.0\""
        "}";

	// mesh adapters
	REGISTER_FECORE_CLASS(FEMMGResetMesh, "mmg_reset");
	
	// criteria
	REGISTER_FECORE_CLASS(FEFixedTimeCriterion, "fixed dt");
	REGISTER_FECORE_CLASS(FEQualityCriterion, "quality");
	REGISTER_FECORE_CLASS(FEAnyCriterion, "min-max");

	// classes derived from FEPlotData
	REGISTER_FECORE_CLASS(FEPlotQuality, "quality");

}
