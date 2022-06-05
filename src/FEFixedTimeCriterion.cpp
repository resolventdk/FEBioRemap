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
#include "FEFixedTimeCriterion.h"
#include <FECore/FEModel.h>

BEGIN_FECORE_CLASS(FEFixedTimeCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_dt, "dt");  // how often to remesh
    ADD_PARAMETER(m_value, "value");  // scale of adaption
END_FECORE_CLASS();

FEFixedTimeCriterion::FEFixedTimeCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_tp = 0.0;
	m_value = 1.0;
}

FEMeshAdaptorSelection FEFixedTimeCriterion::GetElementSelection(FEElementSet* elemSet)
{

	FEMesh& mesh = GetFEModel()->GetMesh();

	double time = GetFEModel()->GetTime().currentTime;
	
	FEMeshAdaptorSelection elemList;	
	if ((time - m_tp) > (m_dt - 1e-6)) {  // if time between remeshes has passed
		// add elements to selection and set scale value
		m_tp = time; // update time				
		for (FEElementIterator it(&mesh, elemSet); it.isValid(); ++it)
		{
			FEElement& el = *it;
			if (el.isActive()) elemList.push_back(el.GetID(), m_value);
		}
	}
	return elemList;
}
