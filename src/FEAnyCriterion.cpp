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
#include "FEAnyCriterion.h"
#include <FECore/FEModel.h>

BEGIN_FECORE_CLASS(FEAnyCriterion, FEMeshAdaptorCriterion)
    ADD_PARAMETER(m_max, "max");  // triggers remeshing of set if any value above this
    ADD_PARAMETER(m_min, "min");  // triggers remeshing of set if any below above this
    ADD_PARAMETER(m_value, "value");  // scale of adaption
	ADD_PARAMETER(m_minStride, "minimum_stride");  // minimum number of steps between consecutive adaptions
	ADD_PROPERTY(m_data, "data");
END_FECORE_CLASS();

FEAnyCriterion::FEAnyCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_value = 1.0;
	m_min = -1.0e37;
	m_max = 1.0e37;
	m_minStride = 1;
	m_last = -999999;
	m_count = 0;
}

FEMeshAdaptorSelection FEAnyCriterion::GetElementSelection(FEElementSet* elemSet)
{

	FEMesh& mesh = GetFEModel()->GetMesh();

	FEMeshAdaptorSelection elemList;	

	// increment eval counter
	m_count++;

	// add elements to selection and set scale value (not data_value)
	bool anyDataInSelection = false;
	if ((m_count - m_last) >= m_minStride) {
		for (FEElementIterator it(&mesh, elemSet); it.isValid(); ++it)
		{
			FEElement& el = *it;
			int nint = el.GaussPoints();
			for (int j = 0; j < nint; ++j)
			{
				double data_value;
				m_data->GetMaterialPointValue(*(el.GetMaterialPoint(j)), data_value);
				if ((data_value > m_min) && (data_value < m_max)) anyDataInSelection = true;

			}
			if (el.isActive()) elemList.push_back(el.GetID(), m_value);
		}
	}
	if (!anyDataInSelection) return FEMeshAdaptorSelection();  // return empty list
	
	// something to adapt
	m_last = m_count;  // note that this evaluation count returned non-empty set
	return elemList;  // return the set
}
