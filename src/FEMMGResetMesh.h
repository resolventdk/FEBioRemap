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
#pragma once
#include <vector>
#include "FEResetMesh.h"
#include <FECore/FEFunction1D.h>
#include <FECore/FEModelParam.h>
#include <FECore/FESolidElement.h>
#include <FECore/FEDomainMap.h>

class FEMMGResetMesh : public FEResetMesh
{
	class MMG;

public:
	FEMMGResetMesh(FEModel* fem);

	bool Init() override;
	bool ResetMesh() override;

private:
	FEMeshAdaptorCriterion* GetCriterion() { return m_criterion; }

private:

	bool    m_bremesh;  // change the mesh or just reset? 

	int		m_maxiter;
	bool	m_relativeSize;
	bool	m_meshCoarsen;
	bool	m_normalizeData;

	double	m_hmin;		// minimum element size
	double	m_hausd;	// Hausdorff value
	double	m_hgrad;	// gradation
	int	    m_nreg;	    // use normal regularization yes/no 1/0

	int m_attachedRid;  // adapted selection contains nodes attached to this rigid body id, we support only a single and subtract rigid body motion
	const int m_nid_offset = 1000;   // mmg point reference offset for to old node ids. Must be out of range of element and surface ids. 

	FEMeshAdaptorCriterion*	m_criterion;

	FEFunction1D*	m_sfunc;	// sizing function

	MMG*	mmg;

	friend class MMG;

	DECLARE_FECORE_CLASS();
};
