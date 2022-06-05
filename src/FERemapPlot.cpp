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
#include "FERemapPlot.h"
#include "FECore/FESolidDomain.h"
#include "FECore/FESolidElement.h"
#include "FECore/FEModel.h"
#include "FECore/FEAnalysis.h"

bool FEPlotQuality::Save(FEDomain& dom, FEDataStream& a)
{
	FESolidDomain* pdom = dynamic_cast<FESolidDomain*>(&dom);
	if (pdom == nullptr) return false;

	FEMesh* mesh = pdom->GetMesh();

	// store
	int NE = pdom->Elements();
	for (int i = 0; i < NE; ++i)
	{
		FESolidElement& e = pdom->Element(i);
		int nint = e.GaussPoints();
		int neln = e.Nodes();

		// nodal coordinates
		vec3d rt[FEElement::MAX_NODES];
		for (int k = 0; k < neln; ++k) rt[k] = mesh->Node(e.m_node[k]).m_rt;

		// average gauss points
		double v = 0.0;
		for (int j = 0; j < nint; ++j)
		{

			// calculate jacobian
			double J[3][3] = { 0 };
			for (int k = 0; k < neln; ++k)
			{
				const double& Gri = e.Gr(j)[k];
				const double& Gsi = e.Gs(j)[k];
				const double& Gti = e.Gt(j)[k];

				const double& x = rt[k].x;
				const double& y = rt[k].y;
				const double& z = rt[k].z;

				J[0][0] += Gri * x; J[0][1] += Gsi * x; J[0][2] += Gti * x;
				J[1][0] += Gri * y; J[1][1] += Gsi * y; J[1][2] += Gti * y;
				J[2][0] += Gri * z; J[2][1] += Gsi * z; J[2][2] += Gti * z;
			}

			// calculate the determinant
			double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
				       + J[0][1] * (J[1][2] * J[2][0] - J[2][2] * J[1][0])
				       + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

			// metric cf J.M.Escobar 2003: Simultaneous untangling and smoothing of tetrahedral meshes
			double norm2 = J[0][0] * J[0][0]
				         + J[0][1] * J[0][1]
		                 + J[0][2] * J[0][2]
				         + J[1][0] * J[1][0]
				         + J[1][1] * J[1][1]
				         + J[1][2] * J[1][2]
				         + J[2][0] * J[2][0]
				         + J[2][1] * J[2][1]
				         + J[2][2] * J[2][2];
			double eta = norm2 / (3.0 * pow(abs(det), 2.0 / 3.0));  // distortion in [1, infty)
			v += 1.0 / eta;  // quality in (0, 1]

		}
		v /= (double)nint;

		a << v;
	}
	return true;
}