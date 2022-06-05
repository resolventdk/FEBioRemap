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
#include "math.h"
#include "FEQualityCriterion.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FESolidElement.h>

BEGIN_FECORE_CLASS(FEQualityCriterion, FEMeshAdaptorCriterion)
END_FECORE_CLASS();

FEQualityCriterion::FEQualityCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
}

bool FEQualityCriterion::GetMaterialPointValue(FEMaterialPoint& mp, double& value)
{

    FEModel& fem = *GetFEModel();

    FEMesh& mesh = fem.GetMesh();

    FESolidElement* el = dynamic_cast<FESolidElement*>(mp.m_elem);
    assert(el);

    // calculate jacobian
    double J[3][3] = { 0 };
    for (int i = 0; i < el->Nodes(); ++i)
    {
        const double& Gri = el->Gr(mp.m_index)[i];
        const double& Gsi = el->Gs(mp.m_index)[i];
        const double& Gti = el->Gt(mp.m_index)[i];

        const double& x = mesh.Node(el->m_node[i]).m_rt.x;
        const double& y = mesh.Node(el->m_node[i]).m_rt.y;
        const double& z = mesh.Node(el->m_node[i]).m_rt.z;

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
    value = 1.0 / eta;  // quality in (0, 1]

	return true;
}
