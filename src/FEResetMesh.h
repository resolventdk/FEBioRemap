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
#include <FECore/FEMeshAdaptor.h>


class FEModel;
class FEMeshDataInterpolator;
class FEMeshTopo;
class FEBoundaryCondition;
class FESurfaceLoad;
class FESurfacePairConstraint;
class FESurface;
class FENodeSet;
class FEDomainMap;
class FESolidElement;

//-----------------------------------------------------------------------------
// 
// Base class for mesh refinement algorithms that also resets mesh deformation:
// 
// Adapts the mesh according to a criteria and further resets the element 
// deformation in the selection.
// To account for the deformation zeroed, the deformation gradient is computed
// before resetting and stored in a user field, which is then mapped to the 
// new mesh so that it may be used with a prestrain material.
//
// The implementation is a modification of the class "FERefineMesh" from
// commit 682d9b14cc5ea09f99e54ccc85a06168b3beffd3.
//
// The mapping here is based on the current mesh configuration, which further
// required that the search trees was be modified to allow this.
// 
// Solution mapping works by interpolating results from nodes in the old
// mesh to either nodes or integration points in the new mesh:
// 
// 1. For integration point variables�the solution variables at the nodes
// of the old mesh is obtained by extrapolating and superpositioning
// values from the integration points to the nodes.
// 
// 2. Then the location of each point in the new mesh is obtained with
// respect to the old mesh.
// 
// 3. The variables are then interpolated from the nodes of the old
// element to the points in the new model.
// 
// TODO: mapping of none nodal fields is diffusive, hence ensure it 
// is only done for elements actually remeshed.
// 
//----------------------------------------------------------------------------

class FEResetMesh : public FEMeshAdaptor
{
protected:
	// Supported transfer methods for mapping data between meshes
	enum TransferMethod {
		TRANSFER_SHAPE,
		TRANSFER_MLQ
	};
	
public:
	FEResetMesh(FEModel* fem);
	~FEResetMesh();

	// Apply mesh refinement
	bool Apply(int iteration);

protected:
	// Derived classes need to override this function
	// The return value should be true if the mesh was refined
	// or false otherwise. 
	virtual bool ResetMesh() = 0;

protected:
	bool BuildMeshTopo();
	void CopyMesh();

	bool BuildMapData();
	void TransferMapData();
	void ClearMapData();

	bool BuildDomainMapData();
	bool BuildDomainMapData(FEDomain& dom, int domIndex);
	void TransferDomainMapData();

	bool BuildUserMapData();
	void TransferUserMapData();

	bool createNodeDataMap(FEDomain& dom, FEDomainMap* map, FEDomainMap* nodeMap);
	FEDomainMap* NodeToElemData(FEModel& fem, FEDomain& dom, FEDomainMap* nodeMap, FEMeshDataInterpolator* dataMapper, Storage_Fmt dataFormat);

	bool UpdateDisplacementMap();
	bool UpdateGradientMap();
	bool RestoreGradientMap();
	void ClearGradientMap();

	// for debugging only
	bool GradientsFromDisplacementMap();
	double defgrad(FESolidElement& el, std::vector<vec3d>& X, std::vector<vec3d>& u, mat3d& F, int n);
	double invjac0(FESolidElement& el, std::vector<vec3d>& r0, double Ji[3][3], int n);

protected:
	FEMeshTopo*	m_topo;		//!< mesh topo structure

	std::vector<vec3d>  m_nodeDispl;    // node displacements before reset

	int		m_maxiter;		// max nr of iterations per time step
	int		m_maxelem;		// max nr of elements

	int		m_transferMethod;		//!< method for transferring data between meshes
	bool	m_bmap_data;			//!< map data flag
	int		m_nnc;					//!< nearest-neighbor-count for MLQ transfer method
	int		m_nsdim;				//!< nearest-neighbor search dimension (2 or 3)

	string m_Dmap_name;     // name of user defined node displacement map 
	string m_Fmap_name;      // name of user defined deformation gradient map

	bool m_bmap_current;   // map field using current mesh configuration?

	double m_atol;         // absolute geometric tolerance for node and element search

	FEDomainMap* m_FmapCopy;  // copy of "old" gradient map before refinement 

	FEMesh*	m_meshCopy;		//!< copy of "old" mesh, before refinement

	std::vector< std::vector<FEDomainMap*> >	m_domainMapList;	// list of nodal data for each domain
	std::vector< FEDomainMap* >	m_userDataList;						// list of nodal data for user-defined mesh data

	DECLARE_FECORE_CLASS();
};
