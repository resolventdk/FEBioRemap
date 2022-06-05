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
#include "FEResetMesh.h"

#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FEEdgeList.h>
#include <FECore/FEElementList.h>
#include <FECore/FEFaceList.h>
#include <FECore/FEFixedBC.h>
#include <FECore/FEPrescribedDOF.h>
#include <FECore/FEMeshTopo.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FESurfacePairConstraint.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FEDomainMap.h>
#include <FECore/FESurfaceMap.h>
#include <FECore/DumpMemStream.h>
#include <FECore/log.h>

#include "FELeastSquaresInterpolator.h"
#include "FEMeshShapeInterpolator.h"
#include "FEDomainShapeInterpolator.h"

#include <FEBioMech/FEConstPrestrain.h>
#include <FEBioMech/FEDeformationMapGenerator.h>

BEGIN_FECORE_CLASS(FEResetMesh, FEMeshAdaptor)
	//ADD_PARAMETER(m_maxiter, "max_iters");
	ADD_PARAMETER(m_maxelem, "max_elements");
	ADD_PARAMETER(m_bmap_data, "map_data");
	ADD_PARAMETER(m_nnc, "nnc");
	ADD_PARAMETER(m_nsdim, "nsdim");
	ADD_PARAMETER(m_transferMethod, "transfer_method");
	ADD_PARAMETER(m_Fmap_name, "gradient_map");
	ADD_PARAMETER(m_Dmap_name, "displacement_map");
END_FECORE_CLASS();

FEResetMesh::FEResetMesh(FEModel* fem) : FEMeshAdaptor(fem), m_topo(nullptr)
{
	m_meshCopy = nullptr;
	m_bmap_data = false;
	m_transferMethod = TRANSFER_SHAPE;
	m_nnc = 8;
	m_nsdim = 3;

	m_maxiter = 1;
	m_maxelem = -1;

	m_bmap_current = false; // do not use current mesh configuration for mapping
}

FEResetMesh::~FEResetMesh()
{
	if (m_meshCopy) delete m_meshCopy;
	m_meshCopy = nullptr;

	ClearMapData();
}

// Apply mesh refinement and reset deformation
bool FEResetMesh::Apply(int iteration)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// check for max iterations
	if ((m_maxiter > 0) && (iteration >= m_maxiter))
	{
		feLog("Skipping refinement: Max iterations reached.");
		return false;
	}

	// see if we reached max elements
	if ((m_maxelem > 0) && (mesh.Elements() >= m_maxelem))
	{
		feLog("Skipping refinement: Element limit reached.\n");
		return false;
	}

	// build the mesh-topo
	if (BuildMeshTopo() == false)
	{
		throw std::runtime_error("Error building topo structure.");
	}

	// update deformation gradient map
	if (UpdateGradientMap() == false)
	{
		throw std::runtime_error("Failed to update gradient map before mesh refinement.");
	}
	
	// build data maps
	feLog("-- Building map data:\n");
	if (BuildMapData() == false)
	{
		throw std::runtime_error("Failed mapping data.");
	}

	// refine the mesh (This is done by sub-classes)
	feLog("-- Starting Mesh refinement.\n");
	if (ResetMesh() == false)
	{

		feLog("Mesh not reset.");

		if (RestoreGradientMap() == false) {
			throw std::runtime_error("Failed to restore gradient map.");
		}

		return false;

	}
	feLog("-- Mesh refinement completed.\n");

	// map data to new mesh
	feLog("-- Transferring map data to new mesh:\n");
	TransferMapData();

	// update displacement map (must be after existing displacements have been transfered)
	if (UpdateDisplacementMap() == false)
	{
		throw std::runtime_error("Failed to update displacement map after mesh refinement.");
	}

	//// only in dummy mode to check mapping error
	//if (GradientsFromDisplacementMap() == false)
	//{
	//	throw std::runtime_error("Failed to generate gradients from displacement map after mesh refinement.");
	//}

	// update the model
	UpdateModel();

	// print some mesh statistics
	int NN = mesh.Nodes();
	int NE = mesh.Elements();
	feLog("\n Mesh Statistics:\n");
	feLog(" \tNumber of nodes    : %d\n", NN);
	feLog(" \tNumber of elements : %d\n", NE);
	feLog("\n");

	// clean up
	ClearGradientMap();  // delete the backup

	// all done!
	return true;
}

void FEResetMesh::ClearMapData()
{
	// clear domain maps
	for (size_t i = 0; i < m_domainMapList.size(); ++i)
	{
		std::vector<FEDomainMap*>& map_i = m_domainMapList[i];
		for (size_t j = 0; j < map_i.size(); ++j) delete map_i[j];
	}
	m_domainMapList.clear();

	// clear user maps
	for (int i = 0; i < m_userDataList.size(); ++i) delete m_userDataList[i];
	m_userDataList.clear();
}

bool FEResetMesh::BuildMeshTopo()
{
	FEModel& fem = *GetFEModel();
	if (m_topo) { delete m_topo; m_topo = nullptr; }
	m_topo = new FEMeshTopo;
	return m_topo->Create(&fem.GetMesh());
}

void FEResetMesh::CopyMesh()
{
	if (m_meshCopy) delete m_meshCopy;

	m_meshCopy = new FEMesh(nullptr);
	FEMesh& mesh = GetFEModel()->GetMesh();
	m_meshCopy->CopyFrom(mesh);
}

bool FEResetMesh::createNodeDataMap(FEDomain& dom, FEDomainMap* map, FEDomainMap* nodeMap)
{
	FEDataType dataType = map->DataType();
	int dataSize = 0;
	switch (dataType)
	{
	case FEDataType::FE_DOUBLE: dataSize = 1; break;
	case FEDataType::FE_VEC3D: dataSize = 3; break;
	case FEDataType::FE_MAT3D: dataSize = 9; break;
	case FEDataType::FE_MAT3DS: dataSize = 6; break;
	default:
		assert(false);
		return false;
	}

	// temp storage 
	double si[FEElement::MAX_INTPOINTS * 9];
	double sn[FEElement::MAX_NODES * 9];

	// allocate node data
	int NN = dom.Nodes();
	vector<double> nodeData(NN*dataSize);

	// build tag list
	vector<int> tag(NN, 0);
	int NE = dom.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int ne = e.Nodes();
		for (int k = 0; k < ne; ++k)
		{
			tag[e.m_lnode[k]]++;
		}
	}

	// get the data format
	int dataFormat = map->StorageFormat();
	//if ((dataFormat != FMT_MATPOINTS) && (dataFormat != FMT_MULT)) return false;

	// loop over all elements
	for (int i = 0; i < NE; ++i)
	{
		FEElement& e = dom.ElementRef(i);
		int ne = e.Nodes();

		int ni;			
		switch (dataFormat) {
		case FMT_MATPOINTS: ni = e.GaussPoints(); break;
		case FMT_MULT: ni = ne; break;
		case FMT_ITEM: ni = 1; break;
		default: feLogError("Unsupported dataFormat for mapping!\n"); return false;  break;
		}

		for (int j = 0; j < dataSize; ++j)
		{
			// get the integration point values
			for (int k = 0; k < ni; ++k)
			{
				switch (dataType)
				{
				case FEDataType::FE_DOUBLE:
					si[k] = map->value<double>(i, k);
					break;
				case FEDataType::FE_VEC3D:
				{
					vec3d v = map->value<vec3d>(i, k);
					if (j == 0) si[k] = v.x;
					if (j == 1) si[k] = v.y;
					if (j == 2) si[k] = v.z;
				}
				break;
				case FEDataType::FE_MAT3D:
				{
					mat3d v = map->value<mat3d>(i, k);
					int LUT[9][2] = { {0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {1,2}, {2,0}, {2,1}, {2,2} };
					si[k] = v(LUT[j][0], LUT[j][1]);
				}
				break;
				case FEDataType::FE_MAT3DS:
				{
					mat3ds v = map->value<mat3ds>(i, k);
					int LUT[6][2] = { {0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2} };
					si[k] = v(LUT[j][0], LUT[j][1]);
				}
				break;
				}
			}

			// project to nodes
			switch (dataFormat) {
			case FMT_MATPOINTS:
				e.project_to_nodes(si, sn); 
				break;
			case FMT_MULT: 
				for (int k = 0; k < ne; ++k) 
					sn[k] = si[k]; 
				break;
			case FMT_ITEM: ni = 1; 
				for (int k = 0; k < ne; ++k)
					sn[k] = si[0];
				break;
			}

			for (int k = 0; k < ne; ++k)
			{
				nodeData[e.m_lnode[k] * dataSize + j] += sn[k];
			}
		}
	}

	// normalize data
	for (int i = 0; i < NN; ++i)
	{
		if (tag[i] > 0)
		{
			for (int j = 0; j < dataSize; ++j)
				nodeData[i*dataSize + j] /= (double)tag[i];
		}
	}

	// check node data map
	if (nodeMap->StorageFormat() != Storage_Fmt::FMT_NODE) return false;
	if (nodeMap->DataType() != dataType) return false;
	if (nodeMap->DataCount() != NN) return false;

	// assign data
	for (int i = 0; i < NN; ++i)
	{
		switch (dataType)
		{
		case FEDataType::FE_DOUBLE: nodeMap->setValue(i, nodeData[i]); break;
		case FEDataType::FE_VEC3D:
		{
			vec3d v;
			v.x = nodeData[i*dataSize];
			v.y = nodeData[i*dataSize + 1];
			v.z = nodeData[i*dataSize + 2];
			nodeMap->setValue(i, v);
		}
		break;
		case FEDataType::FE_MAT3D:
		{
			mat3d v;
			v(0, 0) = nodeData[i*dataSize]; v(0, 1) = nodeData[i*dataSize + 1]; v(0, 2) = nodeData[i*dataSize + 2];
			v(1, 0) = nodeData[i*dataSize + 3]; v(1, 1) = nodeData[i*dataSize + 4]; v(1, 2) = nodeData[i*dataSize + 5];
			v(2, 0) = nodeData[i*dataSize + 6]; v(2, 1) = nodeData[i*dataSize + 7]; v(2, 2) = nodeData[i*dataSize + 8];
			nodeMap->setValue(i, v);
		}
		break;
		case FEDataType::FE_MAT3DS:
		{
			mat3ds v;
			v(0, 0) = nodeData[i*dataSize];
			v(0, 1) = nodeData[i*dataSize + 1];
			v(0, 2) = nodeData[i*dataSize + 2];
			v(1, 1) = nodeData[i*dataSize + 3];
			v(1, 2) = nodeData[i*dataSize + 4];
			v(2, 2) = nodeData[i*dataSize + 5];
			nodeMap->setValue(i, v);
		}
		break;
		}
	}

	return true;
}

FEDomainMap* FEResetMesh::NodeToElemData(FEModel& fem, FEDomain& dom, FEDomainMap* nodeMap, FEMeshDataInterpolator* dataMapper, Storage_Fmt dataFormat)
{
	assert(nodeMap->StorageFormat() == Storage_Fmt::FMT_NODE);
	assert(dataFormat == Storage_Fmt::FMT_MULT || dataFormat == Storage_Fmt::FMT_MATPOINTS);

	FEDataType dataType = nodeMap->DataType();
	int dataSize = 0;
	switch (dataType)
	{
	case FEDataType::FE_DOUBLE: dataSize = 1; break;
	case FEDataType::FE_VEC3D: dataSize = 3; break;
	case FEDataType::FE_MAT3D: dataSize = 9; break;
	case FEDataType::FE_MAT3DS: dataSize = 6; break;
	default:
		assert(false);
		throw std::runtime_error("Error in FEMMGRemesh::MMG::NodeToElemData");
		return nullptr;
	}

	// create new domain map
	FEDomainMap* elemMap = new FEDomainMap(nodeMap->DataType(), dataFormat);
	FEElementSet* eset = new FEElementSet(&fem);
	eset->Create(&dom);
	elemMap->Create(eset);

	// count nr of source points
	int NSourcePoints = nodeMap->DataCount();

	// count nr of target points
	int NTargetPoints = 0;
	switch (elemMap->StorageFormat())
	{
	case Storage_Fmt::FMT_MULT:
		for (int i = 0; i < dom.Elements(); ++i)
		{
			FEElement& el = dom.ElementRef(i);
			NTargetPoints += el.Nodes();
		}
		break;
	case Storage_Fmt::FMT_MATPOINTS:
		for (int i = 0; i < dom.Elements(); ++i)
		{
			FEElement& el = dom.ElementRef(i);
			NTargetPoints += el.GaussPoints();
		}
		break;
	}

	// allocate interpolation buffers
	vector<double> srcData(NSourcePoints);
	vector<double> trgData(NTargetPoints);

	// help buffer(npts, ncomps)
	vector< vector<double> > mappedData(NTargetPoints, vector<double>(9, 0.0));

	// loop components of data
	for (int l = 0; l < dataSize; ++l)
	{
		// pack data from source points
		for (int i = 0; i < NSourcePoints; ++i)
		{
			double vm = 0.0;
			switch (dataType)
			{
			case FEDataType::FE_DOUBLE: vm = nodeMap->value<double>(0, i); break;
			case FEDataType::FE_VEC3D:
			{
				vec3d v = nodeMap->value<vec3d>(0, i);
				if (l == 0) vm = v.x;
				if (l == 1) vm = v.y;
				if (l == 2) vm = v.z;
			}
			break;
			case FEDataType::FE_MAT3D:
			{
				int LUT[9][2] = { {0,0}, {0,1}, {0,2}, {1,0}, {1,1}, {1,2}, {2,0}, {2,1}, {2,2} };
				mat3d v = nodeMap->value<mat3d>(0, i);
				vm = v(LUT[l][0], LUT[l][1]);
			}
			break;
			case FEDataType::FE_MAT3DS:
			{
				int LUT[6][2] = { {0,0}, {0,1}, {0,2}, {1,1}, {1,2}, {2,2} };
				mat3ds v = nodeMap->value<mat3ds>(0, i);
				vm = v(LUT[l][0], LUT[l][1]);
			}
			break;
			default:
				assert(false);
			}
			srcData[i] = vm;
		}

		// interpolate from source to target nodes
		dataMapper->Map(trgData, [&srcData](int sourcePoint) {
			return srcData[sourcePoint];
		});

		// unpack data to target points
		for (int i = 0; i < NTargetPoints; ++i)
		{
			mappedData[i][l] = trgData[i];
		}
	}

	// write mapped data to domain map
	int n = 0;
	for (int i = 0; i < dom.Elements(); ++i)
	{
		FEElement& el = dom.ElementRef(i);

		int NElementDataPoints;
		switch (elemMap->StorageFormat())
		{
		case Storage_Fmt::FMT_MULT:
			NElementDataPoints = el.Nodes();
			break;
		case Storage_Fmt::FMT_MATPOINTS:
			NElementDataPoints = el.GaussPoints();
			break;
		}

		for (int k = 0; k < NElementDataPoints; ++k)
		{
			vector<double>& vj = mappedData[n++];

			switch (dataType)
			{
			case FEDataType::FE_DOUBLE: elemMap->setValue(i, k, vj[0]); break;
			case FEDataType::FE_VEC3D:
			{
				vec3d v;
				v.x = vj[0];
				v.y = vj[1];
				v.z = vj[2];
				elemMap->setValue(i, k, v);
			}
			break;
			case FEDataType::FE_MAT3D:
			{
				mat3d v;
				v(0, 0) = vj[0]; v(0, 1) = vj[1]; v(0, 2) = vj[2];
				v(1, 0) = vj[3]; v(1, 1) = vj[4]; v(1, 2) = vj[5];
				v(2, 0) = vj[6]; v(2, 1) = vj[7]; v(2, 2) = vj[8];
				elemMap->setValue(i, k, v);
			}
			break;
			case FEDataType::FE_MAT3DS:
			{
				mat3ds v;
				v(0, 0) = vj[0];
				v(0, 1) = vj[1];
				v(0, 2) = vj[2];
				v(1, 1) = vj[3];
				v(1, 2) = vj[4];
				v(2, 2) = vj[5];
				elemMap->setValue(i, k, v);
			}
			break;
			default:
				assert(false);
			}
		}
	}

	return elemMap;

}

bool FEResetMesh::BuildMapData()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// make a copy of the old mesh
	// we need it for mapping data
	CopyMesh();

	// clear all map data
	ClearMapData();
	m_domainMapList.clear();
	m_domainMapList.resize(mesh.Domains());

	// only map domain data if requested
	if (m_bmap_data)
	{
		if (BuildDomainMapData() == false)
		{
			return false;
		}
	}

	// do the same thing for the user-defined mesh data
	return BuildUserMapData();
}

bool FEResetMesh::BuildDomainMapData()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// loop over all domains
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);

		feLog(" Processing domain: %s\n", dom.GetName().c_str());

		if (BuildDomainMapData(dom, i) == false)
		{
			return false;
		}
	}

	return true;
}

bool FEResetMesh::BuildDomainMapData(FEDomain& dom, int domIndex)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// write all material point data to a data stream
	DumpMemStream ar(fem);
	ar.Open(true, true);
	ar.WriteTypeInfo(true);

	// loop over all integration points
	int totalPoints = 0;
	for (int j = 0; j < dom.Elements(); ++j)
	{
		FEElement& el = dom.ElementRef(j);
		int nint = el.GaussPoints();
		for (int k = 0; k < nint; ++k)
		{
			FEMaterialPoint* mp = el.GetMaterialPoint(k);
			mp->Serialize(ar);
		}
		totalPoints += nint;
	}

	// figure out how much data was written for each material point
	size_t bytes = ar.bytesSerialized();

	size_t bytesPerPoint = bytes / totalPoints;
	assert((bytes%totalPoints) == 0);

	// re-open for reading
	ar.Open(false, true);

	// create an element set (we need this for the domain map below)
	FEDomain& oldDomain = m_meshCopy->Domain(domIndex);
	FEElementSet* elemSet = new FEElementSet(&fem);
	elemSet->Create(&oldDomain);

	// next, we need to figure out the datamaps for each data item
	vector<FEDomainMap*> mapList;
	DumpStream::DataBlock d;
	while (ar.bytesSerialized() < bytesPerPoint)
	{
		ar.readBlock(d);

		const char* typeStr = nullptr;
		FEDomainMap* map = nullptr;
		switch (d.dataType())
		{
		case TypeID::TYPE_DOUBLE: map = new FEDomainMap(FEDataType::FE_DOUBLE, Storage_Fmt::FMT_MATPOINTS); typeStr = "double"; break;
		case TypeID::TYPE_VEC3D: map = new FEDomainMap(FEDataType::FE_VEC3D, Storage_Fmt::FMT_MATPOINTS); typeStr = "vec3d"; break;
		case TypeID::TYPE_MAT3D: map = new FEDomainMap(FEDataType::FE_MAT3D, Storage_Fmt::FMT_MATPOINTS); typeStr = "mat3d"; break;
		case TypeID::TYPE_MAT3DS: map = new FEDomainMap(FEDataType::FE_MAT3DS, Storage_Fmt::FMT_MATPOINTS); typeStr = "mat3ds"; break;
		default:
			assert(false);
			throw std::runtime_error("Error in mapping data.");
		}

		map->Create(elemSet, 0.0);

		feLog("\tData map %d: %s\n", mapList.size(), typeStr);

		mapList.push_back(map);
	}

	feLog(" %d data maps identified.\n", mapList.size());

	// rewind for processing
	ar.Open(false, true);

	for (int j = 0; j < dom.Elements(); ++j)
	{
		FEElement& el = dom.ElementRef(j);
		int nint = el.GaussPoints();
		for (int k = 0; k < nint; ++k)
		{
			int m = 0;
			size_t bytesRead = 0;
			while (bytesRead < bytesPerPoint)
			{
				size_t size0 = ar.bytesSerialized();
				bool b = ar.readBlock(d); assert(b);
				size_t size1 = ar.bytesSerialized();
				bytesRead += size1 - size0;

				FEDomainMap* map = mapList[m];

				switch (d.dataType())
				{
				case TypeID::TYPE_DOUBLE: { double v = d.value<double>(); map->setValue(j, k, v); } break;
				case TypeID::TYPE_MAT3D: { mat3d  v = d.value<mat3d >(); map->setValue(j, k, v); } break;
				case TypeID::TYPE_MAT3DS: { mat3ds v = d.value<mat3ds>(); map->setValue(j, k, v); } break;
				}

				m++;
				assert(m <= mapList.size());
			}
		}
	}

	// Now, we need to project all the data onto the nodes
	for (int j = 0; j < mapList.size(); ++j)
	{
		feLog("\tProcessing data map %d ...", j);
		FEDomainMap* elemMap = mapList[j];

		FEDataType dataType = elemMap->DataType();
		FEDomainMap* nodeMap = new FEDomainMap(dataType, Storage_Fmt::FMT_NODE);

		FEElementSet* elset = const_cast<FEElementSet*>(elemMap->GetElementSet());
		nodeMap->Create(elset);

		bool bret = createNodeDataMap(dom, elemMap, nodeMap); assert(bret);
		m_domainMapList[domIndex].push_back(nodeMap);
		feLog("done.\n");
	}
    
    return true;
}

bool FEResetMesh::BuildUserMapData()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	m_userDataList.clear();
	int dataMaps = mesh.DataMaps();
	if (dataMaps > 0) feLog(" Processing user data maps:\n");

	for (int i = 0; i < dataMaps; ++i)
	{
		FEDataMap* dataMap = mesh.GetDataMap(i);

		// process domain map
		FEDomainMap* dmap = dynamic_cast<FEDomainMap*>(dataMap);
		if (dmap)
		{
			FEDataType dataType = dmap->DataType();
			if ((dataType != FEDataType::FE_DOUBLE) && 
				(dataType != FEDataType::FE_VEC3D)  &&
				(dataType != FEDataType::FE_MAT3D) ) return false;

			feLog("\tProcessing user data map \"%s\" ...", dmap->GetName().c_str());

			const FEElementSet* elset = dmap->GetElementSet();
			const FEDomainList& domainList = elset->GetDomainList();
			if (domainList.Domains() != 1) return false;
			FEDomain& dom = const_cast<FEDomain&>(*domainList.GetDomain(0));

			FEElementSet* oldElemSet = m_meshCopy->FindElementSet(dom.GetName()); assert(oldElemSet);

			FEDomainMap* nodeMap = new FEDomainMap(dataType, Storage_Fmt::FMT_NODE);
			nodeMap->Create(oldElemSet);

			// create a node data map of this domain map
			createNodeDataMap(dom, dmap, nodeMap);
			m_userDataList.push_back(nodeMap);

			feLog("done.\n");
		}

		FESurfaceMap* smap = dynamic_cast<FESurfaceMap*>(dataMap);
		if (smap)
		{
			assert(false);
		}
	}

	return true;
}

// Transfer data to new mesh
void FEResetMesh::TransferMapData()
{
	// transfer domain data
	if (m_bmap_data) TransferDomainMapData();

	// transfer user-defined maps
	TransferUserMapData();
}

// Transfer domain data back to the new mesh
void FEResetMesh::TransferDomainMapData()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// loop over all domains
	if (m_bmap_data && m_domainMapList.size())
	{
		for (int i = 0; i < mesh.Domains(); ++i)
		{
			FEDomain& dom = mesh.Domain(i);

			std::vector<FEDomainMap*>& nodeMap_i = m_domainMapList[i];
			int mapCount = nodeMap_i.size();

			if (mapCount > 0)
			{
				feLog(" Mapping data for domain \"%s\":\n", dom.GetName().c_str());

				// we need an element set for the domain maps below
				FEElementSet* elemSet = new FEElementSet(&fem);
				elemSet->Create(&dom);

				// build source point list
				FEDomain& oldDomain = m_meshCopy->Domain(i);
				vector<vec3d> srcPoints; srcPoints.reserve(oldDomain.Nodes());
				for (int j = 0; j < oldDomain.Nodes(); ++j)
				{
					vec3d r = m_bmap_current ? oldDomain.Node(j).m_rt : oldDomain.Node(j).m_r0;
					srcPoints.push_back(r);
				}

				// build target node list
				vector<vec3d> trgPoints; trgPoints.reserve(dom.Elements());
				for (int j = 0; j < dom.Elements(); ++j)
				{
					FEElement& el = dom.ElementRef(j);
					int nint = el.GaussPoints();
					for (int k = 0; k < nint; ++k)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(k);
						vec3d r = m_bmap_current ? mp.m_rt : mp.m_r0;
						trgPoints.push_back(r);
					}
				}

				// set up mapper
				FEMeshDataInterpolator* mapper = nullptr;
				switch (m_transferMethod)
				{
				case TRANSFER_SHAPE:
				{
					FEDomain* oldDomain = &m_meshCopy->Domain(i);
					FEDomainShapeInterpolator* dsm = new FEDomainShapeInterpolator(oldDomain);
					dsm->SetTargetPoints(trgPoints);
					mapper = dsm;
				}
				break;
				case TRANSFER_MLQ:
				{
					FELeastSquaresInterpolator* MLQ = new FELeastSquaresInterpolator;
					MLQ->SetNearestNeighborCount(m_nnc);
					MLQ->SetDimension(m_nsdim);
					MLQ->SetSourcePoints(srcPoints);
					MLQ->SetTargetPoints(trgPoints);
					mapper = MLQ;
				}
				break;
				default:
					assert(false);
					return;
				}
				if (mapper->Init() == false)
				{
					assert(false);
					throw std::runtime_error("Failed to initialize LLQ");
				}

				// loop over all the domain maps
				vector<FEDomainMap*> elemMapList(mapCount);
				for (int j = 0; j < mapCount; ++j)
				{
					feLog("\tMapping map %d ...", j);
					FEDomainMap* nodeMap = nodeMap_i[j];

					// map node data to integration points
					FEDomainMap* elemMap = NodeToElemData(fem, dom, nodeMap, mapper, Storage_Fmt::FMT_MATPOINTS);

					elemMapList[j] = elemMap;
					feLog("done.\n");
				}

				// now we need to reconstruct the data stream
				DumpMemStream ar(fem);
				ar.Open(true, true);

				for (int j = 0; j < dom.Elements(); ++j)
				{
					FEElement& el = dom.ElementRef(j);
					int nint = el.GaussPoints();

					for (int k = 0; k < nint; ++k)
					{
						for (int l = 0; l < mapCount; ++l)
						{
							FEDomainMap* map = elemMapList[l];
							switch (map->DataType())
							{
							case FEDataType::FE_DOUBLE: { double v = map->value<double>(j, k); ar << v; } break;
							case FEDataType::FE_VEC3D: { vec3d  v = map->value<vec3d >(j, k); ar << v; } break;
							case FEDataType::FE_MAT3D: { mat3d  v = map->value<mat3d >(j, k); ar << v; } break;
							case FEDataType::FE_MAT3DS: { mat3ds v = map->value<mat3ds>(j, k); ar << v; } break;
							default:
								assert(false);
							}
						}
					}
				}

				// time to serialize everything back to the new integration points
				ar.Open(false, true);
				for (int j = 0; j < dom.Elements(); ++j)
				{
					FEElement& el = dom.ElementRef(j);
					int nint = el.GaussPoints();

					for (int k = 0; k < nint; ++k)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(k);
						mp.Serialize(ar);
					}
				}
			}
		}
	}
}

// Transfer user data maps to new mesh
void FEResetMesh::TransferUserMapData()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// map mesh data
	for (int i = 0; i < m_userDataList.size(); ++i)
	{
		FEDomainMap* elemMap = dynamic_cast<FEDomainMap*>(mesh.GetDataMap(i));
		FEDomainMap* nodeMap = m_userDataList[i];

		Storage_Fmt dataFormat = static_cast<Storage_Fmt>(elemMap->StorageFormat()); // get the storage format
		string name = elemMap->GetName(); // get name of elemMap

		// set target data format
		Storage_Fmt targetFormat;
		switch (dataFormat) {
		case Storage_Fmt::FMT_ITEM:
			feLog("\tChanging storage format of user element map \"%s\" to MATPOINTS\n", name.c_str());
			targetFormat = Storage_Fmt::FMT_MATPOINTS;
			break;
		case Storage_Fmt::FMT_MULT:
			targetFormat = dataFormat;
			break;
		case Storage_Fmt::FMT_MATPOINTS:
			targetFormat = dataFormat;
			break;
		default:
			assert(false);
			throw std::runtime_error("Error in FEResetMesh::TransferUserMapData");
			break;
		}

		feLog("\tMapping user map \"%s\" ...\n", name.c_str());

		// build the source point list
		vector<vec3d> srcPoints;
		const FEElementSet* oldSet = nodeMap->GetElementSet();
		FEMesh* oldMesh = oldSet->GetMesh(); assert(oldMesh == m_meshCopy);
		FENodeList oldNodeList = oldSet->GetNodeList();
		srcPoints.resize(oldNodeList.Size());
		for (int i = 0; i < oldNodeList.Size(); ++i)
		{
			FENode& oldNode_i = oldMesh->Node(oldNodeList[i]);
			vec3d r = m_bmap_current ? oldNode_i.m_rt : oldNode_i.m_r0;
			srcPoints[i] = r;
		}

		// build target points.
		const FEDomainList& domainList = elemMap->GetElementSet()->GetDomainList();
		FEDomain& dom = const_cast<FEDomain&>(*domainList.GetDomain(0));
		int NE = dom.Elements();
		vector<vec3d> trgPoints; trgPoints.reserve(NE);
		switch (targetFormat) {
		case Storage_Fmt::FMT_MATPOINTS:
			for (int n = 0; n < NE; ++n)
			{
				FEElement& el = dom.ElementRef(n);
				int nint = el.GaussPoints();

				// get element node coordinates
				vec3d r[FEElement::MAX_NODES];
				if (m_bmap_current) {
					for (int i = 0; i < el.Nodes(); ++i) r[i] = mesh.Node(el.m_node[i]).m_rt;  // this is not updated at material point
				}
				else {
					for (int i = 0; i < el.Nodes(); ++i) r[i] = mesh.Node(el.m_node[i]).m_r0;  // this could be taken from material point
				}

				// store material point coordinate
				for (int l = 0; l < nint; ++l) trgPoints.push_back(el.Evaluate(r, l));

			}
			break;
		case Storage_Fmt::FMT_MULT:

			for (int n = 0; n < NE; ++n)
			{
				FEElement& el = dom.ElementRef(n);
				int ne = el.Nodes();
				for (int l = 0; l < ne; ++l)
				{
					vec3d r = m_bmap_current ? mesh.Node(el.m_node[l]).m_rt : mesh.Node(el.m_node[l]).m_r0;
					trgPoints.push_back(r);
				}
			}
			break;
		}

		// set up mapper
		FEMeshDataInterpolator* mapper = nullptr;
		switch (m_transferMethod)
		{
		case TRANSFER_SHAPE:
		{
			FEDomain* oldDomain = m_meshCopy->FindDomain(dom.GetName());
			FEDomainShapeInterpolator* dsm = new FEDomainShapeInterpolator(oldDomain, m_bmap_current, m_atol);
			dsm->SetTargetPoints(trgPoints);
			mapper = dsm;
		}
		break;
		case TRANSFER_MLQ:
		{
			FELeastSquaresInterpolator* MLQ = new FELeastSquaresInterpolator;
			MLQ->SetNearestNeighborCount(m_nnc);
			MLQ->SetDimension(m_nsdim);
			MLQ->SetSourcePoints(srcPoints);
			MLQ->SetTargetPoints(trgPoints);
			mapper = MLQ;
		}
		break;
		default:
			assert(false);
			return;
		}
		if (mapper->Init() == false)
		{
			assert(false);
			throw std::runtime_error("Failed to initialize LLQ");
		}

		// Map node data to elements
		// We can interpolate from points to points, so target can either be nodal points or gauss points, not "element"		
		// TODO: does anything point to this map or is it looked up by name? does it matter if we change storage?
		//delete elemMap; // delete existing
		FEDomainMap* updatedMap = NodeToElemData(fem, dom, nodeMap, mapper, targetFormat);
		updatedMap->SetName(name); // set same name
		*elemMap = *updatedMap; // copy

		// clean up
		delete updatedMap;
		delete mapper;

		feLog("\tDone mapping user map \"%s\".\n", name.c_str());
	}
}

void FEResetMesh::ClearGradientMap() {

	delete m_FmapCopy;

}

bool FEResetMesh::RestoreGradientMap() {

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get custom data map used to store deformation gradient
	FEDomainMap* map = dynamic_cast<FEDomainMap*>(mesh.FindDataMap(m_Fmap_name)); // returns a copy of the pointer, so if it is changed it is not stored :(
	if (!map) {
		return false;
	}

	// restore (copy operator)
	*map = *m_FmapCopy;

	return true;
}

bool FEResetMesh::UpdateGradientMap()
{

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get custom data map used to store deformation gradient
	FEDomainMap* map = dynamic_cast<FEDomainMap*>(mesh.FindDataMap(m_Fmap_name)); // returns a copy of the pointer, so if it is changed it is not stored :(
	if (!map) {
		return false;
	}

	// create backup
	m_FmapCopy = new FEDomainMap(*map);		

	// 
	const FEElementSet* set = map->GetElementSet();

	FEDataType dataType = map->DataType();
	if (dataType != FE_MAT3D) return false;

	int storageFormat = map->StorageFormat();
	if (storageFormat != FMT_MATPOINTS) return false;

	// get the domain of the set, works only for a single domain
	const FEDomainList& domainList = set->GetDomainList();
	if (domainList.Domains() != 1) return false;
	FESolidDomain& dom = const_cast<FESolidDomain&>(dynamic_cast<const FESolidDomain&>(*domainList.GetDomain(0)));
	// developers should make methods const whenever possible, so const_cast is not neccesary
	// TODO check valid domain

	// loop all elements compute deformation gradients in integration points and multiply existing prestrain gradient, then store in map
	int N = set->Elements();
	for (int i = 0; i < N; ++i)
	{
		FESolidElement* pel = dynamic_cast<FESolidElement*>(mesh.FindElementFromID((*set)[i]));
		if (pel == nullptr) return false;
		FESolidElement& el = *pel;

		int ni = el.GaussPoints();
		for (int j = 0; j < ni; ++j)
		{
			FEPrestrainMaterialPoint* pp = dynamic_cast<FEPrestrainMaterialPoint*>(el.GetMaterialPoint(j));
			assert(pp);
			mat3d F0 = pp->initialPrestrain();  // initial prestrain, do not include correction
			mat3d F;
			dom.defgrad(el, F, j);  // compute deformation gradient from current mesh displacement relative to reference
			map->setValue(i, j, F * F0);  // product, notice if F0=I, then after remesh 1 you get F*F0=F as intended
		}
	}
	return true;
}

bool FEResetMesh::UpdateDisplacementMap()
{

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// check dimensions
	assert(m_nodeDispl.size() == mesh.Nodes());

	// get custom data map used to store deformation gradient
	FEDomainMap* map = dynamic_cast<FEDomainMap*>(mesh.FindDataMap(m_Dmap_name)); // returns a copy of the pointer, so if it is changed it is not stored :(
	if (!map) {
		return false;
	}

	// get elements
	const FEElementSet* set = map->GetElementSet();

	FEDataType dataType = map->DataType();
	if (dataType != FE_VEC3D) return false;

	int storageFormat = map->StorageFormat();
	if (storageFormat != FMT_MULT) return false;

	// get the domain of the set, works only for a single domain
	const FEDomainList& domainList = set->GetDomainList();
	if (domainList.Domains() != 1) return false;
	FESolidDomain& dom = const_cast<FESolidDomain&>(dynamic_cast<const FESolidDomain&>(*domainList.GetDomain(0)));
	// developers should make methods const whenever possible, so const_cast is not neccesary
	// TODO check valid domain

	// loop all elements get nodal displacements
	int N = set->Elements();
	for (int i = 0; i < N; ++i)
	{
		FESolidElement* pel = dynamic_cast<FESolidElement*>(mesh.FindElementFromID((*set)[i]));
		if (pel == nullptr) return false;
		FESolidElement& el = *pel;

		int neln = el.Nodes();
		for (int j = 0; j < neln; ++j)
		{
			vec3d u = map->value<vec3d>(i, j); // existing value
			map->setValue(i, j, u + m_nodeDispl[el.m_node[j]]);  // increment
		}
	}
	return true;
}

bool FEResetMesh::GradientsFromDisplacementMap()
{

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// get custom data map used to store zeroed displacements
	FEDomainMap* Dmap = dynamic_cast<FEDomainMap*>(mesh.FindDataMap(m_Dmap_name)); // returns a copy of the pointer, so if it is changed it is not stored :(
	if (!Dmap) {
		return false;
	}

	// get custom data map used to deformation gradients
	FEDomainMap* Fmap = dynamic_cast<FEDomainMap*>(mesh.FindDataMap(m_Fmap_name)); // returns a copy of the pointer, so if it is changed it is not stored :(
	if (!Fmap) {
		return false;
	}

	// generate gradients

	// sets should be the same for Dmap and Fmap, TODO check
	const FEElementSet* Fset = Fmap->GetElementSet();
	const FEElementSet* Dset = Dmap->GetElementSet();

	FEDataType dataType = Fmap->DataType();
	if (dataType != FE_MAT3D) return false;

	int storageFormat = Fmap->StorageFormat();
	if (storageFormat != FMT_MATPOINTS) return false;

	double xmin = 1e37, xmax = -1e37;
	double ymin = 1e37, ymax = -1e37;
	double zmin = 1e37, zmax = -1e37;

	int N = Fset->Elements();
	for (int i = 0; i < N; ++i)
	{

		assert((*Dset)[i] == (*Fset)[i]);

		FESolidElement* pel = dynamic_cast<FESolidElement*>(mesh.FindElementFromID((*Fset)[i]));
		if (pel == nullptr) return false;
		FESolidElement& el = *pel;

		int ne = el.Nodes();
		vector<vec3d> u(ne), X(ne);
		for (int j = 0; j < ne; j++)
		{
			u[j] = Dmap->value<vec3d>(i, j);

			if (u[j].x < xmin) xmin = u[j].x;
			if (u[j].y < ymin) ymin = u[j].y;
			if (u[j].z < zmin) zmin = u[j].z;
			if (u[j].x > xmax) xmax = u[j].x;
			if (u[j].y > ymax) ymax = u[j].y;
			if (u[j].z > zmax) zmax = u[j].z;

			X[j] = mesh.Node(el.m_node[j]).m_r0 - u[j];
		}

		int ni = el.GaussPoints();
		for (int j = 0; j < ni; ++j)
		{
			mat3d F;
			defgrad(el, X, u, F, j);
			Fmap->setValue(i, j, F);
		}
	}

	printf("displacement ranges\n  x: %f,%f\n  y: %f,%f\n  z: %f,%f\n", xmin, xmax, ymin, ymax, zmin, zmax);

	return true;
}

double FEResetMesh::invjac0(FESolidElement& el, std::vector<vec3d>& r0, double Ji[3][3], int n)
{
	// calculate Jacobian
	double J[3][3] = { 0 };
	int neln = el.Nodes();
	for (int i = 0; i < neln; ++i)
	{
		const double& Gri = el.Gr(n)[i];
		const double& Gsi = el.Gs(n)[i];
		const double& Gti = el.Gt(n)[i];

		const double& x = r0[i].x;
		const double& y = r0[i].y;
		const double& z = r0[i].z;

		J[0][0] += Gri * x; J[0][1] += Gsi * x; J[0][2] += Gti * x;
		J[1][0] += Gri * y; J[1][1] += Gsi * y; J[1][2] += Gti * y;
		J[2][0] += Gri * z; J[2][1] += Gsi * z; J[2][2] += Gti * z;
	}

	// calculate the determinant
	double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
		+ J[0][1] * (J[1][2] * J[2][0] - J[2][2] * J[1][0])
		+ J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

	if (det != 0.0)
	{
		// calculate the inverse jacobian
		double deti = 1.0 / det;

		Ji[0][0] = deti * (J[1][1] * J[2][2] - J[1][2] * J[2][1]);
		Ji[1][0] = deti * (J[1][2] * J[2][0] - J[1][0] * J[2][2]);
		Ji[2][0] = deti * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

		Ji[0][1] = deti * (J[0][2] * J[2][1] - J[0][1] * J[2][2]);
		Ji[1][1] = deti * (J[0][0] * J[2][2] - J[0][2] * J[2][0]);
		Ji[2][1] = deti * (J[0][1] * J[2][0] - J[0][0] * J[2][1]);

		Ji[0][2] = deti * (J[0][1] * J[1][2] - J[1][1] * J[0][2]);
		Ji[1][2] = deti * (J[0][2] * J[1][0] - J[0][0] * J[1][2]);
		Ji[2][2] = deti * (J[0][0] * J[1][1] - J[0][1] * J[1][0]);
	}

	return det;
}

double FEResetMesh::defgrad(FESolidElement& el, std::vector<vec3d>& X, std::vector<vec3d>& u, mat3d& F, int n)
{
	// calculate inverse jacobian
	double Ji[3][3];
	invjac0(el, X, Ji, n);

	// shape function derivatives
	double* Grn = el.Gr(n);
	double* Gsn = el.Gs(n);
	double* Gtn = el.Gt(n);

	// calculate deformation gradient
	F[0][0] = F[0][1] = F[0][2] = 0;
	F[1][0] = F[1][1] = F[1][2] = 0;
	F[2][0] = F[2][1] = F[2][2] = 0;
	int neln = el.Nodes();
	for (int i = 0; i < neln; ++i)
	{
		double Gri = Grn[i];
		double Gsi = Gsn[i];
		double Gti = Gtn[i];

		double x = u[i].x;
		double y = u[i].y;
		double z = u[i].z;

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		double GX = Ji[0][0] * Gri + Ji[1][0] * Gsi + Ji[2][0] * Gti;
		double GY = Ji[0][1] * Gri + Ji[1][1] * Gsi + Ji[2][1] * Gti;
		double GZ = Ji[0][2] * Gri + Ji[1][2] * Gsi + Ji[2][2] * Gti;

		// calculate deformation gradient F
		F[0][0] += GX * x; F[0][1] += GY * x; F[0][2] += GZ * x;
		F[1][0] += GX * y; F[1][1] += GY * y; F[1][2] += GZ * y;
		F[2][0] += GX * z; F[2][1] += GY * z; F[2][2] += GZ * z;
	}

	F[0][0] += 1.0;
	F[1][1] += 1.0;
	F[2][2] += 1.0;

	double D = F.det();

	return D;
}