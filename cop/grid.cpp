//=======================================================================
// Regional Earth System Model (RegESM)
// Copyright (c) 2013-2019 Ufuk Turuncoglu
// Licensed under the MIT License.
//=======================================================================
#include "grid.h"

#include <iostream>

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkMath.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkIntArray.h"
#include "vtkPointData.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCPPythonAdaptorAPI.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkNew.h"

namespace ESMFAdaptor {
  // constructor
  template<GridType gridType>
  Grid<gridType>::Grid() :
    GridName(NULL),
    NProc(-1), 
    MPIRank(-1), 
    NPoints(0),
    NCells(0),
    NLon(0), 
    NLat(0), 
    NLev(0),
    IndLonL(0), 
    IndLonU(0), 
    IndLatL(0), 
    IndLatU(0), 
    IndLevL(0),
    IndLevU(0),
    Lon(NULL),
    Lat(NULL),
    Lev(NULL),
    grid(NULL) {}

  //
  // deconstructor
  //
  template<GridType gridType>
  Grid<gridType>::~Grid(){}

  //
  // creates the grid and the data used to add attributes to the grid
  //
  template<GridType gridType>
  vtkSmartPointer<vtkMultiBlockDataSet> Grid<gridType>::CreateGrid(const char* name){
    //
    // create structured grid
    // 
    vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();

    //
    // create data structure (points) to store grid coordinates and assign them to grid 
    //
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    grid->SetPoints(points);

    //
    // insert grid coordinate to points
    //
    points->SetDataTypeToFloat();
    points->SetNumberOfPoints(this->NPoints);
    if (this->Lev == NULL) { // 2d
      for (int i = 0; i < this->NPoints; i++) {
        points->SetPoint(i, this->Lon[i], this->Lat[i], 0);
      }
    } else { // 3d
      for (int i = 0; i < this->NPoints; i++) {
        points->SetPoint(i, this->Lon[i], this->Lat[i], this->Lev[i]);
      }
    }
    grid->SetPoints(points);
    if (this->Lev == NULL) { // 2d
      grid->SetExtent(this->IndLonL-1, this->IndLonU-1, 
                      this->IndLatL-1, this->IndLatU-1, 
                      0, 0);
    } else { // 3d
      grid->SetExtent(this->IndLonL-1, this->IndLonU-1, 
                      this->IndLatL-1, this->IndLatU-1, 
                      this->IndLevL-1, this->IndLevU-1);
    }
    
    //
    // create multi-block grid 
    //
    vtkMultiBlockDataSet* mb = vtkMultiBlockDataSet::New();
    vtkNew<vtkMultiPieceDataSet> mpds;
    mpds->SetNumberOfPieces(this->NProc);
    mpds->SetPiece(this->MPIRank, grid.GetPointer());
    mb->SetNumberOfBlocks(1);
    mb->SetBlock(0, mpds.GetPointer());

    return mb;
  }

  template<GridType gridType>
  void Grid<gridType>::SetAttributeValue(vtkCPDataDescription* coprocessorData, float* data, const char* vname, const char* pname, int* size, int* mpiSize, int* mpiRank) {
    //
    // Get grid
    //
    vtkCPInputDataDescription* idd = coprocessorData->GetInputDescriptionByName(pname);
    vtkMultiBlockDataSet *grid = vtkMultiBlockDataSet::SafeDownCast(idd->GetGrid());
    if (!grid) {
      vtkGenericWarningMacro("No adaptor grid to attach field data to.");
      return;
    }
   
    //
    // Create dataset and fill with input data
    //
    vtkMultiPieceDataSet *multiPiece = vtkMultiPieceDataSet::SafeDownCast(grid->GetBlock(0));
    vtkDataSet *dataSet = vtkDataSet::SafeDownCast(multiPiece->GetPiece(*mpiRank));

    vtkSmartPointer<vtkFloatArray> field = vtkSmartPointer<vtkFloatArray>::New();
    field->SetName(vname);
    field->SetNumberOfComponents(1);
      
    field->SetNumberOfValues(*size);
    for (int i = 0; i < *size; i++) {
      field->SetValue(i, data[i]);
    } 

    int i = dataSet->GetPointData()->AddArray(field);
  }

  // Creates a 2D and a 3D grid of the specified 'type'
  template<GridType gridType>
  void Grid<gridType>::Create(const char* name) {
    this->grid = CreateGrid(name);
  }

  template<GridType gridType>
  void Grid<gridType>::SetName(const char* name) {
    char* buf = new char[strlen(name)];
    strcpy(buf, name);    
    this->GridName = buf;
  }

  template<GridType gridType>
  void Grid<gridType>::SetNProc(int nproc) {
    this->NProc = nproc;
  }

  template<GridType gridType>
  void Grid<gridType>::SetMPIRank(int rank) {
    this->MPIRank = rank;
  }

  template<GridType gridType>  
  void Grid<gridType>::SetNCells(int ncells) {
    this->NCells = ncells;
  }

  template<GridType gridType>
  void Grid<gridType>::SetNPoints(int npoints) {
    this->NPoints = npoints;
  }

  template<GridType gridType>  
  void Grid<gridType>::SetLon(int nlon, int npoints, float* lon) {
    this->NLon = nlon;
    if (lon) {
      this->Lon = new float[npoints];
      std::copy(lon, lon + npoints, this->Lon);
    }
  }
  
  template<GridType gridType>  
  void Grid<gridType>::SetLat(int nlat, int npoints, float* lat) { 
    this->NLat = nlat;
    if (lat) {
      this->Lat = new float[npoints];
      std::copy(lat, lat + npoints, this->Lat);
    }
  }

  template<GridType gridType>
  void Grid<gridType>::SetLev(int nlev, int npoints, float* lev) {
    this->NLev = nlev;
    if (lev) {
      this->Lev = new float[npoints];
      std::copy(lev, lev + npoints, this->Lev);
    }
  }

  template<GridType gridType>
  void Grid<gridType>::SetBounds(int* lb, int* ub) {
    if (lb) {
      this->IndLonL = lb[0];
      this->IndLatL = lb[1];
      this->IndLevL = 0;
      if (lb[2] > 0) this->IndLevL = lb[2]; 
    }
    if (ub) {
      this->IndLonU = ub[0];
      this->IndLatU = ub[1];
      this->IndLevU = 0;
      if (ub[2] > 0) this->IndLevU = ub[2]; 
    }
  }

  template<GridType gridType>
  bool Grid<gridType>::SetToCoprocessor(vtkCPDataDescription* coprocessorData, const char* name, int* dims, vtkSmartPointer<vtkMultiBlockDataSet> grid) {
    vtkCPInputDataDescription* idd = coprocessorData->GetInputDescriptionByName(name);
    if (idd) {
      if (dims[2] == 0) {
        idd->SetWholeExtent(0, dims[0]-1, 0, dims[1]-1, 0, 0);
      } else {
        idd->SetWholeExtent(0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1);
      }
      idd->SetGrid(grid);
      return true;
    }
    return false;
  }

  // returns the grid
  template<GridType gridType>
  vtkSmartPointer<vtkMultiBlockDataSet> Grid<gridType>::GetGrid() const {
    return this->grid;
  }

  // returns the grid
  template<GridType gridType>
  char* Grid<gridType>::GetName() {
    return this->GridName;
  }

  template<GridType gridType>
  int Grid<gridType>::GetNPoints() {
    return this->NPoints;
  }

  template<GridType gridType>
  int Grid<gridType>::GetNCells() {
    return this->NCells;
  }

  // instantiations for the template classes
  template class Grid<RECTILINEAR>;
};
