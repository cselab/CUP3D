/*
 *  LayerToVTK.h
 *  Fluids2D
 *
 *  Created by Diego Rossinelli on 7/21/09.
 *  Copyright 2009 CSE Lab, ETH Zurich. All rights reserved.
 *
 */

#pragma once

#include <vtkPoints.h>
#include <vtkCell.h>
//#include <vtkUnstructuredGrid.h>
#include <vtkImageData.h>
#include <vtkImageNoiseSource.h>
#include <vtkFloatArray.h>
//#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

inline void dumpLayer2VTK(int istep, string sFileNamePattern, const Layer& scalar_field, const int nLayers)
{
	const int sizeX = scalar_field.sizeX;
	const int sizeY = scalar_field.sizeY;
	char buf[500];
	
	sprintf(buf, sFileNamePattern.data(), istep);
	string sFileName(buf);
	
	vtkImageData * grid = vtkImageData::New();
	
	if (nLayers==1)
	{
		grid->SetExtent(0,sizeX-1,0,sizeY-1,0,0);
		grid->SetDimensions(sizeX, sizeY, 1);
#ifndef _VTK62_
		grid->SetScalarTypeToFloat();
		grid->SetNumberOfScalarComponents(1);
		grid->AllocateScalars();
#else
		grid->AllocateScalars(VTK_FLOAT,1);
#endif
		grid->SetSpacing(1./sizeX, 1./sizeX, 1);
		grid->SetOrigin(0,0,0);
		
		for(int iy=0; iy<sizeY; iy++)
			for(int ix=0; ix<sizeX; ix++)
				grid->SetScalarComponentFromFloat(ix, iy, 0, 0, scalar_field.read(ix, iy));
	}
	else
	{
		grid->SetExtent(0,sizeX-1,0,sizeY-1,0,0);
		grid->SetDimensions(sizeX, sizeY, 2);
#ifndef _VTK62_
		grid->SetScalarTypeToFloat();
		grid->SetNumberOfScalarComponents(1);
		grid->AllocateScalars();
#else
		grid->AllocateScalars(VTK_FLOAT,1);
#endif
		grid->SetSpacing(1./sizeX, 1./sizeX, 1);
		grid->SetOrigin(0,0,0);
		
		for(int iy=0; iy<sizeY; iy++)
			for(int ix=0; ix<sizeX; ix++)
			{
				grid->SetScalarComponentFromFloat(ix, iy, 0, 0, scalar_field.read(ix, iy, 0));
				grid->SetScalarComponentFromFloat(ix, iy, 1, 0, scalar_field.read(ix, iy, 1));
			}
	}
	
	vtkXMLImageDataWriter * writer = vtkXMLImageDataWriter::New();
	writer->SetFileName(sFileName.c_str());
#ifndef _VTK62_
	writer->SetInput(grid);
#else
	writer->SetInputData(grid);
#endif
	writer->Write();
	
	writer->Delete();
	grid->Delete();
	
	//	printf("done with dumping %s.\n", buf);
}

