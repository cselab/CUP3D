import os
import re

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

#rootDir = '/scratch/daint/cconti/Samara_512_SP_160316_Dora_Euler_DLM1_lowDump/'
#rootDir = '/scratch/daint/cconti/SamaraSP_512_RK2_220316_Dora_lowDump/'
rootDir = '/scratch/daint/cconti/SamaraSP_512_RK2_230316_Dora_DumpViz/'

for dirName, subDirList, fileList in os.walk(rootDir):
	print(dirName)
	for f in fileList:
		if f.endswith('.xmf'):
			#frame = re.search('Samara_512_SP_160316_Dora_Euler_DLM1_lowDump-(.+?).xmf', f)
			#frame = re.search('SamaraSP_512_RK2_220316_Dora_lowDump-(.+?).xmf', f)
			frame = re.search('SamaraSP_512_RK2_230316_Dora_DumpViz-(.+?).xmf', f)
			#pathChi = dirName+'/Chi'+frame.group(1)+'.ply'
			pathVort = dirName+'/Vorticity'+frame.group(1)+'.ply'
			pathVortB = dirName+'/VorticityB'+frame.group(1)+'.ply'
			pathVortC = dirName+'/VorticityC'+frame.group(1)+'.ply'

			if not os.path.isfile(pathVortC):
				# create a new 'XDMF Reader'
				file=dirName+f
				print(file)
				filename = XDMFReader(FileNames=file)

				# create a new 'XDMF Reader'
				filename.PointArrayStatus = ['data']

				# Properties modified on filename
				filename.GridStatus = ['Grid_2']

				# get active view
				renderView1 = GetActiveViewOrCreate('RenderView')
				# uncomment following to set a specific view size
				# renderView1.ViewSize = [2019, 1098]

				# show data in view
				filenameDisplay = Show(filename, renderView1)
				# trace defaults for the display properties.
				filenameDisplay.Representation = 'Outline'
				filenameDisplay.ColorArrayName = [None, '']
				filenameDisplay.GlyphType = 'Arrow'
				filenameDisplay.ScalarOpacityUnitDistance = 0.0016914558667664823
				filenameDisplay.Slice = 511

				# reset view to fit data
				renderView1.ResetCamera()

				## create a new 'Calculator'
				#calculator1 = Calculator(Input=filename)
				#calculator1.Function = ''

				## Properties modified on calculator1
				#calculator1.Function = 'data_4'

				## show data in view
				#calculator1Display = Show(calculator1, renderView1)
				## trace defaults for the display properties.
				#calculator1Display.Representation = 'Outline'
				#calculator1Display.ColorArrayName = ['POINTS', '']
				#calculator1Display.GlyphType = 'Arrow'
				#calculator1Display.ScalarOpacityUnitDistance = 0.0016914558667664823
				#calculator1Display.Slice = 511

				# hide data in view
				#Hide(filename, renderView1)

				## create a new 'Contour'
				#contour1 = Contour(Input=calculator1)
				#contour1.ContourBy = ['POINTS', 'Result']
				#contour1.Isosurfaces = [0.5]
				#contour1.PointMergeMethod = 'Uniform Binning'

				## show data in view
				#contour1Display = Show(contour1, renderView1)
				## trace defaults for the display properties.
				#contour1Display.ColorArrayName = [None, '']
				#contour1Display.GlyphType = 'Arrow'
				#contour1Display.SetScaleArray = [None, '']
				#contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
				#contour1Display.OpacityArray = [None, '']
				#contour1Display.OpacityTransferFunction = 'PiecewiseFunction'

				## set active source
				#SetActiveSource(filename)

				## set active source
				#SetActiveSource(contour1)

				## save data
				#SaveData(pathChi, proxy=contour1)
				#print(pathChi)
				
			#	# set active source
				SetActiveSource(filename)
				
				# create a new 'Calculator'
				calculator2 = Calculator(Input=filename)
				calculator2.Function = 'data_6'
				
#				# Properties modified on calculator2
#				calculator2.Function = 'data_1*iHat+data_2*jHat+data_3*kHat'
#				
#				# show data in view
#				calculator2Display = Show(calculator2, renderView1)
#				# trace defaults for the display properties.
#				calculator2Display.Representation = 'Outline'
#				calculator2Display.ColorArrayName = [None, '']
#				calculator2Display.GlyphType = 'Arrow'
#				calculator2Display.ScalarOpacityUnitDistance = 0.0033096687538939264
#				calculator2Display.Slice = 511
#				
#				# hide data in view
#				Hide(filename, renderView1)
#				
#				# create a new 'Compute Derivatives'
#				computeDerivatives1 = ComputeDerivatives(Input=calculator2)
#				computeDerivatives1.Scalars = [None, '']
#				computeDerivatives1.Vectors = ['POINTS', 'Result']
#				
#				# Properties modified on computeDerivatives1
#				computeDerivatives1.OutputVectorType = 'Vorticity'
#				
#				# show data in view
#				computeDerivatives1Display = Show(computeDerivatives1, renderView1)
#				# trace defaults for the display properties.
#				computeDerivatives1Display.Representation = 'Outline'
#				computeDerivatives1Display.ColorArrayName = [None, '']
#				computeDerivatives1Display.GlyphType = 'Arrow'
#				computeDerivatives1Display.ScalarOpacityUnitDistance = 0.0033096687538939264
#				computeDerivatives1Display.Slice = 511
#
#				# create a new 'Cell Data to Point Data'
#				cellDatatoPointData1 = CellDatatoPointData(Input=computeDerivatives1)
#
#				# show data in view
#				cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1)
#				# trace defaults for the display properties.
#				cellDatatoPointData1Display.Representation = 'Outline'
#				cellDatatoPointData1Display.ColorArrayName = [None, '']
#				cellDatatoPointData1Display.GlyphType = 'Arrow'
#				cellDatatoPointData1Display.ScalarOpacityUnitDistance = 0.0033829117335329646
#				cellDatatoPointData1Display.Slice = 511
#				
#				# hide data in view
#				Hide(calculator2, renderView1)
#				
#				# create a new 'Calculator'
#				calculator3 = Calculator(Input=cellDatatoPointData1)
#				calculator3.Function = 'mag(Vorticity)'
#				
#				# show data in view
#				calculator3Display = Show(calculator3, renderView1)
#				# trace defaults for the display properties.
#				calculator3Display.Representation = 'Outline'
#				calculator3Display.ColorArrayName = ['POINTS', '']
#				calculator3Display.GlyphType = 'Arrow'
#				calculator3Display.ScalarOpacityUnitDistance = 0.0033096687538939264
#				calculator3Display.Slice = 511
#				
#				# hide data in view
#				Hide(computeDerivatives1, renderView1)
#				
				# create a new 'Contour'
#				contour2 = Contour(Input=calculator3)
				contour2 = Contour(Input=calculator2)
				contour2.ContourBy = ['POINTS', 'Result']
				contour2.Isosurfaces = [50.0]
				contour2.PointMergeMethod = 'Uniform Binning'
				
				# show data in view
				contour3Display = Show(contour1, renderView1)
				# trace defaults for the display properties.
				contour3Display.ColorArrayName = [None, '']
				contour3Display.GlyphType = 'Arrow'
				contour3Display.SetScaleArray = ['POINTS', 'Result']
				contour3Display.ScaleTransferFunction = 'PiecewiseFunction'
				contour3Display.OpacityArray = ['POINTS', 'Result']
				contour3Display.OpacityTransferFunction = 'PiecewiseFunction'
				
				# save data
				SaveData(pathVort, proxy=contour2)
				print(pathVort)
				
				## Properties modified on contour1
				#contour2.Isosurfaces = [150.0]
				
				# save data
				#SaveData(pathVortB, proxy=contour2)
				#print(pathVortB)
				
				# Properties modified on contour1
				contour2.Isosurfaces = [300.0]
				
				# save data
				SaveData(pathVortC, proxy=contour2)
				print(pathVortC)
