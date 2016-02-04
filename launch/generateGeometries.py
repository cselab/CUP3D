import os
import re

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

rootDir = '/scratch/daint/cconti/FallingSamaraFixed_24h_1024_310116/'

for dirName, subDirList, fileList in os.walk(rootDir):
	print(dirName)
	for f in fileList:
		if f.endswith('.xmf'):
			frame = re.search('FallingSamaraFixed_24h_1024_310116-(.+?).xmf', f)
			pathChi = dirName+'/Chi'+frame.group(1)+'.ply'
			pathVort = dirName+'/Vorticity'+frame.group(1)+'.ply'
			
			if not os.path.isfile(pathChi):
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

				# create a new 'Calculator'
				calculator1 = Calculator(Input=filename)
				calculator1.Function = ''

				# Properties modified on calculator1
				calculator1.Function = 'data_4'

				# show data in view
				calculator1Display = Show(calculator1, renderView1)
				# trace defaults for the display properties.
				calculator1Display.Representation = 'Outline'
				calculator1Display.ColorArrayName = ['POINTS', '']
				calculator1Display.GlyphType = 'Arrow'
				calculator1Display.ScalarOpacityUnitDistance = 0.0016914558667664823
				calculator1Display.Slice = 511

				# hide data in view
				Hide(filename, renderView1)

				# create a new 'Contour'
				contour1 = Contour(Input=calculator1)
				contour1.ContourBy = ['POINTS', 'Result']
				contour1.Isosurfaces = [0.5]
				contour1.PointMergeMethod = 'Uniform Binning'

				# show data in view
				contour1Display = Show(contour1, renderView1)
				# trace defaults for the display properties.
				contour1Display.ColorArrayName = [None, '']
				contour1Display.GlyphType = 'Arrow'
				contour1Display.SetScaleArray = [None, '']
				contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
				contour1Display.OpacityArray = [None, '']
				contour1Display.OpacityTransferFunction = 'PiecewiseFunction'

				# set active source
				SetActiveSource(filename)

				# set active source
				SetActiveSource(contour1)

				# save data
				SaveData(pathChi, proxy=contour1)
				print(pathChi)

				# set active source
				SetActiveSource(filename)

				# create a new 'Calculator'
				calculator2 = Calculator(Input=filename)
				calculator2.Function = ''

				# Properties modified on calculator2
				calculator2.Function = 'data_6'

				# show data in view
				calculator2Display = Show(calculator2, renderView1)
				# trace defaults for the display properties.
				calculator2Display.Representation = 'Outline'
				calculator2Display.ColorArrayName = ['POINTS', '']
				calculator2Display.GlyphType = 'Arrow'
				calculator2Display.ScalarOpacityUnitDistance = 0.0016772642744068997
				calculator2Display.Slice = 511

				# hide data in view
				Hide(filename, renderView1)

				# create a new 'Contour'
				contour2 = Contour(Input=calculator2)
				contour2.ContourBy = ['POINTS', 'Result']
				contour2.Isosurfaces = [977.6132202198063]
				contour2.PointMergeMethod = 'Uniform Binning'

				# Properties modified on contour2
				contour2.Isosurfaces = [50.0]

				# show data in view
				contour2Display = Show(contour2, renderView1)
				# trace defaults for the display properties.
				contour2Display.ColorArrayName = [None, '']
				contour2Display.GlyphType = 'Arrow'
				contour2Display.SetScaleArray = [None, '']
				contour2Display.ScaleTransferFunction = 'PiecewiseFunction'
				contour2Display.OpacityArray = [None, '']
				contour2Display.OpacityTransferFunction = 'PiecewiseFunction'

				# save data
				SaveData(pathVort, proxy=contour2)
				print(pathVort)
