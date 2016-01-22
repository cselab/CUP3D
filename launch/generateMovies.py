Version = "4.4"

import os

#### import the simple module from the paraview
#from paraview.simple import *
try: paraview.simple
except: from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

rootDir = '/scratch/daint/cconti/'

for dirName, subDirList, fileList in os.walk(rootDir):
	print(dirName)
	if "FallingSamaraFixed_MPI32" in dirName and "150116" in dirName and any(fname.endswith('.xmf') for fname in fileList):

		# create a new 'XDMF Reader'
		files=[dirName+'/'+f for f in fileList if f.endswith('.xmf')]
		print(files)
		filename = XDMFReader(FileNames=[dirName+'/'+f for f in fileList if f.endswith('.xmf')])

		filename.PointArrayStatus = ['data']
		filename.CellArrayStatus = []
		filename.SetStatus = []
		filename.GridStatus = []
		filename.Stride = [1, 1, 1]

		# get animation scene
		animationScene1 = GetAnimationScene()

		# update animation scene based on data timesteps
		animationScene1.UpdateAnimationUsingDataTimeSteps()

		# Properties modified on filename
		filename.GridStatus = ['Grid_479']

		# get active view
		renderView1 = GetActiveViewOrCreate('RenderView')
		# uncomment following to set a specific view size
		# renderView1.ViewSize = [1945, 1099]

		# show data in view
		filenameDisplay = Show(filename, renderView1)
		# trace defaults for the display properties.
		filenameDisplay.Representation = 'Outline'
		filenameDisplay.AmbientColor = [1.0, 1.0, 1.0]
		filenameDisplay.ColorArrayName = [None, '']
		filenameDisplay.DiffuseColor = [1.0, 1.0, 1.0]
		filenameDisplay.LookupTable = None
		filenameDisplay.MapScalars = 1
		filenameDisplay.InterpolateScalarsBeforeMapping = 1
		filenameDisplay.Opacity = 1.0
		filenameDisplay.PointSize = 2.0
		filenameDisplay.LineWidth = 1.0
		filenameDisplay.Interpolation = 'Gouraud'
		filenameDisplay.Specular = 0.0
		filenameDisplay.SpecularColor = [1.0, 1.0, 1.0]
		filenameDisplay.SpecularPower = 100.0
		filenameDisplay.Ambient = 0.0
		filenameDisplay.Diffuse = 1.0
		filenameDisplay.EdgeColor = [0.0, 0.0, 0.5]
		filenameDisplay.BackfaceRepresentation = 'Follow Frontface'
		filenameDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
		filenameDisplay.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
		filenameDisplay.BackfaceOpacity = 1.0
		filenameDisplay.Position = [0.0, 0.0, 0.0]
		filenameDisplay.Scale = [1.0, 1.0, 1.0]
		filenameDisplay.Orientation = [0.0, 0.0, 0.0]
		filenameDisplay.Origin = [0.0, 0.0, 0.0]
		filenameDisplay.Pickable = 1
		filenameDisplay.Texture = None
		filenameDisplay.Triangulate = 0
		filenameDisplay.NonlinearSubdivisionLevel = 1
		filenameDisplay.CubeAxesColor = [1.0, 1.0, 1.0]
		filenameDisplay.CubeAxesCornerOffset = 0.0
		filenameDisplay.CubeAxesFlyMode = 'Closest Triad'
		filenameDisplay.CubeAxesInertia = 1
		filenameDisplay.CubeAxesTickLocation = 'Inside'
		filenameDisplay.CubeAxesXAxisMinorTickVisibility = 1
		filenameDisplay.CubeAxesXAxisTickVisibility = 1
		filenameDisplay.CubeAxesXAxisVisibility = 1
		filenameDisplay.CubeAxesXGridLines = 0
		filenameDisplay.CubeAxesXTitle = 'X-Axis'
		filenameDisplay.CubeAxesUseDefaultXTitle = 1
		filenameDisplay.CubeAxesYAxisMinorTickVisibility = 1
		filenameDisplay.CubeAxesYAxisTickVisibility = 1
		filenameDisplay.CubeAxesYAxisVisibility = 1
		filenameDisplay.CubeAxesYGridLines = 0
		filenameDisplay.CubeAxesYTitle = 'Y-Axis'
		filenameDisplay.CubeAxesUseDefaultYTitle = 1
		filenameDisplay.CubeAxesZAxisMinorTickVisibility = 1
		filenameDisplay.CubeAxesZAxisTickVisibility = 1
		filenameDisplay.CubeAxesZAxisVisibility = 1
		filenameDisplay.CubeAxesZGridLines = 0
		filenameDisplay.CubeAxesZTitle = 'Z-Axis'
		filenameDisplay.CubeAxesUseDefaultZTitle = 1
		filenameDisplay.CubeAxesGridLineLocation = 'All Faces'
		filenameDisplay.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		filenameDisplay.CustomBoundsActive = [0, 0, 0]
		filenameDisplay.OriginalBoundsRangeActive = [0, 0, 0]
		filenameDisplay.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		filenameDisplay.CustomRangeActive = [0, 0, 0]
		filenameDisplay.UseAxesOrigin = 0
		filenameDisplay.AxesOrigin = [0.0, 0.0, 0.0]
		filenameDisplay.CubeAxesXLabelFormat = '%-#6.3g'
		filenameDisplay.CubeAxesYLabelFormat = '%-#6.3g'
		filenameDisplay.CubeAxesZLabelFormat = '%-#6.3g'
		filenameDisplay.StickyAxes = 0
		filenameDisplay.CenterStickyAxes = 0
		filenameDisplay.SelectionCellLabelBold = 0
		filenameDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
		filenameDisplay.SelectionCellLabelFontFamily = 'Arial'
		filenameDisplay.SelectionCellLabelFontSize = 18
		filenameDisplay.SelectionCellLabelItalic = 0
		filenameDisplay.SelectionCellLabelJustification = 'Left'
		filenameDisplay.SelectionCellLabelOpacity = 1.0
		filenameDisplay.SelectionCellLabelShadow = 0
		filenameDisplay.SelectionPointLabelBold = 0
		filenameDisplay.SelectionPointLabelColor = [0.5, 0.5, 0.5]
		filenameDisplay.SelectionPointLabelFontFamily = 'Arial'
		filenameDisplay.SelectionPointLabelFontSize = 18
		filenameDisplay.SelectionPointLabelItalic = 0
		filenameDisplay.SelectionPointLabelJustification = 'Left'
		filenameDisplay.SelectionPointLabelOpacity = 1.0
		filenameDisplay.SelectionPointLabelShadow = 0
		filenameDisplay.ScalarOpacityUnitDistance = 0.006765823467065929
		filenameDisplay.VolumeRenderingMode = 'Smart'
		filenameDisplay.Shade = 0
		filenameDisplay.SliceMode = 'XY Plane'
		filenameDisplay.Slice = 127

		# reset view to fit data
		renderView1.ResetCamera()

		# create a new 'Calculator'
		calculator1 = Calculator(Input=filename)
		calculator1.AttributeMode = 'Point Data'
		calculator1.CoordinateResults = 0
		calculator1.ResultNormals = 0
		calculator1.ResultTCoords = 0
		calculator1.ResultArrayName = 'Result'
		calculator1.Function = ''
		calculator1.ReplaceInvalidResults = 1
		calculator1.ReplacementValue = 0.0

		# Properties modified on calculator1
		calculator1.Function = 'data_4'

		# show data in view
		calculator1Display = Show(calculator1, renderView1)
		# trace defaults for the display properties.
		calculator1Display.Representation = 'Outline'
		calculator1Display.AmbientColor = [1.0, 1.0, 1.0]
		calculator1Display.ColorArrayName = ['POINTS', '']
		calculator1Display.DiffuseColor = [1.0, 1.0, 1.0]
		calculator1Display.LookupTable = None
		calculator1Display.MapScalars = 1
		calculator1Display.InterpolateScalarsBeforeMapping = 1
		calculator1Display.Opacity = 1.0
		calculator1Display.PointSize = 2.0
		calculator1Display.LineWidth = 1.0
		calculator1Display.Interpolation = 'Gouraud'
		calculator1Display.Specular = 0.0
		calculator1Display.SpecularColor = [1.0, 1.0, 1.0]
		calculator1Display.SpecularPower = 100.0
		calculator1Display.Ambient = 0.0
		calculator1Display.Diffuse = 1.0
		calculator1Display.EdgeColor = [0.0, 0.0, 0.5]
		calculator1Display.BackfaceRepresentation = 'Follow Frontface'
		calculator1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
		calculator1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
		calculator1Display.BackfaceOpacity = 1.0
		calculator1Display.Position = [0.0, 0.0, 0.0]
		calculator1Display.Scale = [1.0, 1.0, 1.0]
		calculator1Display.Orientation = [0.0, 0.0, 0.0]
		calculator1Display.Origin = [0.0, 0.0, 0.0]
		calculator1Display.Pickable = 1
		calculator1Display.Texture = None
		calculator1Display.Triangulate = 0
		calculator1Display.NonlinearSubdivisionLevel = 1
		calculator1Display.CubeAxesColor = [1.0, 1.0, 1.0]
		calculator1Display.CubeAxesCornerOffset = 0.0
		calculator1Display.CubeAxesFlyMode = 'Closest Triad'
		calculator1Display.CubeAxesInertia = 1
		calculator1Display.CubeAxesTickLocation = 'Inside'
		calculator1Display.CubeAxesXAxisMinorTickVisibility = 1
		calculator1Display.CubeAxesXAxisTickVisibility = 1
		calculator1Display.CubeAxesXAxisVisibility = 1
		calculator1Display.CubeAxesXGridLines = 0
		calculator1Display.CubeAxesXTitle = 'X-Axis'
		calculator1Display.CubeAxesUseDefaultXTitle = 1
		calculator1Display.CubeAxesYAxisMinorTickVisibility = 1
		calculator1Display.CubeAxesYAxisTickVisibility = 1
		calculator1Display.CubeAxesYAxisVisibility = 1
		calculator1Display.CubeAxesYGridLines = 0
		calculator1Display.CubeAxesYTitle = 'Y-Axis'
		calculator1Display.CubeAxesUseDefaultYTitle = 1
		calculator1Display.CubeAxesZAxisMinorTickVisibility = 1
		calculator1Display.CubeAxesZAxisTickVisibility = 1
		calculator1Display.CubeAxesZAxisVisibility = 1
		calculator1Display.CubeAxesZGridLines = 0
		calculator1Display.CubeAxesZTitle = 'Z-Axis'
		calculator1Display.CubeAxesUseDefaultZTitle = 1
		calculator1Display.CubeAxesGridLineLocation = 'All Faces'
		calculator1Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		calculator1Display.CustomBoundsActive = [0, 0, 0]
		calculator1Display.OriginalBoundsRangeActive = [0, 0, 0]
		calculator1Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		calculator1Display.CustomRangeActive = [0, 0, 0]
		calculator1Display.UseAxesOrigin = 0
		calculator1Display.AxesOrigin = [0.0, 0.0, 0.0]
		calculator1Display.CubeAxesXLabelFormat = '%-#6.3g'
		calculator1Display.CubeAxesYLabelFormat = '%-#6.3g'
		calculator1Display.CubeAxesZLabelFormat = '%-#6.3g'
		calculator1Display.StickyAxes = 0
		calculator1Display.CenterStickyAxes = 0
		calculator1Display.SelectionCellLabelBold = 0
		calculator1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
		calculator1Display.SelectionCellLabelFontFamily = 'Arial'
		calculator1Display.SelectionCellLabelFontSize = 18
		calculator1Display.SelectionCellLabelItalic = 0
		calculator1Display.SelectionCellLabelJustification = 'Left'
		calculator1Display.SelectionCellLabelOpacity = 1.0
		calculator1Display.SelectionCellLabelShadow = 0
		calculator1Display.SelectionPointLabelBold = 0
		calculator1Display.SelectionPointLabelColor = [0.5, 0.5, 0.5]
		calculator1Display.SelectionPointLabelFontFamily = 'Arial'
		calculator1Display.SelectionPointLabelFontSize = 18
		calculator1Display.SelectionPointLabelItalic = 0
		calculator1Display.SelectionPointLabelJustification = 'Left'
		calculator1Display.SelectionPointLabelOpacity = 1.0
		calculator1Display.SelectionPointLabelShadow = 0
		calculator1Display.ScalarOpacityUnitDistance = 0.006765823467065929
		calculator1Display.VolumeRenderingMode = 'Smart'
		calculator1Display.Shade = 0
		calculator1Display.SliceMode = 'XY Plane'
		calculator1Display.Slice = 127

		# hide data in view
		Hide(filename, renderView1)

		# create a new 'Contour'
		contour1 = Contour(Input=calculator1)
		contour1.ContourBy = ['POINTS', 'Result']
		contour1.ComputeNormals = 1
		contour1.ComputeGradients = 0
		contour1.ComputeScalars = 0
		contour1.GenerateTriangles = 1
		contour1.Isosurfaces = [0.5]
		contour1.PointMergeMethod = 'Uniform Binning'

		# init the 'Uniform Binning' selected for 'PointMergeMethod'
		contour1.PointMergeMethod.Divisions = [50, 50, 50]
		contour1.PointMergeMethod.Numberofpointsperbucket = 8

		# show data in view
		contour1Display = Show(contour1, renderView1)
		# trace defaults for the display properties.
		contour1Display.Representation = 'Surface'
		contour1Display.AmbientColor = [1.0, 1.0, 1.0]
		contour1Display.ColorArrayName = [None, '']
		contour1Display.DiffuseColor = [1.0, 1.0, 1.0]
		contour1Display.LookupTable = None
		contour1Display.MapScalars = 1
		contour1Display.InterpolateScalarsBeforeMapping = 1
		contour1Display.Opacity = 1.0
		contour1Display.PointSize = 2.0
		contour1Display.LineWidth = 1.0
		contour1Display.Interpolation = 'Gouraud'
		contour1Display.Specular = 0.0
		contour1Display.SpecularColor = [1.0, 1.0, 1.0]
		contour1Display.SpecularPower = 100.0
		contour1Display.Ambient = 0.0
		contour1Display.Diffuse = 1.0
		contour1Display.EdgeColor = [0.0, 0.0, 0.5]
		contour1Display.BackfaceRepresentation = 'Follow Frontface'
		contour1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
		contour1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
		contour1Display.BackfaceOpacity = 1.0
		contour1Display.Position = [0.0, 0.0, 0.0]
		contour1Display.Scale = [1.0, 1.0, 1.0]
		contour1Display.Orientation = [0.0, 0.0, 0.0]
		contour1Display.Origin = [0.0, 0.0, 0.0]
		contour1Display.Pickable = 1
		contour1Display.Texture = None
		contour1Display.Triangulate = 0
		contour1Display.NonlinearSubdivisionLevel = 1
		contour1Display.CubeAxesColor = [1.0, 1.0, 1.0]
		contour1Display.CubeAxesCornerOffset = 0.0
		contour1Display.CubeAxesFlyMode = 'Closest Triad'
		contour1Display.CubeAxesInertia = 1
		contour1Display.CubeAxesTickLocation = 'Inside'
		contour1Display.CubeAxesXAxisMinorTickVisibility = 1
		contour1Display.CubeAxesXAxisTickVisibility = 1
		contour1Display.CubeAxesXAxisVisibility = 1
		contour1Display.CubeAxesXGridLines = 0
		contour1Display.CubeAxesXTitle = 'X-Axis'
		contour1Display.CubeAxesUseDefaultXTitle = 1
		contour1Display.CubeAxesYAxisMinorTickVisibility = 1
		contour1Display.CubeAxesYAxisTickVisibility = 1
		contour1Display.CubeAxesYAxisVisibility = 1
		contour1Display.CubeAxesYGridLines = 0
		contour1Display.CubeAxesYTitle = 'Y-Axis'
		contour1Display.CubeAxesUseDefaultYTitle = 1
		contour1Display.CubeAxesZAxisMinorTickVisibility = 1
		contour1Display.CubeAxesZAxisTickVisibility = 1
		contour1Display.CubeAxesZAxisVisibility = 1
		contour1Display.CubeAxesZGridLines = 0
		contour1Display.CubeAxesZTitle = 'Z-Axis'
		contour1Display.CubeAxesUseDefaultZTitle = 1
		contour1Display.CubeAxesGridLineLocation = 'All Faces'
		contour1Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		contour1Display.CustomBoundsActive = [0, 0, 0]
		contour1Display.OriginalBoundsRangeActive = [0, 0, 0]
		contour1Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		contour1Display.CustomRangeActive = [0, 0, 0]
		contour1Display.UseAxesOrigin = 0
		contour1Display.AxesOrigin = [0.0, 0.0, 0.0]
		contour1Display.CubeAxesXLabelFormat = '%-#6.3g'
		contour1Display.CubeAxesYLabelFormat = '%-#6.3g'
		contour1Display.CubeAxesZLabelFormat = '%-#6.3g'
		contour1Display.StickyAxes = 0
		contour1Display.CenterStickyAxes = 0
		contour1Display.SelectionCellLabelBold = 0
		contour1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
		contour1Display.SelectionCellLabelFontFamily = 'Arial'
		contour1Display.SelectionCellLabelFontSize = 18
		contour1Display.SelectionCellLabelItalic = 0
		contour1Display.SelectionCellLabelJustification = 'Left'
		contour1Display.SelectionCellLabelOpacity = 1.0
		contour1Display.SelectionCellLabelShadow = 0
		contour1Display.SelectionPointLabelBold = 0
		contour1Display.SelectionPointLabelColor = [0.5, 0.5, 0.5]
		contour1Display.SelectionPointLabelFontFamily = 'Arial'
		contour1Display.SelectionPointLabelFontSize = 18
		contour1Display.SelectionPointLabelItalic = 0
		contour1Display.SelectionPointLabelJustification = 'Left'
		contour1Display.SelectionPointLabelOpacity = 1.0
		contour1Display.SelectionPointLabelShadow = 0

		# hide data in view
		Hide(calculator1, renderView1)

		# reset view to fit data
		renderView1.ResetCamera()

		# set active source
		SetActiveSource(calculator1)

		# show data in view
		calculator1Display = Show(calculator1, renderView1)

		# set active source
		SetActiveSource(filename)

		# create a new 'Calculator'
		calculator2 = Calculator(Input=filename)
		calculator2.AttributeMode = 'Point Data'
		calculator2.CoordinateResults = 0
		calculator2.ResultNormals = 0
		calculator2.ResultTCoords = 0
		calculator2.ResultArrayName = 'Result'
		calculator2.Function = ''
		calculator2.ReplaceInvalidResults = 1
		calculator2.ReplacementValue = 0.0

		# Properties modified on calculator2
		calculator2.Function = 'data_6'

		# show data in view
		calculator2Display = Show(calculator2, renderView1)
		# trace defaults for the display properties.
		calculator2Display.Representation = 'Outline'
		calculator2Display.AmbientColor = [1.0, 1.0, 1.0]
		calculator2Display.ColorArrayName = ['POINTS', '']
		calculator2Display.DiffuseColor = [1.0, 1.0, 1.0]
		calculator2Display.LookupTable = None
		calculator2Display.MapScalars = 1
		calculator2Display.InterpolateScalarsBeforeMapping = 1
		calculator2Display.Opacity = 1.0
		calculator2Display.PointSize = 2.0
		calculator2Display.LineWidth = 1.0
		calculator2Display.Interpolation = 'Gouraud'
		calculator2Display.Specular = 0.0
		calculator2Display.SpecularColor = [1.0, 1.0, 1.0]
		calculator2Display.SpecularPower = 100.0
		calculator2Display.Ambient = 0.0
		calculator2Display.Diffuse = 1.0
		calculator2Display.EdgeColor = [0.0, 0.0, 0.5]
		calculator2Display.BackfaceRepresentation = 'Follow Frontface'
		calculator2Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
		calculator2Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
		calculator2Display.BackfaceOpacity = 1.0
		calculator2Display.Position = [0.0, 0.0, 0.0]
		calculator2Display.Scale = [1.0, 1.0, 1.0]
		calculator2Display.Orientation = [0.0, 0.0, 0.0]
		calculator2Display.Origin = [0.0, 0.0, 0.0]
		calculator2Display.Pickable = 1
		calculator2Display.Texture = None
		calculator2Display.Triangulate = 0
		calculator2Display.NonlinearSubdivisionLevel = 1
		calculator2Display.CubeAxesColor = [1.0, 1.0, 1.0]
		calculator2Display.CubeAxesCornerOffset = 0.0
		calculator2Display.CubeAxesFlyMode = 'Closest Triad'
		calculator2Display.CubeAxesInertia = 1
		calculator2Display.CubeAxesTickLocation = 'Inside'
		calculator2Display.CubeAxesXAxisMinorTickVisibility = 1
		calculator2Display.CubeAxesXAxisTickVisibility = 1
		calculator2Display.CubeAxesXAxisVisibility = 1
		calculator2Display.CubeAxesXGridLines = 0
		calculator2Display.CubeAxesXTitle = 'X-Axis'
		calculator2Display.CubeAxesUseDefaultXTitle = 1
		calculator2Display.CubeAxesYAxisMinorTickVisibility = 1
		calculator2Display.CubeAxesYAxisTickVisibility = 1
		calculator2Display.CubeAxesYAxisVisibility = 1
		calculator2Display.CubeAxesYGridLines = 0
		calculator2Display.CubeAxesYTitle = 'Y-Axis'
		calculator2Display.CubeAxesUseDefaultYTitle = 1
		calculator2Display.CubeAxesZAxisMinorTickVisibility = 1
		calculator2Display.CubeAxesZAxisTickVisibility = 1
		calculator2Display.CubeAxesZAxisVisibility = 1
		calculator2Display.CubeAxesZGridLines = 0
		calculator2Display.CubeAxesZTitle = 'Z-Axis'
		calculator2Display.CubeAxesUseDefaultZTitle = 1
		calculator2Display.CubeAxesGridLineLocation = 'All Faces'
		calculator2Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		calculator2Display.CustomBoundsActive = [0, 0, 0]
		calculator2Display.OriginalBoundsRangeActive = [0, 0, 0]
		calculator2Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		calculator2Display.CustomRangeActive = [0, 0, 0]
		calculator2Display.UseAxesOrigin = 0
		calculator2Display.AxesOrigin = [0.0, 0.0, 0.0]
		calculator2Display.CubeAxesXLabelFormat = '%-#6.3g'
		calculator2Display.CubeAxesYLabelFormat = '%-#6.3g'
		calculator2Display.CubeAxesZLabelFormat = '%-#6.3g'
		calculator2Display.StickyAxes = 0
		calculator2Display.CenterStickyAxes = 0
		calculator2Display.SelectionCellLabelBold = 0
		calculator2Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
		calculator2Display.SelectionCellLabelFontFamily = 'Arial'
		calculator2Display.SelectionCellLabelFontSize = 18
		calculator2Display.SelectionCellLabelItalic = 0
		calculator2Display.SelectionCellLabelJustification = 'Left'
		calculator2Display.SelectionCellLabelOpacity = 1.0
		calculator2Display.SelectionCellLabelShadow = 0
		calculator2Display.SelectionPointLabelBold = 0
		calculator2Display.SelectionPointLabelColor = [0.5, 0.5, 0.5]
		calculator2Display.SelectionPointLabelFontFamily = 'Arial'
		calculator2Display.SelectionPointLabelFontSize = 18
		calculator2Display.SelectionPointLabelItalic = 0
		calculator2Display.SelectionPointLabelJustification = 'Left'
		calculator2Display.SelectionPointLabelOpacity = 1.0
		calculator2Display.SelectionPointLabelShadow = 0
		calculator2Display.ScalarOpacityUnitDistance = 0.006765823467065929
		calculator2Display.VolumeRenderingMode = 'Smart'
		calculator2Display.Shade = 0
		calculator2Display.SliceMode = 'XY Plane'
		calculator2Display.Slice = 127

		# hide data in view
		Hide(filename, renderView1)

		# set active source
		SetActiveSource(filename)

		# hide data in view
		Hide(calculator2, renderView1)

		# show data in view
		filenameDisplay = Show(filename, renderView1)

		# destroy calculator2
		Delete(calculator2)
		del calculator2

		# create a new 'Calculator'
		calculator2 = Calculator(Input=filename)
		calculator2.AttributeMode = 'Point Data'
		calculator2.CoordinateResults = 0
		calculator2.ResultNormals = 0
		calculator2.ResultTCoords = 0
		calculator2.ResultArrayName = 'Result'
		calculator2.Function = ''
		calculator2.ReplaceInvalidResults = 1
		calculator2.ReplacementValue = 0.0

		# Properties modified on calculator2
		calculator2.Function = 'data_5'

		# show data in view
		calculator2Display = Show(calculator2, renderView1)
		# trace defaults for the display properties.
		calculator2Display.Representation = 'Outline'
		calculator2Display.AmbientColor = [1.0, 1.0, 1.0]
		calculator2Display.ColorArrayName = ['POINTS', '']
		calculator2Display.DiffuseColor = [1.0, 1.0, 1.0]
		calculator2Display.LookupTable = None
		calculator2Display.MapScalars = 1
		calculator2Display.InterpolateScalarsBeforeMapping = 1
		calculator2Display.Opacity = 1.0
		calculator2Display.PointSize = 2.0
		calculator2Display.LineWidth = 1.0
		calculator2Display.Interpolation = 'Gouraud'
		calculator2Display.Specular = 0.0
		calculator2Display.SpecularColor = [1.0, 1.0, 1.0]
		calculator2Display.SpecularPower = 100.0
		calculator2Display.Ambient = 0.0
		calculator2Display.Diffuse = 1.0
		calculator2Display.EdgeColor = [0.0, 0.0, 0.5]
		calculator2Display.BackfaceRepresentation = 'Follow Frontface'
		calculator2Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
		calculator2Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
		calculator2Display.BackfaceOpacity = 1.0
		calculator2Display.Position = [0.0, 0.0, 0.0]
		calculator2Display.Scale = [1.0, 1.0, 1.0]
		calculator2Display.Orientation = [0.0, 0.0, 0.0]
		calculator2Display.Origin = [0.0, 0.0, 0.0]
		calculator2Display.Pickable = 1
		calculator2Display.Texture = None
		calculator2Display.Triangulate = 0
		calculator2Display.NonlinearSubdivisionLevel = 1
		calculator2Display.CubeAxesColor = [1.0, 1.0, 1.0]
		calculator2Display.CubeAxesCornerOffset = 0.0
		calculator2Display.CubeAxesFlyMode = 'Closest Triad'
		calculator2Display.CubeAxesInertia = 1
		calculator2Display.CubeAxesTickLocation = 'Inside'
		calculator2Display.CubeAxesXAxisMinorTickVisibility = 1
		calculator2Display.CubeAxesXAxisTickVisibility = 1
		calculator2Display.CubeAxesXAxisVisibility = 1
		calculator2Display.CubeAxesXGridLines = 0
		calculator2Display.CubeAxesXTitle = 'X-Axis'
		calculator2Display.CubeAxesUseDefaultXTitle = 1
		calculator2Display.CubeAxesYAxisMinorTickVisibility = 1
		calculator2Display.CubeAxesYAxisTickVisibility = 1
		calculator2Display.CubeAxesYAxisVisibility = 1
		calculator2Display.CubeAxesYGridLines = 0
		calculator2Display.CubeAxesYTitle = 'Y-Axis'
		calculator2Display.CubeAxesUseDefaultYTitle = 1
		calculator2Display.CubeAxesZAxisMinorTickVisibility = 1
		calculator2Display.CubeAxesZAxisTickVisibility = 1
		calculator2Display.CubeAxesZAxisVisibility = 1
		calculator2Display.CubeAxesZGridLines = 0
		calculator2Display.CubeAxesZTitle = 'Z-Axis'
		calculator2Display.CubeAxesUseDefaultZTitle = 1
		calculator2Display.CubeAxesGridLineLocation = 'All Faces'
		calculator2Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		calculator2Display.CustomBoundsActive = [0, 0, 0]
		calculator2Display.OriginalBoundsRangeActive = [0, 0, 0]
		calculator2Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		calculator2Display.CustomRangeActive = [0, 0, 0]
		calculator2Display.UseAxesOrigin = 0
		calculator2Display.AxesOrigin = [0.0, 0.0, 0.0]
		calculator2Display.CubeAxesXLabelFormat = '%-#6.3g'
		calculator2Display.CubeAxesYLabelFormat = '%-#6.3g'
		calculator2Display.CubeAxesZLabelFormat = '%-#6.3g'
		calculator2Display.StickyAxes = 0
		calculator2Display.CenterStickyAxes = 0
		calculator2Display.SelectionCellLabelBold = 0
		calculator2Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
		calculator2Display.SelectionCellLabelFontFamily = 'Arial'
		calculator2Display.SelectionCellLabelFontSize = 18
		calculator2Display.SelectionCellLabelItalic = 0
		calculator2Display.SelectionCellLabelJustification = 'Left'
		calculator2Display.SelectionCellLabelOpacity = 1.0
		calculator2Display.SelectionCellLabelShadow = 0
		calculator2Display.SelectionPointLabelBold = 0
		calculator2Display.SelectionPointLabelColor = [0.5, 0.5, 0.5]
		calculator2Display.SelectionPointLabelFontFamily = 'Arial'
		calculator2Display.SelectionPointLabelFontSize = 18
		calculator2Display.SelectionPointLabelItalic = 0
		calculator2Display.SelectionPointLabelJustification = 'Left'
		calculator2Display.SelectionPointLabelOpacity = 1.0
		calculator2Display.SelectionPointLabelShadow = 0
		calculator2Display.ScalarOpacityUnitDistance = 0.006765823467065929
		calculator2Display.VolumeRenderingMode = 'Smart'
		calculator2Display.Shade = 0
		calculator2Display.SliceMode = 'XY Plane'
		calculator2Display.Slice = 127

		# hide data in view
		Hide(filename, renderView1)

		# set active source
		SetActiveSource(filename)

		# hide data in view
		Hide(calculator2, renderView1)

		# show data in view
		filenameDisplay = Show(filename, renderView1)

		# destroy calculator2
		Delete(calculator2)
		del calculator2

		# create a new 'Calculator'
		calculator2 = Calculator(Input=filename)
		calculator2.AttributeMode = 'Point Data'
		calculator2.CoordinateResults = 0
		calculator2.ResultNormals = 0
		calculator2.ResultTCoords = 0
		calculator2.ResultArrayName = 'Result'
		calculator2.Function = ''
		calculator2.ReplaceInvalidResults = 1
		calculator2.ReplacementValue = 0.0

		# Properties modified on calculator2
		calculator2.Function = 'data_6'

		# show data in view
		calculator2Display = Show(calculator2, renderView1)
		# trace defaults for the display properties.
		calculator2Display.Representation = 'Outline'
		calculator2Display.AmbientColor = [1.0, 1.0, 1.0]
		calculator2Display.ColorArrayName = ['POINTS', '']
		calculator2Display.DiffuseColor = [1.0, 1.0, 1.0]
		calculator2Display.LookupTable = None
		calculator2Display.MapScalars = 1
		calculator2Display.InterpolateScalarsBeforeMapping = 1
		calculator2Display.Opacity = 1.0
		calculator2Display.PointSize = 2.0
		calculator2Display.LineWidth = 1.0
		calculator2Display.Interpolation = 'Gouraud'
		calculator2Display.Specular = 0.0
		calculator2Display.SpecularColor = [1.0, 1.0, 1.0]
		calculator2Display.SpecularPower = 100.0
		calculator2Display.Ambient = 0.0
		calculator2Display.Diffuse = 1.0
		calculator2Display.EdgeColor = [0.0, 0.0, 0.5]
		calculator2Display.BackfaceRepresentation = 'Follow Frontface'
		calculator2Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
		calculator2Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
		calculator2Display.BackfaceOpacity = 1.0
		calculator2Display.Position = [0.0, 0.0, 0.0]
		calculator2Display.Scale = [1.0, 1.0, 1.0]
		calculator2Display.Orientation = [0.0, 0.0, 0.0]
		calculator2Display.Origin = [0.0, 0.0, 0.0]
		calculator2Display.Pickable = 1
		calculator2Display.Texture = None
		calculator2Display.Triangulate = 0
		calculator2Display.NonlinearSubdivisionLevel = 1
		calculator2Display.CubeAxesColor = [1.0, 1.0, 1.0]
		calculator2Display.CubeAxesCornerOffset = 0.0
		calculator2Display.CubeAxesFlyMode = 'Closest Triad'
		calculator2Display.CubeAxesInertia = 1
		calculator2Display.CubeAxesTickLocation = 'Inside'
		calculator2Display.CubeAxesXAxisMinorTickVisibility = 1
		calculator2Display.CubeAxesXAxisTickVisibility = 1
		calculator2Display.CubeAxesXAxisVisibility = 1
		calculator2Display.CubeAxesXGridLines = 0
		calculator2Display.CubeAxesXTitle = 'X-Axis'
		calculator2Display.CubeAxesUseDefaultXTitle = 1
		calculator2Display.CubeAxesYAxisMinorTickVisibility = 1
		calculator2Display.CubeAxesYAxisTickVisibility = 1
		calculator2Display.CubeAxesYAxisVisibility = 1
		calculator2Display.CubeAxesYGridLines = 0
		calculator2Display.CubeAxesYTitle = 'Y-Axis'
		calculator2Display.CubeAxesUseDefaultYTitle = 1
		calculator2Display.CubeAxesZAxisMinorTickVisibility = 1
		calculator2Display.CubeAxesZAxisTickVisibility = 1
		calculator2Display.CubeAxesZAxisVisibility = 1
		calculator2Display.CubeAxesZGridLines = 0
		calculator2Display.CubeAxesZTitle = 'Z-Axis'
		calculator2Display.CubeAxesUseDefaultZTitle = 1
		calculator2Display.CubeAxesGridLineLocation = 'All Faces'
		calculator2Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		calculator2Display.CustomBoundsActive = [0, 0, 0]
		calculator2Display.OriginalBoundsRangeActive = [0, 0, 0]
		calculator2Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		calculator2Display.CustomRangeActive = [0, 0, 0]
		calculator2Display.UseAxesOrigin = 0
		calculator2Display.AxesOrigin = [0.0, 0.0, 0.0]
		calculator2Display.CubeAxesXLabelFormat = '%-#6.3g'
		calculator2Display.CubeAxesYLabelFormat = '%-#6.3g'
		calculator2Display.CubeAxesZLabelFormat = '%-#6.3g'
		calculator2Display.StickyAxes = 0
		calculator2Display.CenterStickyAxes = 0
		calculator2Display.SelectionCellLabelBold = 0
		calculator2Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
		calculator2Display.SelectionCellLabelFontFamily = 'Arial'
		calculator2Display.SelectionCellLabelFontSize = 18
		calculator2Display.SelectionCellLabelItalic = 0
		calculator2Display.SelectionCellLabelJustification = 'Left'
		calculator2Display.SelectionCellLabelOpacity = 1.0
		calculator2Display.SelectionCellLabelShadow = 0
		calculator2Display.SelectionPointLabelBold = 0
		calculator2Display.SelectionPointLabelColor = [0.5, 0.5, 0.5]
		calculator2Display.SelectionPointLabelFontFamily = 'Arial'
		calculator2Display.SelectionPointLabelFontSize = 18
		calculator2Display.SelectionPointLabelItalic = 0
		calculator2Display.SelectionPointLabelJustification = 'Left'
		calculator2Display.SelectionPointLabelOpacity = 1.0
		calculator2Display.SelectionPointLabelShadow = 0
		calculator2Display.ScalarOpacityUnitDistance = 0.006765823467065929
		calculator2Display.VolumeRenderingMode = 'Smart'
		calculator2Display.Shade = 0
		calculator2Display.SliceMode = 'XY Plane'
		calculator2Display.Slice = 127

		# hide data in view
		Hide(filename, renderView1)

		# create a new 'Contour'
		contour2 = Contour(Input=calculator2)
		contour2.ContourBy = ['POINTS', 'Result']
		contour2.ComputeNormals = 1
		contour2.ComputeGradients = 0
		contour2.ComputeScalars = 0
		contour2.GenerateTriangles = 1
		contour2.Isosurfaces = [0.0]
		contour2.PointMergeMethod = 'Uniform Binning'

		# init the 'Uniform Binning' selected for 'PointMergeMethod'
		contour2.PointMergeMethod.Divisions = [50, 50, 50]
		contour2.PointMergeMethod.Numberofpointsperbucket = 8

		# set active source
		SetActiveSource(calculator2)

		# destroy contour2
		Delete(contour2)
		del contour2

		animationScene1.GoToLast()

		# create a new 'Contour'
		contour2 = Contour(Input=calculator2)
		contour2.ContourBy = ['POINTS', 'Result']
		contour2.ComputeNormals = 1
		contour2.ComputeGradients = 0
		contour2.ComputeScalars = 0
		contour2.GenerateTriangles = 1
		contour2.Isosurfaces = [624.5271606445312]
		contour2.PointMergeMethod = 'Uniform Binning'

		# init the 'Uniform Binning' selected for 'PointMergeMethod'
		contour2.PointMergeMethod.Divisions = [50, 50, 50]
		contour2.PointMergeMethod.Numberofpointsperbucket = 8

		# Properties modified on contour2
		contour2.Isosurfaces = [0.0, 10.101010101010102, 20.202020202020204, 30.303030303030305, 40.40404040404041, 50.505050505050505, 60.60606060606061, 70.7070707070707, 80.80808080808082, 90.9090909090909, 101.01010101010101, 111.1111111111111, 121.21212121212122, 131.31313131313132, 141.4141414141414, 151.51515151515153, 161.61616161616163, 171.7171717171717, 181.8181818181818, 191.91919191919192, 202.02020202020202, 212.12121212121212, 222.2222222222222, 232.32323232323233, 242.42424242424244, 252.52525252525254, 262.62626262626264, 272.7272727272727, 282.8282828282828, 292.92929292929296, 303.03030303030306, 313.13131313131316, 323.23232323232327, 333.3333333333333, 343.4343434343434, 353.5353535353535, 363.6363636363636, 373.73737373737373, 383.83838383838383, 393.93939393939394, 404.04040404040404, 414.14141414141415, 424.24242424242425, 434.34343434343435, 444.4444444444444, 454.5454545454545, 464.64646464646466, 474.74747474747477, 484.8484848484849, 494.949494949495, 505.0505050505051, 515.1515151515151, 525.2525252525253, 535.3535353535353, 545.4545454545454, 555.5555555555555, 565.6565656565656, 575.7575757575758, 585.8585858585859, 595.959595959596, 606.0606060606061, 616.1616161616162, 626.2626262626263, 636.3636363636364, 646.4646464646465, 656.5656565656566, 666.6666666666666, 676.7676767676768, 686.8686868686868, 696.969696969697, 707.070707070707, 717.1717171717171, 727.2727272727273, 737.3737373737373, 747.4747474747475, 757.5757575757576, 767.6767676767677, 777.7777777777778, 787.8787878787879, 797.979797979798, 808.0808080808081, 818.1818181818182, 828.2828282828283, 838.3838383838383, 848.4848484848485, 858.5858585858585, 868.6868686868687, 878.7878787878788, 888.8888888888888, 898.989898989899, 909.090909090909, 919.1919191919193, 929.2929292929293, 939.3939393939395, 949.4949494949495, 959.5959595959596, 969.6969696969697, 979.7979797979798, 989.89898989899, 1000.0]

		# show data in view
		contour2Display = Show(contour2, renderView1)
		# trace defaults for the display properties.
		contour2Display.Representation = 'Surface'
		contour2Display.AmbientColor = [1.0, 1.0, 1.0]
		contour2Display.ColorArrayName = [None, '']
		contour2Display.DiffuseColor = [1.0, 1.0, 1.0]
		contour2Display.LookupTable = None
		contour2Display.MapScalars = 1
		contour2Display.InterpolateScalarsBeforeMapping = 1
		contour2Display.Opacity = 1.0
		contour2Display.PointSize = 2.0
		contour2Display.LineWidth = 1.0
		contour2Display.Interpolation = 'Gouraud'
		contour2Display.Specular = 0.0
		contour2Display.SpecularColor = [1.0, 1.0, 1.0]
		contour2Display.SpecularPower = 100.0
		contour2Display.Ambient = 0.0
		contour2Display.Diffuse = 1.0
		contour2Display.EdgeColor = [0.0, 0.0, 0.5]
		contour2Display.BackfaceRepresentation = 'Follow Frontface'
		contour2Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
		contour2Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
		contour2Display.BackfaceOpacity = 1.0
		contour2Display.Position = [0.0, 0.0, 0.0]
		contour2Display.Scale = [1.0, 1.0, 1.0]
		contour2Display.Orientation = [0.0, 0.0, 0.0]
		contour2Display.Origin = [0.0, 0.0, 0.0]
		contour2Display.Pickable = 1
		contour2Display.Texture = None
		contour2Display.Triangulate = 0
		contour2Display.NonlinearSubdivisionLevel = 1
		contour2Display.CubeAxesColor = [1.0, 1.0, 1.0]
		contour2Display.CubeAxesCornerOffset = 0.0
		contour2Display.CubeAxesFlyMode = 'Closest Triad'
		contour2Display.CubeAxesInertia = 1
		contour2Display.CubeAxesTickLocation = 'Inside'
		contour2Display.CubeAxesXAxisMinorTickVisibility = 1
		contour2Display.CubeAxesXAxisTickVisibility = 1
		contour2Display.CubeAxesXAxisVisibility = 1
		contour2Display.CubeAxesXGridLines = 0
		contour2Display.CubeAxesXTitle = 'X-Axis'
		contour2Display.CubeAxesUseDefaultXTitle = 1
		contour2Display.CubeAxesYAxisMinorTickVisibility = 1
		contour2Display.CubeAxesYAxisTickVisibility = 1
		contour2Display.CubeAxesYAxisVisibility = 1
		contour2Display.CubeAxesYGridLines = 0
		contour2Display.CubeAxesYTitle = 'Y-Axis'
		contour2Display.CubeAxesUseDefaultYTitle = 1
		contour2Display.CubeAxesZAxisMinorTickVisibility = 1
		contour2Display.CubeAxesZAxisTickVisibility = 1
		contour2Display.CubeAxesZAxisVisibility = 1
		contour2Display.CubeAxesZGridLines = 0
		contour2Display.CubeAxesZTitle = 'Z-Axis'
		contour2Display.CubeAxesUseDefaultZTitle = 1
		contour2Display.CubeAxesGridLineLocation = 'All Faces'
		contour2Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		contour2Display.CustomBoundsActive = [0, 0, 0]
		contour2Display.OriginalBoundsRangeActive = [0, 0, 0]
		contour2Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		contour2Display.CustomRangeActive = [0, 0, 0]
		contour2Display.UseAxesOrigin = 0
		contour2Display.AxesOrigin = [0.0, 0.0, 0.0]
		contour2Display.CubeAxesXLabelFormat = '%-#6.3g'
		contour2Display.CubeAxesYLabelFormat = '%-#6.3g'
		contour2Display.CubeAxesZLabelFormat = '%-#6.3g'
		contour2Display.StickyAxes = 0
		contour2Display.CenterStickyAxes = 0
		contour2Display.SelectionCellLabelBold = 0
		contour2Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
		contour2Display.SelectionCellLabelFontFamily = 'Arial'
		contour2Display.SelectionCellLabelFontSize = 18
		contour2Display.SelectionCellLabelItalic = 0
		contour2Display.SelectionCellLabelJustification = 'Left'
		contour2Display.SelectionCellLabelOpacity = 1.0
		contour2Display.SelectionCellLabelShadow = 0
		contour2Display.SelectionPointLabelBold = 0
		contour2Display.SelectionPointLabelColor = [0.5, 0.5, 0.5]
		contour2Display.SelectionPointLabelFontFamily = 'Arial'
		contour2Display.SelectionPointLabelFontSize = 18
		contour2Display.SelectionPointLabelItalic = 0
		contour2Display.SelectionPointLabelJustification = 'Left'
		contour2Display.SelectionPointLabelOpacity = 1.0
		contour2Display.SelectionPointLabelShadow = 0

		# set active source
		SetActiveSource(calculator2)

		# hide data in view
		Hide(contour2, renderView1)

		# show data in view
		calculator2Display = Show(calculator2, renderView1)

		# destroy contour2
		Delete(contour2)
		del contour2

		# create a new 'Contour'
		contour2 = Contour(Input=calculator2)
		contour2.ContourBy = ['POINTS', 'Result']
		contour2.ComputeNormals = 1
		contour2.ComputeGradients = 0
		contour2.ComputeScalars = 0
		contour2.GenerateTriangles = 1
		contour2.Isosurfaces = [624.5271606445312]
		contour2.PointMergeMethod = 'Uniform Binning'

		# init the 'Uniform Binning' selected for 'PointMergeMethod'
		contour2.PointMergeMethod.Divisions = [50, 50, 50]
		contour2.PointMergeMethod.Numberofpointsperbucket = 8

		# Properties modified on contour2
		contour2.Isosurfaces = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 79.99999999999999, 90.00000000000001, 100.0, 110.0, 119.99999999999999, 130.0, 140.0, 149.99999999999997, 160.0, 170.00000000000003, 180.0, 190.0, 200.0, 210.0, 220.0, 229.99999999999997, 240.0, 250.0, 260.0, 270.0, 279.99999999999994, 289.99999999999994, 300.0, 310.00000000000006, 320.00000000000006, 330.00000000000006, 340.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0, 410.0, 420.0, 430.0, 440.0, 449.99999999999994, 459.99999999999994, 470.0, 480.0, 490.0, 500.0, 510.00000000000006, 520.0, 530.0, 540.0, 549.9999999999999, 560.0, 569.9999999999999, 580.0, 590.0, 600.0, 610.0000000000001, 620.0, 630.0000000000001, 640.0, 650.0000000000001, 660.0, 670.0, 680.0, 690.0, 700.0, 710.0, 719.9999999999999, 730.0, 739.9999999999999, 750.0, 760.0, 770.0, 780.0, 790.0, 800.0, 810.0, 820.0000000000001, 830.0, 840.0, 850.0, 860.0, 870.0, 880.0, 889.9999999999999, 900.0, 909.9999999999999, 920.0000000000001, 930.0, 940.0000000000001, 950.0, 960.0, 970.0, 980.0, 990.0, 1000.0]

		# show data in view
		contour2Display = Show(contour2, renderView1)
		# trace defaults for the display properties.
		contour2Display.Representation = 'Surface'
		contour2Display.AmbientColor = [1.0, 1.0, 1.0]
		contour2Display.ColorArrayName = [None, '']
		contour2Display.DiffuseColor = [1.0, 1.0, 1.0]
		contour2Display.LookupTable = None
		contour2Display.MapScalars = 1
		contour2Display.InterpolateScalarsBeforeMapping = 1
		contour2Display.Opacity = 1.0
		contour2Display.PointSize = 2.0
		contour2Display.LineWidth = 1.0
		contour2Display.Interpolation = 'Gouraud'
		contour2Display.Specular = 0.0
		contour2Display.SpecularColor = [1.0, 1.0, 1.0]
		contour2Display.SpecularPower = 100.0
		contour2Display.Ambient = 0.0
		contour2Display.Diffuse = 1.0
		contour2Display.EdgeColor = [0.0, 0.0, 0.5]
		contour2Display.BackfaceRepresentation = 'Follow Frontface'
		contour2Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
		contour2Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
		contour2Display.BackfaceOpacity = 1.0
		contour2Display.Position = [0.0, 0.0, 0.0]
		contour2Display.Scale = [1.0, 1.0, 1.0]
		contour2Display.Orientation = [0.0, 0.0, 0.0]
		contour2Display.Origin = [0.0, 0.0, 0.0]
		contour2Display.Pickable = 1
		contour2Display.Texture = None
		contour2Display.Triangulate = 0
		contour2Display.NonlinearSubdivisionLevel = 1
		contour2Display.CubeAxesColor = [1.0, 1.0, 1.0]
		contour2Display.CubeAxesCornerOffset = 0.0
		contour2Display.CubeAxesFlyMode = 'Closest Triad'
		contour2Display.CubeAxesInertia = 1
		contour2Display.CubeAxesTickLocation = 'Inside'
		contour2Display.CubeAxesXAxisMinorTickVisibility = 1
		contour2Display.CubeAxesXAxisTickVisibility = 1
		contour2Display.CubeAxesXAxisVisibility = 1
		contour2Display.CubeAxesXGridLines = 0
		contour2Display.CubeAxesXTitle = 'X-Axis'
		contour2Display.CubeAxesUseDefaultXTitle = 1
		contour2Display.CubeAxesYAxisMinorTickVisibility = 1
		contour2Display.CubeAxesYAxisTickVisibility = 1
		contour2Display.CubeAxesYAxisVisibility = 1
		contour2Display.CubeAxesYGridLines = 0
		contour2Display.CubeAxesYTitle = 'Y-Axis'
		contour2Display.CubeAxesUseDefaultYTitle = 1
		contour2Display.CubeAxesZAxisMinorTickVisibility = 1
		contour2Display.CubeAxesZAxisTickVisibility = 1
		contour2Display.CubeAxesZAxisVisibility = 1
		contour2Display.CubeAxesZGridLines = 0
		contour2Display.CubeAxesZTitle = 'Z-Axis'
		contour2Display.CubeAxesUseDefaultZTitle = 1
		contour2Display.CubeAxesGridLineLocation = 'All Faces'
		contour2Display.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		contour2Display.CustomBoundsActive = [0, 0, 0]
		contour2Display.OriginalBoundsRangeActive = [0, 0, 0]
		contour2Display.CustomRange = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
		contour2Display.CustomRangeActive = [0, 0, 0]
		contour2Display.UseAxesOrigin = 0
		contour2Display.AxesOrigin = [0.0, 0.0, 0.0]
		contour2Display.CubeAxesXLabelFormat = '%-#6.3g'
		contour2Display.CubeAxesYLabelFormat = '%-#6.3g'
		contour2Display.CubeAxesZLabelFormat = '%-#6.3g'
		contour2Display.StickyAxes = 0
		contour2Display.CenterStickyAxes = 0
		contour2Display.SelectionCellLabelBold = 0
		contour2Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
		contour2Display.SelectionCellLabelFontFamily = 'Arial'
		contour2Display.SelectionCellLabelFontSize = 18
		contour2Display.SelectionCellLabelItalic = 0
		contour2Display.SelectionCellLabelJustification = 'Left'
		contour2Display.SelectionCellLabelOpacity = 1.0
		contour2Display.SelectionCellLabelShadow = 0
		contour2Display.SelectionPointLabelBold = 0
		contour2Display.SelectionPointLabelColor = [0.5, 0.5, 0.5]
		contour2Display.SelectionPointLabelFontFamily = 'Arial'
		contour2Display.SelectionPointLabelFontSize = 18
		contour2Display.SelectionPointLabelItalic = 0
		contour2Display.SelectionPointLabelJustification = 'Left'
		contour2Display.SelectionPointLabelOpacity = 1.0
		contour2Display.SelectionPointLabelShadow = 0

		# Properties modified on contour2Display
		contour2Display.Opacity = 0.3

		# set scalar coloring
		ColorBy(contour2Display, ('POINTS', 'data'))

		# rescale color and/or opacity maps used to include current data range
		contour2Display.RescaleTransferFunctionToDataRange(True)

		# show color bar/color legend
		contour2Display.SetScalarBarVisibility(renderView1, True)

		# get color transfer function/color map for 'data'
		dataLUT = GetColorTransferFunction('data')
		dataLUT.InterpretValuesAsCategories = 0
		dataLUT.EnableOpacityMapping = 0
		dataLUT.RGBPoints = [10.04987562112089, 0.231373, 0.298039, 0.752941, 514.5891600107368, 0.865003, 0.865003, 0.865003, 1019.1284444003527, 0.705882, 0.0156863, 0.14902]
		dataLUT.UseLogScale = 0
		dataLUT.LockScalarRange = 0
		dataLUT.ColorSpace = 'Diverging'
		dataLUT.UseBelowRangeColor = 0
		dataLUT.BelowRangeColor = [0.0, 0.0, 0.0]
		dataLUT.UseAboveRangeColor = 0
		dataLUT.AboveRangeColor = [1.0, 1.0, 1.0]
		dataLUT.NanColor = [0.5, 0.0, 0.0]
		dataLUT.Discretize = 1
		dataLUT.NumberOfTableValues = 256
		dataLUT.ScalarRangeInitialized = 1.0
		dataLUT.HSVWrap = 0
		dataLUT.VectorComponent = 0
		dataLUT.VectorMode = 'Magnitude'
		dataLUT.AllowDuplicateScalars = 1
		dataLUT.Annotations = []
		dataLUT.IndexedColors = []

		# get opacity transfer function/opacity map for 'data'
		dataPWF = GetOpacityTransferFunction('data')
		dataPWF.Points = [10.04987562112089, 0.0, 0.5, 0.0, 1019.1284444003527, 1.0, 0.5, 0.0]
		dataPWF.AllowDuplicateScalars = 1
		dataPWF.ScalarRangeInitialized = 1

		#change array component used for coloring
		dataLUT.RGBPoints = [-6.707550525665283, 0.231373, 0.298039, 0.752941, -2.8330719470977783, 0.865003, 0.865003, 0.865003, 1.0414066314697266, 0.705882, 0.0156863, 0.14902]
		dataLUT.VectorComponent = 2
		dataLUT.VectorMode = 'Component'

		# Properties modified on dataPWF
		dataPWF.Points = [-6.707550525665283, 0.0, 0.5, 0.0, 1.0414066314697266, 1.0, 0.5, 0.0]

		# Rescale transfer function
		dataLUT.RescaleTransferFunction(-10.0, 10.0)

		# Rescale transfer function
		dataPWF.RescaleTransferFunction(-10.0, 10.0)

		# Properties modified on dataLUT
		dataLUT.RGBPoints = [-1.0, 0.0, 1.0, 1.0, -0.09999999999999998, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.10000000000000009, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0]

		# Rescale transfer function
		dataLUT.RescaleTransferFunction(-0.1, 0.1)

		# Rescale transfer function
		dataPWF.RescaleTransferFunction(-0.1, 0.1)

		# hide color bar/color legend
		contour2Display.SetScalarBarVisibility(renderView1, False)

		# show color bar/color legend
		contour2Display.SetScalarBarVisibility(renderView1, True)

		# set active source
		SetActiveSource(filename)

		# current camera placement for renderView1
		renderView1.CameraPosition = [0.016102336606662604, 1.4945698217578052, 1.503567933991281]
		renderView1.CameraFocalPoint = [0.5975143916106522, 0.6122517855185808, 0.4516851337866846]
		renderView1.CameraViewUp = [0.030733487995079693, 0.7740952461459013, -0.6323227044871731]
		renderView1.CameraParallelScale = 0.14877755006544496
        
		renderView1.CameraPosition = [1.192582862882915, 2.193633913591958, 1.8591615301070805]
		renderView1.CameraFocalPoint = [0.4990234375000003, 0.49902343749999967, 0.49902343749999967]
		renderView1.CameraViewUp = [-0.4100779495870859, 0.6670580667347654, -0.6219884330648993]
		renderView1.CameraParallelScale = 0.8643339479176722

		dirLocal = '/scratch/daint/cconti/PVRender/'+os.path.relpath(dirName,rootDir)
		if not os.path.exists(dirLocal):
			os.makedirs(dirLocal)
		path = dirLocal+'/samara.png'
		# save animation images/movie
		WriteAnimation(path, Magnification=1, FrameRate=15.0, Compression=True)
