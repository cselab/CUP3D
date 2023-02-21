import argparse
from paraview.simple import *
parser = argparse.ArgumentParser()
parser.add_argument('--frames', required=True, type=int)
args = vars(parser.parse_args())
frameStart = args['frames']

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XDMF Reader'
cf4Re16000xmf = XDMFReader(registrationName='cf8.xmf', FileNames=['all_nodedata_files.xmf'])
cf4Re16000xmf.PointArrayStatus = ['chi']#, 'vorticity_magnitude']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
cf4Re16000xmfDisplay = Show(cf4Re16000xmf, renderView1, 'UniformGridRepresentation')

# trace defaults for the display properties.
cf4Re16000xmfDisplay.Representation = 'Outline'
cf4Re16000xmfDisplay.ColorArrayName = ['POINTS', '']
cf4Re16000xmfDisplay.SelectTCoordArray = 'None'
cf4Re16000xmfDisplay.SelectNormalArray = 'None'
cf4Re16000xmfDisplay.SelectTangentArray = 'None'
cf4Re16000xmfDisplay.OSPRayScaleArray = 'chi'
cf4Re16000xmfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
cf4Re16000xmfDisplay.SelectOrientationVectors = 'None'
cf4Re16000xmfDisplay.ScaleFactor = 0.7996093750000001
cf4Re16000xmfDisplay.SelectScaleArray = 'chi'
cf4Re16000xmfDisplay.GlyphType = 'Arrow'
cf4Re16000xmfDisplay.GlyphTableIndexArray = 'chi'
cf4Re16000xmfDisplay.GaussianRadius = 0.03998046875
cf4Re16000xmfDisplay.SetScaleArray = ['POINTS', 'chi']
cf4Re16000xmfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
cf4Re16000xmfDisplay.OpacityArray = ['POINTS', 'chi']
cf4Re16000xmfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
cf4Re16000xmfDisplay.DataAxesGrid = 'GridAxesRepresentation'
cf4Re16000xmfDisplay.PolarAxes = 'PolarAxesRepresentation'
cf4Re16000xmfDisplay.ScalarOpacityUnitDistance = 0.007595617541536202
cf4Re16000xmfDisplay.OpacityArrayName = ['POINTS', 'chi']
cf4Re16000xmfDisplay.IsosurfaceValues = [0.1]
cf4Re16000xmfDisplay.SliceFunction = 'Plane'
cf4Re16000xmfDisplay.Slice = 511

# init the 'Plane' selected for 'SliceFunction'
cf4Re16000xmfDisplay.SliceFunction.Origin = [4.0, 2.0, 2.0]

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=cf4Re16000xmf)
contour1.ContourBy = ['POINTS', 'chi']
contour1.Isosurfaces = [0.1]
contour1.PointMergeMethod = 'Uniform Binning'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'chi'
chiLUT = GetColorTransferFunction('chi')

# trace defaults for the display properties.
contour1Display.Representation = 'Surface'
contour1Display.ColorArrayName = ['POINTS', 'chi']
contour1Display.LookupTable = chiLUT
contour1Display.SelectTCoordArray = 'None'
contour1Display.SelectNormalArray = 'Normals'
contour1Display.SelectTangentArray = 'None'
contour1Display.OSPRayScaleArray = 'chi'
contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
contour1Display.SelectOrientationVectors = 'None'
contour1Display.ScaleFactor = 0.39654275178909304
contour1Display.SelectScaleArray = 'chi'
contour1Display.GlyphType = 'Arrow'
contour1Display.GlyphTableIndexArray = 'chi'
contour1Display.GaussianRadius = 0.01982713758945465
contour1Display.SetScaleArray = ['POINTS', 'chi']
contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
contour1Display.OpacityArray = ['POINTS', 'chi']
contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
contour1Display.DataAxesGrid = 'GridAxesRepresentation'
contour1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
contour1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
contour1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# show color bar/color legend
contour1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'chi'
chiPWF = GetOpacityTransferFunction('chi')

# create a new 'Smooth'
smooth1 = Smooth(registrationName='Smooth1', Input=contour1)

# Properties modified on smooth1
smooth1.NumberofIterations = 10

# show data in view
smooth1Display = Show(smooth1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
smooth1Display.Representation = 'Surface'
smooth1Display.ColorArrayName = ['POINTS', 'chi']
smooth1Display.LookupTable = chiLUT
smooth1Display.SelectTCoordArray = 'None'
smooth1Display.SelectNormalArray = 'Normals'
smooth1Display.SelectTangentArray = 'None'
smooth1Display.OSPRayScaleArray = 'chi'
smooth1Display.OSPRayScaleFunction = 'PiecewiseFunction'
smooth1Display.SelectOrientationVectors = 'None'
smooth1Display.ScaleFactor = 0.3954362988471985
smooth1Display.SelectScaleArray = 'chi'
smooth1Display.GlyphType = 'Arrow'
smooth1Display.GlyphTableIndexArray = 'chi'
smooth1Display.GaussianRadius = 0.019771814942359925
smooth1Display.SetScaleArray = ['POINTS', 'chi']
smooth1Display.ScaleTransferFunction = 'PiecewiseFunction'
smooth1Display.OpacityArray = ['POINTS', 'chi']
smooth1Display.OpacityTransferFunction = 'PiecewiseFunction'
smooth1Display.DataAxesGrid = 'GridAxesRepresentation'
smooth1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
smooth1Display.ScaleTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
smooth1Display.OpacityTransferFunction.Points = [0.5, 0.0, 0.5, 0.0, 0.5001220703125, 1.0, 0.5, 0.0]

# hide data in view
Hide(contour1, renderView1)

# show color bar/color legend
smooth1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# turn off scalar coloring
ColorBy(smooth1Display, None)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(chiLUT, renderView1)

# Properties modified on smooth1Display
smooth1Display.Specular = 1.0
#smooth1Display.Opacity = 0.5

# Properties modified on smooth1Display
smooth1Display.Luminosity = 60.0

# change solid color
smooth1Display.AmbientColor = [0.8, 0.8, 0.8]
smooth1Display.DiffuseColor = [0.8, 0.8, 0.8]

'''
# set active source
SetActiveSource(cf4Re16000xmf)
ColorBy(cf4Re16000xmfDisplay, ('POINTS', 'vorticity_magnitude'))
cf4Re16000xmfDisplay.RescaleTransferFunctionToDataRange(True, False)
vorticity_magnitudeLUT = GetColorTransferFunction('vorticity_magnitude')
vorticity_magnitudePWF = GetOpacityTransferFunction('vorticity_magnitude')
vorticity_magnitudeLUT.RescaleTransferFunction(0.0, 1.0)
vorticity_magnitudePWF.RescaleTransferFunction(0.0, 1.0)
cf4Re16000xmfDisplay.SetRepresentationType('Volume')
cf4Re16000xmfDisplay.ScalarOpacityUnitDistance = 0.01
cf4Re16000xmfDisplay.SetScalarBarVisibility(renderView1, False)
#vorticity_magnitudeLUT.ApplyPreset('Yellow - Gray - Blue', True)
#vorticity_magnitudePWF.ApplyPreset('Yellow - Gray - Blue', True)
vorticity_magnitudeLUT.RGBPoints = 		[
			2.6010972889456418e-11,
			0.231373,
			0.298039,
			0.75294099999999997,
			2.4765881388036952,
			0.0,
			0.60784313725490191,
			0.9137254901960784,
			6.0,
			0.81176470588235294,
			0.81176470588235294,
			0.40392156862745099,
			55.0,
			0.77254901960784317,
			0.23529411764705882,
			0.23529411764705882,
			79.971961975097656,
			0.54509803921568623,
			0.011764705882352941,
			0.11764705882352941
		]
vorticity_magnitudeLUT.ColorSpace = 'RGB'
vorticity_magnitudeLUT.ScalarRangeInitialized = 1.0
vorticity_magnitudePWF.Points = 		[
			2.6010972889456418e-11,
			0.0,
			0.5,
			0.0,
			79.971961975097656,
			1.0,
			0.5,
			0.0
		]
'''

LoadPalette(paletteName='BlackBackground')
#LoadPalette(paletteName='WhiteBackground')
# Hide orientation axes
renderView1.OrientationAxesVisibility = 0


# get camera animation track for the view
Cues = GetCameraTrack(view=renderView1)
Cues.Mode =  'Interpolate Camera'
Cues.Interpolation =  'Spline'
k00 = CameraKeyFrame()
k00.KeyTime =  0.0
k00.Position =  [7.27744, 5.39281, 3.25396]
k00.FocalPoint =  [6.72919, 4.89184, 3.09188]
k00.ViewUp =  [0.0, 0.0, 1.0]
k00.ViewAngle =  30.0
k00.ParallelScale =  1.0
k01 = CameraKeyFrame()
k01.KeyTime =  0.25
k01.Position =  [4.77413, 4.24098, 2.61751]
k01.FocalPoint =  [4.39389, 3.59783, 2.47752]
k01.ViewUp =  [-0.138711, -0.131634, 0.981546]
k01.ViewAngle =  30.0
k01.ParallelScale =  1.0
k02 = CameraKeyFrame()
k02.KeyTime =  0.5
k02.Position =  [2.00452, 4.64707, 2.28156]
k02.FocalPoint =  [2.36701, 3.98134, 2.22473]
k02.ViewUp =  [0.0586561, -0.0531575, 0.996862]
k02.ViewAngle =  30.0
k02.ParallelScale =  1.0
k03 = CameraKeyFrame()
k03.KeyTime =  0.75
k03.Position =  [1.72864, 2.80101, 2.15553]
k03.FocalPoint =  [2.44399, 2.54669, 2.11775]
k03.ViewUp =  [0.0440094, -0.0245645, 0.998729]
k03.ViewAngle =  30.0
k03.ParallelScale =  1.0
k04 = CameraKeyFrame()
k04.KeyTime =  1.0
k04.Position =  [1.811, 2.14885, 2.38225]
k04.FocalPoint =  [2.54601, 2.23332, 2.20779]
k04.ViewUp =  [0.230079, 0.00788421, 0.97314]
k04.ViewAngle =  30.0
k04.ParallelScale =  1.0
Cues.KeyFrames = [k00,k01,k02,k03,k04]
	

# save animation
SaveAnimation('./fish.png', renderView1, ImageResolution=[1920, 1080], FrameWindow=[frameStart, frameStart+9])
