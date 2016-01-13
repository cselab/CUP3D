import os

#### import the simple module from the paraview
#from paraview.simple import *
#### disable automatic camera reset on 'Show'
#paraview.simple._DisableFirstRenderCameraReset()

rootDir = '/cluster/scratch_xp/public/cconti/CubismUP/'

for dirName, subDirList, fileList in os.walk(rootDir):
	if "FallingSamaraFixed_tdump_MPI4_2812" in dirName:
		print(os.path.relpath(dirName,rootDir))
		filename = [f for f in fileList if f.endswith('.xmf')]
		#print(filename)

