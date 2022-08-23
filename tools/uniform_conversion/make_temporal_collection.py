import glob, string, sys
import h5py
import re

baseDirectory = "./"
dt = 0.50
h  = 0.001953 * 4
files1 = glob.glob(baseDirectory + "chi_*-uniform.h5")
files1.sort()

# open first file to inspect dimensions
Cfile = h5py.File(files1[0], "r")
dims = Cfile["/data"].shape
Cfile.close()

# for the origin, we add a half-cell length h/2
#
#open new file to dump all timesteps descriptions
f2 = open("all_nodedata_files.xmf", "w")
f2.write('<?xml version="1.0" ?>\n')
f2.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
f2.write('<Xdmf Version="2.0">\n')
f2.write('<Domain>\n')
f2.write('    <Topology TopologyType="3DCoRectMesh"\n')
f2.write(format("        Dimensions=\"%d %d %d\">\n" % dims))
f2.write('    </Topology>\n')
f2.write('    <Geometry Type="ORIGIN_DXDYDZ">\n')
f2.write(format("        <DataItem Dimensions=\"3\" NumberType=\"Double\" Precision=\"8\" Format=\"XML\">\n"))
f2.write(format("        %s %s %s \n" % (0.5*h,0.5*h,0.5*h)))
f2.write('        </DataItem>\n')
f2.write(format("        <DataItem Dimensions=\"3\" NumberType=\"Double\" Precision=\"8\" Format=\"XML\">\n"))
f2.write(format("        %s %s %s \n" % (h,h,h)))
f2.write('        </DataItem>\n')
f2.write('    </Geometry>\n')
f2.write('\n')
f2.write('    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n')
f2.write('      <Time TimeType="List">\n')
f2.write('        <DataItem Format="XML" NumberType="Float" Dimensions="')
f2.write(str(len(files1)))
f2.write('">\n')
for t in range(len(files1)):
  f2.write(format("%s " % (t*dt)))
f2.write('\n        </DataItem>\n')
f2.write('      </Time>\n')
for f in files1:
  f2.write('    <Grid Name="TS" GridType="Uniform">\n')
  f2.write('      <Topology Reference="/Xdmf/Domain/Topology[1]"/>\n')
  f2.write('      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>\n')
  fdata = f[:-3] + '.h5'
  fdata2 = f[:-3].replace("chi", "tmp") + '.h5'
  f2.write('      <Attribute Name="' + 'chi' + '" AttributeType="Scalar" Center="Node">\n')
  f2.write(format("        <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n" % (dims)))
  f2.write('        ' + fdata + ':/data\n')
  f2.write('        </DataItem>\n')
  f2.write('      </Attribute>\n')
#  f2.write('      <Attribute Name="' + 'vorticity_magnitude' + '" AttributeType="Scalar" Center="Node">\n')
#  f2.write(format("        <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n" % (dims)))
#  f2.write('        ' + fdata2 + ':/data\n')
#  f2.write('        </DataItem>\n')
#  f2.write('      </Attribute>\n')
  f2.write('    </Grid>\n')
f2.write('\n    </Grid>\n  </Domain>\n</Xdmf>\n')
f2.close()
