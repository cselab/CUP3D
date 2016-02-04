import sys
import os

n = int(sys.argv[1])

f = open("file_list.txt", "w")
wd = os.getcwd()

pdir = sys.argv[2]
gdir = sys.argv[3]

for i in range(n+1):
    f.write(wd + "/" + pdir + "/p%05d.vpcache\n" % (i))
    f.write(wd + "/" + gdir + "/g%05d.vpcache\n" % (i))
f.close()
