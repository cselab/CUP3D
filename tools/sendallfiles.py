import os
import glob
dst='michaich@hippo.ethz.ch:/data/users/michaich/'
parts = 128
res = []
for path in os.listdir("."):
    if os.path.isfile(os.path.join(".", path)):
        res.append(path)
step = int (len(res)/parts)
number = 0
for i in range(0,len(res),step):
    number += 1
    s = " "
    smax = min(i+step,len(res))
    for j in range(i,smax):
        s += res[j] + " "
    os.system("nohup rsync -r -v --progress --ignore-existing " + s + dst + " > syncfile" + str(number) + " &")
    os.system("sleep 1")
