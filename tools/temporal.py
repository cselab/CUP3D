import glob

def ReadAppend(name,lines,lines2):
    f = open(name,'r')
    data = f.read()
    f.close()
    data = data.split("\n",lines)[lines]
    data = ("\n".join(data.split("\n")[:-lines2]))
    return data

def CreateTemporal(name):
    all_files = glob.glob(name+"*.xmf")
    all_files.sort()

    mystring  = "<?xml version=\"1.0\" ?>\n"
    mystring += "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
    mystring += "<Xdmf Version=\"2.0\">\n"
    mystring += "<Domain>\n"
    mystring += "  <Grid Name=\"" + name + "--temporal" + "\" GridType=\"Collection\" CollectionType=\"Temporal\" >\n"
    k = 0
    for i in all_files:
        mystring += "\n\n\n"
        mystring += "<Grid GridType=\"Uniform\">\n"
        mystring += "   <Time Value=\""  + str(k)  +  "\"/>\n"
        mystring += ReadAppend(i,7,3)
        mystring += "\n\n\n"
        k = k + 1
    mystring += " </Grid>\n"
    mystring += "</Domain>\n"
    mystring += "</Xdmf>\n"

    f = open(name + "-temporalCollection.xmf", 'a+')
    f.write(mystring) 
    f.close()

if __name__ == "__main__":
    CreateTemporal("chi")
    CreateTemporal("tmp")
