from lxml import etree as ET
import numpy as np
import os

baseconfig = "helpers/precice-config.xml"
basepath = os.path.dirname(baseconfig)
basename = os.path.basename(baseconfig)

varied_variable = "timestep-length"

tau0 = 1
varied_values = tau0 * .5**np.arange(0,5)

tree = ET.parse(baseconfig)

elem = tree.findall('**/'+varied_variable)[0]

outfolder = 'configs'

if os.path.exists(outfolder):
    os.rmdir(outfolder)

os.mkdir(outfolder)

for value in varied_values:
    elem.attrib['value'] = str(value)
    newname = os.path.splitext(basename)[0]+"_"+str(value).replace('.','_')+".xml"
    tree.write(os.path.join(outfolder,newname))
