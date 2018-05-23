from lxml import etree as ET
import numpy as np
import os
import configuration_file as config

template = "helpers/precice-config.xml"
basepath = os.path.dirname(template)
basename = os.path.basename(template)

tree = ET.parse(template)

## make sure that template and configuration_file.py are consistent

assert np.abs(float(tree.findall('**/max-time')[0].attrib['value']) - config.T_max) < 10**-10

## vary timestep length and generate multiple precice-config.xml from template

varied_variable = "timestep-length"

tau0 = config.tau0
N_experiments = config.n_tau
varied_values = tau0 * .5**np.arange(0,N_experiments)

elem = tree.findall('**/'+varied_variable)[0]

outfolder = 'configs'

if os.path.exists(outfolder):
    for file in os.listdir(outfolder):
        os.remove(os.path.join(outfolder,file))
    os.rmdir(outfolder)

os.mkdir(outfolder)

for value in varied_values:
    elem.attrib['value'] = str(value)
    newname = os.path.splitext(basename)[0]+"_"+str(value).replace('.','_')+".xml"
    tree.write(os.path.join(outfolder,newname))
