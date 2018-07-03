from lxml import etree as ET
import numpy as np
import os
import configuration_file as config

template = "helpers/precice-config.xml"  # template precice-config.xml used to generate series of configs
basepath = os.path.dirname(template)
basename = os.path.basename(template)

tree = ET.parse(template)

## make sure that template and configuration_file.py are consistent

assert np.abs(float(tree.findall('**/max-time')[0].attrib['value']) - config.T_max) < 10**-10

## vary timestep length and generate multiple precice-config.xml from template

varied_variable = "timestep-length"  # corresponding tag in precice-config.xml template
elem = tree.findall('**/'+varied_variable)[0]  # find tag in template

## create list holding the values that you want to use for generating series of configs
tau0 = config.tau0  # get smallest timestep length from python config
N_experiments = config.n_tau
varied_values = tau0 * .5**np.arange(0,N_experiments)  # create list of considered timestep sizes

## create folder to write output to 
outfolder = 'configs'
if os.path.exists(outfolder):
    for file in os.listdir(outfolder):
        os.remove(os.path.join(outfolder,file))
    os.rmdir(outfolder)
os.mkdir(outfolder)

## iterate through list varied_values and create a precice-config.xml based on the given template with modified varied_variable
for value in varied_values:
    elem.attrib['value'] = str(value)
    newname = os.path.splitext(basename)[0]+"_"+str(value).replace('.','_')+".xml"
    tree.write(os.path.join(outfolder,newname))
