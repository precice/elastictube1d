import vtk
import numpy as np
import os
import netCDF4 as nc

import configuration_file as config
import tubePlotting

import matplotlib.pyplot as plt


def numpyDataToVTKPointData(grid, numpy_data, dataname):

    if numpy_data.shape.__len__() == 1:
        new_shape = np.array([numpy_data.shape[0],1,1])
        new_data = np.ones(new_shape)
        new_data[:,0,0] = numpy_data[:]
        numpy_data = new_data
    elif numpy_data != 3:
        print "cannot handle shape = %i" % (numpy_data.shape)
        quit()

    # add point dataset
    vtk_array = vtk.vtkDoubleArray()
    vtk_array.SetNumberOfComponents(1)
    vtk_array.SetNumberOfTuples(grid.GetNumberOfPoints())
    vtk_array.SetName(dataname)
    for z_id in range(grid.GetDimensions()[2]):
        for y_id in range(grid.GetDimensions()[1]):
            for x_id in range(grid.GetDimensions()[0]):
                id = grid.ComputePointId((x_id , y_id , z_id))
                vtk_array.SetValue(id, numpy_data[x_id, y_id, z_id])

    return vtk_array


def numpyDataToVTKCellData(grid, numpy_data, dataname):

    if numpy_data.shape.__len__() == 1:
        new_shape = np.array([numpy_data.shape[0],1,1])
        new_data = np.ones(new_shape)
        new_data[:,0,0] = numpy_data[:]
        numpy_data = new_data
    elif numpy_data != 3:
        print "cannot handle shape = %i" % (numpy_data.shape)
        quit()

    # add cell dataset
    vtk_array = vtk.vtkDoubleArray()
    vtk_array.SetNumberOfComponents(1)
    vtk_array.SetNumberOfTuples(grid.GetNumberOfCells())
    vtk_array.SetName(dataname)
    for z_id in range(grid.GetDimensions()[2]-1):
        for y_id in range(grid.GetDimensions()[1]-1):
            for x_id in range(grid.GetDimensions()[0]-1):
                id = grid.ComputeCellId((x_id , y_id , z_id))
                vtk_array.SetValue(id, numpy_data[x_id, y_id, z_id])

    return vtk_array


def writeOutputToVTK(time, xpos, filename, datasets, datanames, metadata = dict()):

    writer = vtk.vtkStructuredPointsWriter()

    if type(datasets) is not list:
        datasets = list(datasets)
    if type(datanames) is not list:
        datanames = list(datanames)

    n_datasets = datasets.__len__()
    assert n_datasets == datanames.__len__()

    dy = dz = 0
    ny = nz = 1

    outpath = os.path.join(os.getcwd(), 'VTK')

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    filename = filename + "." + str(time) + ".vtk"
    filepath = os.path.join(outpath , filename)

    dx = xpos[1] - xpos[0]
    nx = xpos.shape[0]

    # initialize  vtk  grid
    grid = vtk.vtkImageData()
    grid.SetOrigin(0, 0, 0)
    grid.SetSpacing(dx, dy, dz)
    grid.SetDimensions(nx+1, ny+1, nz+1)

    for i in range(n_datasets):
        numpy_data = datasets[i]
        dataname = datanames[i]
        vtk_array = numpyDataToVTKCellData(grid, numpy_data, dataname)
        grid.GetCellData().AddArray(vtk_array)

    writer.SetInputData(grid)
    writer.SetFileName(filepath)
    writer.Write()


def writeOutputToNetCDF(time, xpos, filename, datasets, datanames, metadata = dict()):

    if type(datasets) is not list:
        datasets = list(datasets)
    if type(datanames) is not list:
        datanames = list(datanames)

    n_datasets = datasets.__len__()
    assert n_datasets == datanames.__len__()

    outpath = os.path.join(os.getcwd(), 'NCDF')

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    filepath = os.path.join(outpath, filename+'.nc')

    if not os.path.isfile(filepath):  # file does not exist -> create new
        nc_dataset = nc.Dataset(filepath, 'w', format='NETCDF4_CLASSIC')
    else:  # file already exists -> append
        nc_dataset = nc.Dataset(filepath, 'r+', format='NETCDF4_CLASSIC')

    if not nc_dataset.dimensions.has_key('X'):
        nc_dataset.createDimension('X', xpos.shape[0])
    if not nc_dataset.dimensions.has_key('time'):
        nc_dataset.createDimension('time', None)

    if not nc_dataset.variables.has_key('xcoords'):
        nc_dataset.createVariable('xcoords', np.float64, ('X'))
        nc_dataset.variables['xcoords'][:] = xpos

    if not nc_dataset.variables.has_key('times'):
        nc_dataset.createVariable('times', np.float64, ('time'))

    step_id = nc_dataset.variables['times'].__len__()
    nc_dataset.variables['times'][step_id] = time

    for i_data in range(n_datasets):
        dataname = datanames[i_data]
        dataset = datasets[i_data]
        if not nc_dataset.variables.has_key(dataname):
            nc_dataset.createVariable(dataname, np.float64, ('X', 'time'))

        nc_dataset.variables[dataname][:,step_id] = dataset[:]

    for key in metadata.keys():
        if not hasattr(nc_dataset, key):
            nc_dataset.setncattr(str(key), metadata[key])

    nc_dataset.close()


def create_output(t, velocity, pressure, crossSectionLength, metadata, output_mode):
    dx = metadata["dx"]
    N = metadata["n_elem"]

    if output_mode is not config.OutputModes.OFF:
        x = np.linspace(0,dx*N,N+1,endpoint=True)
        filename = "fluid_"+str(metadata["created_on"])
        if output_mode is config.OutputModes.VTK:
            writeOutputToVTK(t, x, filename, datanames=["velocity", "pressure", "crossSection"], datasets=[velocity, pressure, crossSectionLength])
        if output_mode is config.OutputModes.NETCDF:
            writeOutputToNetCDF(t, x, filename, datanames=["velocity", "pressure", "crossSection"], datasets=[velocity, pressure, crossSectionLength], metadata=metadata)

def create_video(t, velocity, pressure, crossSectionLength, metadata, writer):
    tubePlotting.doPlotting(plt.gca(), crossSectionLength, velocity, pressure, metadata['dx'], t)
    writer.grab_frame()
    plt.gca().cla()
