from __future__ import print_function
from vtk import vtkStructuredPointsWriter, vtkDoubleArray, vtkImageData
import numpy as np
import os

writer = vtkStructuredPointsWriter()

def numpyDataToVTKPointData(grid, numpy_data, dataname):

    if numpy_data.shape.__len__() == 1:
        new_shape = np.array([numpy_data.shape[0],1,1])
        new_data = np.ones(new_shape)
        new_data[:,0,0] = numpy_data[:]
        numpy_data = new_data
    elif numpy_data != 3:
        print("cannot handle shape = %i" % (numpy_data.shape))
        quit()

    # add point dataset
    vtk_array = vtkDoubleArray()
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
        print("cannot handle shape = %i" % (numpy_data.shape))
        quit()

    # add point dataset
    vtk_array = vtkDoubleArray()
    vtk_array.SetNumberOfComponents(1)
    vtk_array.SetNumberOfTuples(grid.GetNumberOfCells())
    vtk_array.SetName(dataname)
    for z_id in range(grid.GetDimensions()[2]-1):
        for y_id in range(grid.GetDimensions()[1]-1):
            for x_id in range(grid.GetDimensions()[0]-1):
                id = grid.ComputeCellId((x_id , y_id , z_id))
                vtk_array.SetValue(id, numpy_data[x_id, y_id, z_id])

    return vtk_array


def writeOutputToVTK(time, name, dx, nx, data, datanames):

    if type(data) is not list:
        data = list(data)
    if type(datanames) is not list:
        datanames = list(datanames)

    n_datasets = data.__len__()
    assert n_datasets == datanames.__len__()

    dy = dz = 0
    ny = nz = 1

    outpath = os.path.join(os.getcwd(), 'VTK')

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    filename = "{name}.{time:.2f}.vtk".format(name=name,time=time)
    # filename = name+"."+str(time)+".vtk"
    filepath = os.path.join(outpath , filename)

    # initialize  vtk  grid
    grid = vtkImageData()
    grid.SetOrigin(0, 0, 0)
    grid.SetSpacing(dx, dy, dz)
    grid.SetDimensions(nx+1, ny+1, nz+1)

    for i in range(n_datasets):
        numpy_data = data[i]
        dataname = datanames[i]
        vtk_array = numpyDataToVTKCellData(grid, numpy_data, dataname)
        grid.GetCellData().AddArray(vtk_array)

    writer.SetInputData(grid)
    writer.SetFileName(filepath)
    writer.Write()
