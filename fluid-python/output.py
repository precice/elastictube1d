from __future__ import print_function
import vtk
import numpy as np
import os

writer = vtk.vtkStructuredPointsWriter()

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
        print("cannot handle shape = %i" % (numpy_data.shape))
        quit()

    # add point dataset
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


def writeOutputToVTK1(time, name, dx, nx, data, datanames):

    if type(data) is not list:
        data = list(data)
    if type(datanames) is not list:
        datanames = list(datanames)

    n_datasets = data.__len__()
    assert n_datasets == datanames.__len__()

    dy = dz = 0
    ny = nz = 1

    outpath = os.path.join(os.getcwd(), 'Postproc')

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    filename = name+"."+str(time)+".vtk"
    filepath = os.path.join(outpath , filename)

    # initialize  vtk  grid
    grid = vtk.vtkImageData()
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

def writeOutputToVTK(time, name, dx, nx, data, datanames):

    if type(data) is not list:
        data = list(data)
    if type(datanames) is not list:
        datanames = list(datanames)

    n_datasets = data.__len__()
    assert n_datasets == datanames.__len__()

    dy = dz = 0
    ny = nz = 1

    outpath = os.path.join(os.getcwd(), 'Postproc')

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    filename = name+str(time)+".vtk"
    filepath = os.path.join(outpath , filename)

    i = 0

    f=open(filepath,'w')

    f.write("# vtk DataFile Version 2.0")
    f.write("\n")
    f.write("\n")

    f.write("ASCII")
    f.write("\n")
    f.write("\n")

    f.write("DATASET UNSTRUCTURED_GRID")
    f.write("\n")
    f.write("\n")

    f.write("POINTS ")
    f.write(str(len(data[i])))
    f.write (" float")
    f.write("\n")
    f.write("\n")

    for k in range(len(data[i])):
        f.write(str("{:.16e}".format(0.0+ k*dx)))
        f.write(" 0.0000000000000000e+00 0.0000000000000000e+00")
        f.write("\n")
    f.write("\n")

    f.write("POINT_DATA ")
    f.write(str(len(data[i])))
    f.write("\n")
    f.write("\n")

    for dataname in datanames:

        if (i==0):
            f.write("VECTORS ")
        else:
            f.write("SCALARS ")
           
        f.write(dataname)
        f.write(" float")
        f.write("\n")
        if (i!=0):
            f.write("LOOKUP_TABLE default")
            f.write("\n")
        for element in data[i]:
            f.write(str("{:.16e}".format(element)))
            if (i==0):
                f.write(" 0.0000000000000000e+00 0.0000000000000000e+00")
            f.write("\n")
        f.write("\n")
        f.write("\n")
        i = i+1
    f.close()