import sys
import os
import datetime

import xarray as xr
import subprocess
import numpy as np

from os.path import basename
from osgeo import gdal, osr, gdal_array, gdalnumeric
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *

def files_mean(file_array, filename_out):
    myFiles=[]
    myNDV=[]
    DimensionsCheck=None
    allBandsCount=1
    myOutNDV=-32767

    for (i, filename) in enumerate(file_array):
        if os.path.isfile(filename):
            myFile = gdal.Open(filename, gdal.GA_ReadOnly)
            if not myFile:
                raise IOError("No such file or directory: '%s'" % myFile)
            myFiles.append(myFile)
            myNDV.append(myFile.GetRasterBand(1).GetNoDataValue())
            if DimensionsCheck:
                if DimensionsCheck != [myFile.RasterXSize, myFile.RasterYSize]:
                    raise Exception("Error! Dimensions of file %s (%i, %i) are different from other files (%i, %i).  Cannot proceed" % (myF, myFile.RasterXSize, myFile.RasterYSize, DimensionsCheck[0], DimensionsCheck[1]))
            else:
                DimensionsCheck = [myFile.RasterXSize, myFile.RasterYSize]

    if os.path.isfile(filename_out):
        os.remove(filename_out)

    myBlockSize=myFiles[0].GetRasterBand(1).GetBlockSize();
    # store these numbers in variables that may change later
    nXValid = myBlockSize[0]
    nYValid = myBlockSize[1]
    # find total x and y blocks to be read
    nXBlocks = (int)((DimensionsCheck[0] + myBlockSize[0] - 1) / myBlockSize[0]);
    nYBlocks = (int)((DimensionsCheck[1] + myBlockSize[1] - 1) / myBlockSize[1]);
    myBufSize = myBlockSize[0]*myBlockSize[1]

    myOutDrv = gdal.GetDriverByName('GTiff')
    myOut = myOutDrv.Create(filename_out, DimensionsCheck[0], DimensionsCheck[1], 1,gdal.GDT_Float32)

    # set output geo info based on first input layer
    myOut.SetGeoTransform(myFiles[0].GetGeoTransform())
    myOut.SetProjection(myFiles[0].GetProjection())

    for bandNo in range(1,allBandsCount+1):
        for X in range(0,nXBlocks):
            if X==nXBlocks-1:
                nXValid = DimensionsCheck[0] - X * myBlockSize[0]
                myBufSize = nXValid*nYValid

            # find X offset
            myX=X*myBlockSize[0]

            # reset buffer size for start of Y loop
            nYValid = myBlockSize[1]
            myBufSize = nXValid*nYValid

            # loop through Y lines
            for Y in range(0,nYBlocks):
                # change the block size of the final piece
                if Y==nYBlocks-1:
                    nYValid = DimensionsCheck[1] - Y * myBlockSize[1]
                    myBufSize = nXValid*nYValid
                # find Y offset
                myY=Y*myBlockSize[1]
                # create empty buffer to mark where nodata occurs
                myNDVs = None
                # make local namespace for calculation
                local_namespace = []
                for (i, myFile) in enumerate(myFiles):
                    myval=gdalnumeric.BandReadAsArray(myFile.GetRasterBand(1),
                                            xoff=myX, yoff=myY,
                                            win_xsize=nXValid, win_ysize=nYValid)
                    # fill in nodata values
                    if myNDV[i] is not None:
                        if myNDVs is None:
                            myNDVs = numpy.zeros(myBufSize)
                            myNDVs.shape=(nYValid,nXValid)
                        myNDVs=1*numpy.logical_or(myNDVs==1, myval==myNDV[i])

                    local_namespace.append(myval)
                    myval=None

                # calculate mean between bands using Numpy
                myResult = numpy.mean(local_namespace, axis=0)

                if myNDVs is not None:
                    myResult = ((1*(myNDVs==0))*myResult) + (myOutNDV*myNDVs)
                elif not isinstance(myResult, numpy.ndarray):
                    myResult = numpy.ones( (nYValid,nXValid) ) * myResult

                myOutB=myOut.GetRasterBand(bandNo)
                gdalnumeric.BandWriteArray(myOutB, myResult, xoff=myX, yoff=myY)
    return

    # how to use it
    files_mean(['./input_01.tif', './input_02.tif', './input_02.tif', ...  './input_n.tif',],'./output.tif')