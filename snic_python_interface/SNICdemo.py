#--------------------------------------------------------------
# This is a demo file intended to show the use of the SNIC algorithm
# Please compile the C files of snic.h and snic.c using:
# "python snic.c" on the command prompt prior to using this file.
#
# To see the demo use: "python SNICdemo.py" on the command prompt
#------------------------------------------------------------------

import os
import subprocess
from PIL import Image
from skimage.io import imread,imshow
import numpy as np
from timeit import default_timer as timer
from _snic.lib import SNIC_main
from cffi import FFI


def segment(imgname,numsuperpixels,compactness,doRGBtoLAB):
	#--------------------------------------------------------------
	# read image and change image shape from (h,w,c) to (c,h,w)
	#--------------------------------------------------------------
	img = Image.open(imgname)
	# img = imread(imgname)
	img = np.asarray(img)
	print(img.shape)

	dims = img.shape
	h,w,c = dims[0],dims[1],1
	if len(dims) > 1:
		c = dims[2]
		img = img.transpose(2,0,1)
		print(c, "channels")
	
	#--------------------------------------------------------------
	# Reshape image to a single dimensional vector
	#--------------------------------------------------------------
	img = img.reshape(-1).astype(np.double)
	labels = np.zeros((h,w), dtype = np.int32)
	numlabels = np.zeros(1,dtype = np.int32)
	#--------------------------------------------------------------
	# Prepare the pointers to pass to the C function
	#--------------------------------------------------------------
	ffibuilder = FFI()
	pinp = ffibuilder.cast("double*", ffibuilder.from_buffer(img))
	plabels = ffibuilder.cast("int*", ffibuilder.from_buffer(labels.reshape(-1)))
	pnumlabels = ffibuilder.cast("int*", ffibuilder.from_buffer(numlabels))

	
	start = timer()
	SNIC_main(pinp,w,h,c,numsuperpixels,compactness,doRGBtoLAB,plabels,pnumlabels)
	end = timer()

	#--------------------------------------------------------------
	# Collect labels
	#--------------------------------------------------------------
	print("number of superpixels: ", numlabels[0])
	print("time taken in seconds: ", end-start)

	return labels.reshape(h,w),numlabels[0]


	# lib.SNICmain.argtypes = [np.ctypeslib.ndpointer(dtype=POINTER(c_double),ndim=2)]+[c_int]*4 +[c_double,c_bool,ctypes.data_as(POINTER(c_int)),ctypes.data_as(POINTER(c_int))]

def drawBoundaries(imgname,labels,numlabels):

	img = Image.open(imgname)
	# img = imread(imgname)
	img = np.array(img)
	print(img.shape)

	ht,wd = labels.shape

	for y in range(1,ht-1):
		for x in range(1,wd-1):
			if labels[y,x-1] != labels[y,x+1] or labels[y-1,x] != labels[y+1,x]:
				img[y,x,:] = 0

	return img
	
# Before calling this function, please compile the C code using
# "python compile.py" on the command line
def snicdemo():
	#--------------------------------------------------------------
	# Set parameters and call the C function
	#--------------------------------------------------------------
	numsuperpixels = 500
	compactness = 20.0
	doRGBtoLAB = True # only works if it is a three channel image
	# imgname = "/Users/achanta/Pictures/classics/lena.png"
	imgname = "bee.png"
	labels,numlabels = segment(imgname,numsuperpixels,compactness,doRGBtoLAB)
	#--------------------------------------------------------------
	# Display segmentation result
	#------------------------------------------------------------
	segimg = drawBoundaries(imgname,labels,numlabels)
	# Image.fromarray(segimg).show()
	Image.fromarray(segimg).save("bee_snic.png")
	return

snicdemo()



