

import os
import sys
import gdal, osr, ogr

sitename='sa1'
path='C:\\sync\\wrk\\Appia0216\\GeoRes\\sa1'
pathTif=path+ '\\tif\\'
pathCsv=path+'\\raw\\'

def getSites():
	import os, sys
	dir= 'C:\\sync\\wrk\\Appia0216\\GeoRes\\'
	sites = os.listdir(dir)
	if (len(sites)>0):
		for i in range (0,len(sites)):
			print("%s\t%s" %(i+1, sites[i]))
		while(True):
			siteSelect=input()
			if siteSelect in range (1, len(sites)+1):
				path=dir+sites[siteSelect-1]
				print (sites[siteSelect-1]+" selected." + path)
				pathTif= path + '\\tif\\'
				pathCsv= path + '\\raw\\'
				break
			else:
				print('Enter valid number')
				continue;
	else:
		print('No site dircetory found')
	return(path, pathTif, pathCsv)
	
def downloadRM85():
	import os.path
	import serial
	ser = serial.Serial('COM3', 9600, timeout=1) #Tried with and without the last 3 parameters, and also at 1Mbps, same happens.
	ser.flushInput()
	ser.flushOutput()
	overwrite="na"
	init=int(raw_input('first grid number:'))
	print("Waiting for data...")
	started = 0
	ii=init-1
	data_raw=""
	while True:
		bytesToRead = ser.readline()
		if (bytesToRead!=""):
			started=1
			ii+=1
			if (ii<10):
				countStr=str("0%s" %ii)
			else: 
				countStr=ii
			datafiltered=""
			for x in range (0,20):
				sys.stdout.write("\rLoading Grid %s\tL%s" %(ii,(x+1)))
				sys.stdout.flush()
				line=""
				for y in range (0,20):
					
					if (x==0 and y==0):
						dataline=str(bytesToRead)
					else:
						dataline=str(ser.readline())
					metadataline=str(ser.readline())
					if ("4095" in dataline):  
						value="0"
					elif ("4094" in dataline):
						value="0.01"
					else:
						if ('11' in metadataline):
							value=int(dataline.strip())*10
						else:
							value=int(dataline.strip())
						value=value/float(255)
						value= str("%.2f" % round(value, 2))#convert to resistance values
					if (x%2!=0):#Zigzag mode
						if (y >0):
							value+=","
						line=str(value)+line
					else:
						if (y < 19):
							value+=","
						line+=str(value)#/(float(255)))
						
				if (x<20):
					datafiltered+=line+"\n"
			fname='%ssa1_%s.csv' %(pathCsv,countStr)
			if (os.path.isfile(fname) and overwrite=="ask"):
				while True:
					overwrite=raw_input('\nGrid already downloaded. Overwite?\n YA:Yes for all\n Y:Yes\n N: Not this one\n NA: Only add new grids\n').lower()
					if (overwrite=="y" or overwrite=="ya"):
						fileout= open(fname, 'wb')
						fileout.write(datafiltered)
						fileout.close()
						sys.stdout.write(" Existed (Overwritten)")
						break
					elif (overwrite=="n" or overwrite=="na"):
						sys.stdout.write(" Existed (Not Overwritten)")
						break
				if (overwrite=="y" or overwrite=="n"):
					overwrite="ask"
			elif (os.path.isfile(fname) and overwrite=="ya"):
				fileout= open(fname, 'wb')
				fileout.write(datafiltered)
				fileout.close()
				sys.stdout.write(" Existed (Overwritten)")
			elif (os.path.isfile(fname) and overwrite=="na"):
				sys.stdout.write(" Existed (Not Overwritten)")
			else:
				fileout= open(fname, 'wb')
				fileout.write(datafiltered)
				fileout.close()
			
			
			
			sys.stdout.write("\n")
			sys.stdout.flush()
		if (started==1 and bytesToRead==''):
			print('Finished')
			break


	
def testFilterdata():
	inputName=raw_input('Enter file name:')
	file= open('C:\\sync\\wrk\\Appia0216\\RM85_data\\raw\\'+inputName+'.txt', 'r')
	
	data_split= file.readlines()
	numberOfGrids= int(len(data_split)/float(800))
	print("%s grids loaded" % numberOfGrids)
	out=""
	datafiltered=""
	for ii in range (0, numberOfGrids):
		out+="----Grid %s----\n" %(ii+1)
		datafiltered=""
		for x in range (0,20):
			for y in range (0,20):
				if (x%2!=0):#Zigzag mode
					y_effect= 19-y
				else:
					y_effect=y
				i= 2*(ii*400+(x*20)+(y_effect))
				if "4095" in data_split[i]:  
					datafiltered+="0"
				else:
					datafiltered+='%.2f' %(int(data_split[i].strip())/(float(10)))
					
				if (y < 19):
					datafiltered+=","
			datafiltered+="\n"
		datafiltered+="\n"
		out=datafiltered
		fileout= open('C:\\sync\\wrk\\Appia0216\\RM85_data\\csv\\sa1_10%s_.csv' %ii, 'wb')
		fileout.write(out)
		fileout.close()
	file.close()
	
def Excel2CSV():
	import os
	import xlrd
	import csv
	dir='C:\\sync\\wrk\\Appia0216\\RM85_data\\matrix_excel\\'
	listXls=os.listdir(dir)
	for i in range (0, len(listXls)):
		ExcelFile=dir+listXls[i]
		outputname=((listXls[i].split('grid'))[1].split('.'))[0]
		print outputname
		CSVFile='C:\\sync\\wrk\\Appia0216\\RM85_data\\csv\\sa1_10%s.csv' %outputname
		workbook = xlrd.open_workbook(ExcelFile)
		worksheet = workbook.sheet_by_index(0)
		csvfile = open(CSVFile, 'wb')
		wr = csv.writer(csvfile, quoting=csv.QUOTE_NONE)

		for rownum in xrange(worksheet.nrows):
			row= worksheet.row_values(rownum)
			newrow=lst = [None] * len(row)
			for ii in range (0,len(row)):
				row[ii]=float(row[ii])
				if (row[ii]=='4095.0'):
					newrow[ii]=0.0
				else:
					newrow[ii]=float(row[ii])/10
			wr.writerow(newrow)

		csvfile.close()

def histeq(im,nbr_bins=256):
	from PIL import Image
	from numpy import *
	#get image histogram
	imhist,bins = histogram(im.flatten(),nbr_bins,normed=True)
	cdf = imhist.cumsum() #cumulative distribution function
	cdf = 255 * cdf / cdf[-1] #normalize

	#use linear interpolation of cdf to find new pixel values
	im2 = interp(im.flatten(),bins[:-1],cdf)

	return im2.reshape(im.shape), cdf
   
def despike(array):
	cols=array.shape[0]
	rows=array.shape[1]
	despikeArray=array
	count=0
	for x in range (1, cols-2):
		for y in range (1, rows-2):
			value=array[x][y]
			neighbours=[array[x-1][y-1],array[x][y-1],array[x+1][y-1],array[x-1][y],array[x+1][y],array[x-1][y+1],array[x][y+1],array[x+1][y+1]]
			if (neighbours.count(0.0)<3):
				mean=sum(neighbours)/float(len(neighbours))			
				for n in (0,7):
					if (neighbours[n]/mean>float(1.5) or neighbours[n]/mean<float(0.66) ):
						oddNeighbours=1
					else:
						oddNeighbours=0
				if (oddNeighbours==0):
					if (value/mean>3 or value/mean<0.33 ):
						despikeArray[x][y]=int(mean)
						count+=1
					else:
						despikeArray[x][y]=value
			else:
				despikeArray[x][y]=value
	print("%s values despiked." %count)
	return despikeArray
	
def dataInputRes2d():
	import numpy as np
	#file=raw_input('Enter file')
	file='C:\\Data\\Tapino\\S.Andrea_TAP23\\cross-section\\Test2\\tap23_cs1_dd.csv'
	gridname=raw_input('Enter grid name')
	mat=np.genfromtxt(file,delimiter=',');
	
	arrayNumber=3 #3 for dipole dipole
	length=19
	spacing=0.5
	levels=10
	numberOfPoints=((length-2)*levels)-(levels*(levels-1)/2)
	xpos=1	#0for 1st electrode, 1 for midpoint
	
	
	line1=int(length-3)
	
	count=-1
	fileString=('%s\n%s\n%s\n%s\n%s\n0\n' %(gridname,spacing,arrayNumber,numberOfPoints,xpos))
	print(mat.shape)
	for i in range (0,(levels)):
		for ii in range (0, (line1)):
			if mat[i][ii]!=0:
				count+=1
				value=mat[i][ii]*(i+1)*(i+2)*(i+3)
				y=i+1
				x=i*0.5*spacing+ii*spacing
				fileString+=('%s,   %s,   %s,   %s\n' %(x,spacing,y,value))
	
	print(fileString)
	
	text_file = open("C:\\Data\\Tapino\\S.Andrea_TAP23\\cross-section\\22082015\\"+gridname+"_output.dat", "w")
	text_file.write(fileString)
	text_file.close()
	
	
def processDipoleDipole(tempfile):
	import numpy as np
	import matplotlib
	import matplotlib.mlab as mlab
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	while True:
		dest=raw_input('Enter your destination:')
		if os.path.exists(dest):
			listFiles=os.listdir(dest)
			for i in range (0,len(listFiles)):
				if '.csv' in listFiles[i]:
					if 'files' in locals():
						files.extend([listFiles[i]])
					else:
						files=list([listFiles[i]])
			if 'files' not in locals():
				print('No valid files found')
				break;
			else:
				while True:
					print('Select your file:')
					for i in range (0,len(files)):
						print('%s   %s' %(i, files[i]))
					fileId=input()
					if fileId in range(0,len(files)):
						break;
				tempfile=dest+'\\'+files[fileId]
				gridname=files[fileId]
				
				#mergeDipoleDipole(files, dest)
				
				
				print(gridname+' is being processed')
				
				
				
				mat=np.genfromtxt(tempfile,delimiter=',');
				length=25
				spacing=1
				
				line1=int(((length+spacing)/spacing)-3)
				levels=6
				count=-1
				max=(line1)+(line1-1)+(line1-2)+(line1-3)+(line1-4)+(line1-5)
				value=np.zeros(max)
				y=np.zeros(max)
				x=np.zeros(max)
				print(mat.shape)
				for i in range (0,(levels)):
					for ii in range (0, (line1)):
						if mat[i][ii]!=0:
							count+=1
							value[count]=mat[i][ii]*(i+1)*(i+2)*(i+3)
							y[count]=-i
							x[count]=i*0.5+ii
				print('Length = %s m'%length)
				print('Depth levels = %s'%levels)
				
				plt.figure()
				plt.tricontourf(x,y,value, 100)
				plt.title('Pseudo Section - dipole dipole \n ' + gridname)
				
				plt.show()
				
			break;
		else:
			print('Invalid input')
			break;

def mergeDipoleDipole(arraylist,dir):
	import numpy as np
	overlap=input('Overlap?')
	depthlevels=input('Depth levels?')
	if len(arraylist)<=1:
		print('You need at least two grids')
	else:
		tempfile=('%s\\%s'%(dir,arraylist[0]))
		print(tempfile)	
		mat=np.genfromtxt(tempfile,delimiter=',');
		PS=mat[0:depthlevels][0:-(overlap-1)]
		print(PS)
		overlap=mat[0:depthlevels][-(overlap-1):]
		print(overlap)
		for i in range(1,len(arraylist)):
			tempfile=('%s\\%s'%(dir,arraylist[i]))
			mat=np.genfromtxt(tempfile,delimiter=',');
			overlap2=mat[0:depthlevels][:(overlap-1)]
			averaged=(overlap+overlap2)/2
			PS.extend(averaged)
			if i==len(arraylist):
				PS.extend(mat[0:depthlevels][(overlap-1):])
			else:
				PS.extend(mat[0:depthlevels][(overlap-1):-(overlap-1)])
	print(PS.shape)
	print(PS)	

def findNewGrids():
	listTif=os.listdir(pathTif)
	listCsv=os.listdir(pathCsv)
	for i in range (0,len(listTif)):
		str=listTif[i].split(".")
		if (str[1]=="tif"):
			a=str[0].split("_")
			site=a[0]
			gridId=a[1]
			if 'tifs' in locals():
				tifs.extend([gridId])
			else:
				tifs=list([gridId])
	if 'tifs' not in locals():
			tifs=""

	for i in range (0,len(listCsv)):
		str=listCsv[i].split(".")
		if (str[1]=="csv"):
			a=str[0].split("_")
			site=a[0]
			gridId=a[1]
			try:
				rotate=a[2]
			except:
				rotate=""
			if 'csvs' in locals():
				csvs.extend([gridId])
			else:
				csvs=list([gridId])
		
	for i in range (0,len(csvs)):
		found=1
		for ii in range (0,len(tifs)):
			if csvs[i]==tifs[ii]:
				found=0
	
		if (found!=0):	
			if 'newGrids' in locals():
				newGrids.extend([site+'_'+csvs[i]])
			else:
				newGrids=list([site+'_'+csvs[i]])

	if 'newGrids' in locals():
		print ('%d unprocessed grids found:\n %s' %(len(newGrids),newGrids))
		csv2tif(newGrids)
	else:
		print ('All grids already converted to GeoTiffs')
		print('\n')
		
	del (csvs ,tifs ,listCsv, listTif)

def csv2tif(newGrids):
	import numpy as np
	batchChoice=raw_input("Do you need to rotate(enter to skip, y for yes)? ")
	
	for tempfile in newGrids:
		mat=np.genfromtxt(pathCsv+tempfile+'.csv',delimiter=',');
		if "yes" in batchChoice:
			rotate=raw_input(tempfile + ": Rotate? (rot: 90 180 flip: ud lr -- enter to skip\n\
 _____   _____   _____   _____   _____   _____   _____   _____\n\
||    | |-    | |    || |    -| |     | |     | |     | |     |\n\
|     | |     | |     | |     | |     | |     | |     | |     |\n\
|_____| |_____| |_____| |_____| |____-| |____|| |-____| ||____|\n\
   1       2       3       4       5       6       7       8\n")
			if '90' in rotate:
				mat=np.rot90(mat,1)
			elif '180' in rotate:
				mat=np.rot90(mat,2)
			elif '270' in rotate:
				mat=np.rot90(mat,3)
			if 'ud' in rotate:
				mat=np.flipud(mat)
			elif 'lr' in rotate:
				mat=np.fliplr(mat)
		mat=np.flipud(mat)
		x=mat.shape[1]
		y=mat.shape[0]
		step=1
		mat=(100*mat)
		filename=tempfile
		exportGeoTiff(filename, mat);

def getGeometry():
	import numpy as np
	index=np.genfromtxt(pathTif+'geometry.txt',dtype='string');
	return index;
	
def convertCoordinates(inEPSG, outEPSG, x, y):
	from pyproj import Proj, transform
	inProj = Proj(init='epsg:'+inEPSG)#whatever system in metric units
	outProj = Proj(init='epsg:'+outEPSG)
	x2,y2 = transform(inProj,outProj,x,y)
	return x2,y2
	
def getGeoRef():
	import math
	try:
		text_file = open(pathTif+'georef.txt', "r")
		for line in text_file: 
			if "EPSG" in line:
				lineSplit=line.split("=")
				epsg=lineSplit[1].strip('\n')
			elif "GRIDPOINT1" in line:
				lineSplit=line.split("=")
				lineSplit=lineSplit[1].split(",")
				gridPoint1=(lineSplit[0],lineSplit[1].strip('\n'))
			elif "GRIDPOINT2" in line:
				lineSplit=line.split("=")
				lineSplit=lineSplit[1].split(",")
				gridPoint2=(lineSplit[0],lineSplit[1].strip('\n'))
			elif "COORDPOINT1" in line:
				lineSplit=line.split("=")
				lineSplit=lineSplit[1].split(",")
				coordPoint1=(lineSplit[0],lineSplit[1].strip('\n'))
			elif "COORDPOINT2" in line:
				lineSplit=line.split("=")
				lineSplit=lineSplit[1].split(",")
				coordPoint2=(lineSplit[0],lineSplit[1].strip('\n'))
		print('Georeferencing from file...')
		georef=[epsg,gridPoint1,coordPoint1,gridPoint2,coordPoint2]
		print(georef)
	except:
		print('No georeferencing file found')
		georef=""
	if (georef!=""):
		#define angle of grids
		dXGrid=float(gridPoint2[0])-float(gridPoint1[0])
		dYGrid=float(gridPoint2[1])-float(gridPoint1[1])
		if (dYGrid==0):
			anglegridAxix=0
		else:
			anglegridAxix=math.atan(dXGrid/dYGrid)
		print (anglegridAxix)
		
		x1,y1 = coordPoint1[0], coordPoint1[1]#convertCoordinates('3004', '32633', coordPoint1[0], coordPoint1[1])
		x2,y2 = coordPoint2[0], coordPoint2[1]#convertCoordinates('3004', '32633', coordPoint2[0], coordPoint2[1])
		
		dX=float(x2)-float(x1)
		dY=float(y2)-float(y1)
		angle=math.atan(dX/dY)
		print(angle)
		dAngle=angle-anglegridAxix
		#Need to define the size of the grid in the coordinate system
		
		
		
		
		#find origin of the grid
		dYtoOrigin= math.sin(dAngle)*20*float(gridPoint1[1])
		print(dYtoOrigin, float(y1))
		dXtoOrigin= math.cos(dAngle)*20*float(gridPoint1[0])
		print(dXtoOrigin, float(x1))
		print (dXtoOrigin, dYtoOrigin)
		origin=(float(x1)+dXtoOrigin, float(y1)+dYtoOrigin)
		pi = math.pi
		angleDegrees = (180 * dAngle / pi)-90
		print(angleDegrees)
		xorig, yorig=convertCoordinates('3004', '3004', origin[0], origin[1])
		georef= [angleDegrees, xorig, yorig, 3004]
		print (georef)
	return georef
	
def csv2matrix():
	import numpy as np
	import os
	dir='C:\\sync\\wrk\\Appia0216\\GeoRes\\frank\\convert\\'
	list=os.listdir(dir)
	print('Converting %s files' % len(list))
	for a in range(0,len(list)):
		filename=dir+list[a]
		mat=np.genfromtxt(filename,delimiter=',')
		for i in range (0,mat.shape[0]):
			for ii in range (0,mat.shape[1]):
				if (mat[ii][i]==409.5):
					mat[ii][i]=0
				elif (mat[ii][i]!=0.0):
					mat[ii][i]=round(mat[ii][i]/2.55, 2)
		np.savetxt(filename, mat, delimiter=",", fmt='%.2f')

def hillshade(array, azimuth, angle_altitude, z): 
    from numpy import gradient
    from numpy import pi
    from numpy import arctan
    from numpy import arctan2
    from numpy import sin
    from numpy import cos
    from numpy import sqrt
    from numpy import zeros
    from numpy import uint8
	
    nodata=array!=0
    array=z*array
    x, y = gradient(array)
    slope = pi/2. - arctan(sqrt(x*x + y*y))
    aspect = arctan2(-x, y)
    azimuthrad = azimuth*pi / 180.
    altituderad = angle_altitude*pi / 180.
     
 
    shaded = sin(altituderad) * sin(slope)\
     + cos(altituderad) * cos(slope)\
     * cos(azimuthrad - aspect)
    shaded[shaded==0]=1
    shaded=255*(shaded + 1)/2
    
    
    return shaded*nodata
   				
def exportGeoTiff( filename, mat):
	import numpy as np
	#get coordinates
	mat=np.int16(mat)
	mat=np.flipud(mat)
	georef=getGeoRef()
	if (georef!=''):
		xCoord=georef[1]
		yCoord=georef[2]
		rotY=georef[0]+90
		rotX=georef[0]+90
		epsg=georef[3]
		print('Autogeoref...')
	else:
		xCoord=0
		yCoord=0
		rotX=0
		rotY=0
		epsg=32633
		print('No georeferencing')
	xSize=mat.shape[1]
	ySize=mat.shape[0]
	step=1
	output_file=pathTif+filename+'.tif'
	# Create gtif (uncomment when using gdal)
	# 
	DataType= gdal.GDT_UInt16
	driver = gdal.GetDriverByName("GTiff")
	dst_ds = driver.Create(output_file, xSize, ySize, 1, DataType)

	#top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
	dst_ds.SetGeoTransform( [ xCoord , step , rotX, yCoord ,  rotY, step ] )

	#set the reference info 
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(epsg)
	dst_ds.SetProjection( srs.ExportToWkt() )
	


	# write the band
	dst_ds.GetRasterBand(1).WriteArray(mat)
	band = dst_ds.GetRasterBand(1)
	band.SetNoDataValue(0)
	band.FlushCache()

	
def findHoles():
	import numpy as np
	import matplotlib
	import matplotlib.mlab as mlab
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from PIL import Image
	import scipy as sp
	import scipy.ndimage.morphology
	image = Image.open(pathTif+'ciocca_mosaic.tif')
	array=np.array(image)
	nulls=(array==0)
	objects, num_objects = sp.ndimage.label(nulls)
	plt.subplot(221), plt.imshow(image)
	plt.title('degraded image')
	plt.subplot(222), plt.imshow(nulls, 'gray')
	plt.title('mask image')
	plt.subplot(223), plt.imshow(objects, 'gray')
	plt.title('holes')
	
	plt.tight_layout()
	plt.show()

	
def convertData():
	import numpy as np
	listTif=os.listdir(pathTif)
	for i in range (0,len(listTif)):
		print (pathTif+listTif[i])
		try:
			im = Image.open(pathTif+listTif[i])
		except:
			print ('Warning: Grid '+index[i][ii]+'  not found.')
		
		A=np.array(im)
		x=A.shape[1]
		y=A.shape[0]
		step=1
		mat=(10*A)
		filename=listTif[i].split(".")
		filename=filename[0].split('grid')
		filename="sa3_"+filename[1]
		print (x +y)
		print (pathTif+filename)
		exportGeoTiff( filename, mat);

def equalizeHorz(A,B):
	import numpy as np
	from PIL import Image
	import numpy as np
	if (A.shape[0]!=B.shape[0]):
		print('Warning sizes ara different!')
	ALastLine=A[:,A.shape[1]-1]
	BLastLine=B[:,B.shape[1]-1]
	if ((BLastLine==0.0).sum()<B.shape[1] and (BLastLine==0.0).sum()<A.shape[1]):
		diffArray=np.empty(A.shape[0])
		for i in range (0, A.shape[0]):
			if (ALastLine[i]!=0 and BLastLine[i]!=0):
				diffArray[i]=float(ALastLine[i])-float(BLastLine[i])
		preDiffMean=np.median(diffArray)
		count=0
		ii=0
		diffArrayExtremesOut=0
		print(len(diffArray))
		for i in range (0, len(diffArray)):
			if (diffArray[i]/preDiffMean<float(0.2) or diffArray[i]/preDiffMean>float(5)):
				count+=1
			else:
				diffArrayExtremesOut+=diffArray[i]
				ii+=1
		diffMean=diffArrayExtremesOut/float(ii)
		print(diffMean)
		if (ii>10):
			BEqualized=np.empty(B.shape, dtype='float64')
			BEqualized=B+diffMean#/float(2)
		else:
			print('not enough correspondances')
			BEqualized=B
		print('%s values used. %s ignored. Average: %s' %(ii, count, diffMean))
	else:
		print ('empty grid neighbour')
		BEqualized=B
	return BEqualized;
	
def makeMosaic():
	import numpy as np
	from PIL import Image
	import numpy as np
	gridSize=(20,20)
	index=getGeometry()
	zero=np.zeros(gridSize)
	x=index.shape[0]
	y=index.shape[1]
	for i in range (0,x):
		for ii in range (0,y):
			
			if (index[i][ii]=='0'):
				A=zero
			else:
				try:
					im = Image.open(pathTif+sitename+'_'+index[i][ii]+'.tif')			
				except:
					print ('Warning: Grid '+index[i][ii]+'  not found.')
					im=zero
				A=np.array(im)
				A=np.flipud(im)
			if (ii==0):
				line=A
			else:
				#print('--------Grid: %s - %s -------' %(i, ii))
				grid=A
				#grid=equalizeHorz(line,A)
				line=np.concatenate((line,grid),axis=1)
				
		if (i==0):
			mos=line
		else:
			mos=np.concatenate((mos,line),axis=0)
	filename=sitename+'_mosaic2'
	step=1
	xtif=y*gridSize[0]
	ytif=x*gridSize[1]
	mos=despike(mos)
	print('%s x %s' %(x,y))
	exportGeoTiff(filename, mos)
	hs_array = hillshade(mos,45, 315, 0.0025)
	exportGeoTiff(filename+"_hs", hs_array)
	mos,cdf = histeq(mos)
	makePreview(mos,filename)
	

def makePreview(array, filename):
	import numpy as np
	import matplotlib
	import matplotlib.mlab as mlab
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from PIL import Image
	import cv2
	import scipy.ndimage
	x=array.shape[0]
	y=array.shape[1]
	min= 99999
	max=0

	for i in range (0,x):
		for ii in range (0,y):
			if (array[i][ii]!=0 and array[i][ii]<min):
				min=array[i][ii]
			if (array[i][ii]!=0 and array[i][ii]>max):
				max=array[i][ii]
	print('%s %s %s %s' %(min,max, x, y))

	

	mask=np.empty([x, y], dtype='float64')
	for i in range (0,x):
		for ii in range (0,y):
			if (array[i][ii]!=0):
				mask[i][ii]=255
			if (array[i][ii]==0):
				mask[i][ii]=0
	
	
	if (min!=max and max!=0):
		array=(1-(array-min)/(max-min))
		array = cv2.resize(array, (0,0), fx=5, fy=5) 
		mask = cv2.resize(mask, (0,0), fx=5, fy=5) 
		array= scipy.ndimage.median_filter(array, 4)	
		values = Image.fromarray(np.uint8(cm.jet_r(array)*255)).convert('RGB')
		
		
		
			
		hs_array = Image.fromarray(np.uint8(hillshade(array*255,45, 315, 0.5))).convert('RGB')
		new_img = Image.blend(values, hs_array, 0.5).convert('RGBA')
		mask = Image.fromarray(np.uint8(mask)).convert('L')
		new_img.putalpha(mask)
		new_img.save(pathTif+filename+'prev.png')
		img = Image.open(pathTif+filename+'prev.png')
		img.show() 
	else:
		print('error in reading image')
		
def main():

	while(True):
		menu=raw_input(" l     Process All\n m     Create mosaic \n d     Download data RM85\n f     Find Holes regions \n e     Extract WMS\n c     Create Cross section \n q     Quit\n\n")
		if menu == 'q':
			exit();

		elif menu == "l":
			findNewGrids()
			makeMosaic()
			continue;
	
		elif menu == "m":
			makeMosaic()
			continue;
		
		elif menu == "c":
			csv2matrix()
			continue;
		
		elif menu == "e":
			#Excel2CSV()
			import numpy as np
			from PIL import Image
			import numpy as np
			A=np.array(Image.open(pathTif+sitename+'_'+'1054'+'.tif'))
			B=np.array(Image.open(pathTif+sitename+'_'+'1045'+'.tif'))
			equalizeHorz(A,B)
			continue;
			
		elif menu == "d":
			downloadRM85()
			#testFilterdata()
			continue;
		
		elif menu == "f":
			im = Image.open('C:\\Data\\LERC\\Isernia\\GeoRes\\olive\\tif\\a145_mosaic.tif')
			A=np.array(im)
			y=A.shape[0]
			x=A.shape[1]
			filename="a405_mosaic2"
			exportGeoTiff( x, 1,y, filename, A)
			hs_array = hillshade(A,45, 180, 0.0025)
			exportGeoTiff( x, 1,y, filename+"_hs", hs_array)
			makePreview(A,filename)
			makePreview(hs_array,filename+"_hs")
			
			

			continue;

main()


