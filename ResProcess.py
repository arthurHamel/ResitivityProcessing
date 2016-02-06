

import os
import sys
import gdal, osr, ogr

sitename='sa24a'
path= 'C:\\sync\\wrk\\Appia0216\\Hekje_nieuw\\'
pathTif= path + '\\tif\\'
pathCsv= path + '\\raw\\'

def downloadRM85():
	import os.path
	import serial
	ser = serial.Serial('COM3', 9600, timeout=1) #Tried with and without the last 3 parameters, and also at 1Mbps, same happens.
	ser.flushInput()
	ser.flushOutput()
	overwrite="ask"
	print("Waiting for data...")
	started = 0
	ii=0
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
			fname='%ssa24a_1%s.csv' %(pathCsv,countStr)
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
	file= open('C:\\sync\\wrk\\Appia0216\\Hekje_niew\\test2.txt', 'r')
	
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
					datafiltered+=str(int(data_split[i].strip()))#/(float(255)))
					
				if (y < 19):
					datafiltered+=","
			datafiltered+="\n"
		datafiltered+="\n"
		out=datafiltered
		fileout= open('C:\\sync\\wrk\\Appia0216\\Hekje_niew\\raw\\sa24a_10%s.csv' %ii, 'wb')
		fileout.write(out)
		fileout.close()
	file.close()
	
	
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
	batchChoice=raw_input("NewGrids to Tiffs? ")
	for tempfile in newGrids:
		mat=np.genfromtxt(pathCsv+tempfile+'.csv',delimiter=',');
 		rotate=raw_input(tempfile + ": Rotate? (rot: 90 180 flip: ud lr -- enter to skip")
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
		exportGeoTiff( x, step,y, filename, mat,'no');

def getGeometry():
	import numpy as np
	index=np.genfromtxt(pathTif+'geometry.txt',dtype='string');
	return index;
def csv2matrix(filename):
	import numpy as np
	mat=np.genfromtxt(pathCsv+filename+'*.csv',delimiter=',');

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
   				
def exportGeoTiff( x, step,y, filename, mat,georef):
	import numpy as np
	#get coordinates
	mat=np.int16(mat)
	if (georef=='yes'):
		shapefile = path+"grid\\fishnet.shp"
		driver = ogr.GetDriverByName("ESRI Shapefile")
		dataSource = driver.Open(shapefile, 0)
		layer = dataSource.GetLayer()
		sref=layer.GetSpatialRef()
		feature_elmt = layer.GetNextFeature()
		field_vals = []
		while feature_elmt:
			field_vals.append(feature_elmt.GetFieldAsString('Grid_ID'))
			feature_elmt = layer.GetNextFeature()
			while True:
				gridnum=(filename.split("_"))[1];
				gridID=gridnum;
				if gridID in field_vals:
					#format the request
					querry="Grid_ID = '"+gridID+"'";
					layer.SetAttributeFilter(querry)
					for feature in layer:
						xCoord=float(feature.GetField("X"))
						yCoord=float(feature.GetField("Y"))

				else: 
					print('No grid with this reference.')
				break;
	else:
		xCoord=0
		yCoord=0
	xSize=int(x/step)
	ySize=int(y/step)
	output_file=pathTif+filename+'.tif'
	# Create gtif (uncomment when using gdal)
	# 
	DataType= gdal.GDT_UInt16
	driver = gdal.GetDriverByName("GTiff")
	dst_ds = driver.Create(output_file, xSize, ySize, 1, DataType)

	#top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
	dst_ds.SetGeoTransform( [ xCoord , step , 0, yCoord , 0, step ] )

	#set the reference info 
	srs = osr.SpatialReference()
	srs.ImportFromEPSG(32633)
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
		exportGeoTiff( x, step,y, filename, mat,'no');
		
	
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
				line=np.concatenate((line,A),axis=1)
				
		if (i==0):
			mos=line
		else:
			mos=np.concatenate((mos,line),axis=0)
	filename=sitename+'_mosaic2'
	step=1
	xtif=y*gridSize[0]
	ytif=x*gridSize[1]
	print('%s x %s' %(x,y))
	exportGeoTiff( xtif, step,ytif, filename, mos,'no')
	hs_array = hillshade(mos,45, 315, 0.0025)
	exportGeoTiff( xtif, step,ytif, filename+"_hs", hs_array,'no')
	makePreview(mos,filename)
	
def extractWMS():
	from PIL import Image
	from owslib.wms import WebMapService
	urlList=['http://servizi.geo.regione.molise.it/arcgis/services/Uso_suolo/MapServer/WMSServer', 
		'http://servizi.geo.regione.molise.it/arcgis/services/Viabilita_Teleatlas/MapServer/WMSServer',
		'http://servizi.geo.regione.molise.it/arcgis/services/Ortofoto_Molise_2007/MapServer/WMSServer',
		'http://wms.pcn.minambiente.it/ogc?map=/ms_ogc/WMS_v1.3/raster/ortofoto_colore_12.map'
	]
	
	dest='\\VUW\Personal$\Homes\H\hamelac1\Downloads'
	filename='WMS_output'
	print('Available WMS servers (enter number or empty for custom)')
	for i in range(0,len(urlList)):
		print ('%s   %s' %(i+1,urlList[i]))
	choice=input()
	if (not choice):
		wmsUrl=raw_input('Enter custom WMS address:')
	else:
		wmsUrl=urlList[choice-1]
	print('Getting Info...')
	wms = WebMapService(wmsUrl, version='1.1.1')
	print('Service Name : '+wms.identification.title)
	layerList=list(wms.contents)
	print('Available layers')
	for i in range(0,len(layerList)):
		print ('%s   %s' %(i+1,layerList[i]))
	choiceLayer=input()
	layer=layerList[choiceLayer-1]
	boundingBox=(14.50, 41.65, 14.55, 41.68) #xmin,ymin,xmax, ymax   
	#a205: 14.264, 41.582 , 14.273, 41.591
	#poi4026:  14.346355, 41.578670, 14.354092, 41.585850
	#a405: 14.166, 41.619, 14.70, 41.623
	#laromana: 14.179, 41.618, 14.192, 41.624 
	#Tap03: 14.731, 41.572, 14.738, 41.575
	#CSA:  14.522, 41.656 14.532, 41.666
	#Isernia: 14.071, 41.509, 14.361, 41.679
	
	sizeImg=(2000,2000)
	print ('Getting '+layer)
	xextend=boundingBox[2]-boundingBox[0]
	yextend=boundingBox[3]-boundingBox[1]
	
	if (xextend>0.01 and yextend>0.01):
		rows=yextend/0.01;
		cols=xextend/0.01;
		
		blank_im = Image.new('RGB', (cols*sizeImg[0],rows*sizeImg[1]))
		
		for i in range (0, rows):
			windowBBox=(0, boundingBox[1] + (i*0.01), 0, boundingBox[3] + (i*0.01))
			for ii in range (0, cols):
				windowBBox[0]=boundingBox[0] + (i*0.01)
				windowBBox[2]=boundingBox[2] + (i*0.01)
				img = wms.getmap(layers=[layer],
				srs = 'EPSG:4326',
				bbox = windowBBox,
				size = sizeImg,
				format = 'image/png',
				transparent = True
				)
				blank_im.paste(img, (i*sizeImg[0], i*sizeImg[1]))
				blank_im.show()
		img=blank_im
		
	else: 	
		img = wms.getmap(layers=[layer],
        srs = 'EPSG:4326',
        bbox = boundingBox,
        size = sizeImg,
        format = 'image/png',
        transparent = True
		)
	
	#project
	
	xCoord=boundingBox[0]
	yCoord=boundingBox[1]
	ppx=(boundingBox[2]-boundingBox[0])/sizeImg[0]
	ppy=(boundingBox[3]-boundingBox[1])/sizeImg[1]
	
	
	
	
	
	out = open(dest+filename+'.tif', 'wb')
	out.write(img.read())
	out.close() 
	
	image = Image.open(dest+filename+'.tif')
	imarray=np.array(image)
	imarray=np.flipud(imarray)
	
	while True:
		try:
			output_file=dest+filename+'.tif'
			# Create gtif (uncomment when using gdal)
			DataType= gdal.GDT_Byte
			driver = gdal.GetDriverByName("GTiff")
			dst_ds = driver.Create(output_file, sizeImg[0], sizeImg[1], 3, DataType)
			#top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
			dst_ds.SetGeoTransform( [ xCoord , ppx , 0, yCoord , 0, ppy ] )
			#set the reference info 
			srs = osr.SpatialReference()
			srs.ImportFromEPSG(4326)
			dst_ds.SetProjection( srs.ExportToWkt() )	
			# write the band
			dst_ds.GetRasterBand(1).WriteArray(imarray[:,:,0])
			dst_ds.GetRasterBand(2).WriteArray(imarray[:,:,1])
			dst_ds.GetRasterBand(3).WriteArray(imarray[:,:,2])
			dst_ds.FlushCache()
			break;
		except:
			print ('error wrinting the file\n Delete/close previous version and press enter')
			raw_input()
	
	image.show()
	print('Done')

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
			dataInputRes2d()
			continue;
		
		elif menu == "e":
			convertData()
			#extractWMS()
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
			exportGeoTiff( x, 1,y, filename, A,'no')
			hs_array = hillshade(A,45, 180, 0.0025)
			exportGeoTiff( x, 1,y, filename+"_hs", hs_array,'no')
			makePreview(A,filename)
			makePreview(hs_array,filename+"_hs")
			
			

			continue;

main()


