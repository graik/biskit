def dbPre (str ,type='int',default=0,mandatory = False):
	if   'X' in str:
		return default

	if type == 'int':
		return int(str)
	elif type== 'float':
		return float(str)
	elif type == 'int_range':
		aux = str.split(',')
		return [int(aux[0]),int(aux[1])]
	else:
		print "[ERROR DBpre] Type not known ( "+type+" )."
		return default

def overlapCalc ( wl,acc_spectra,don_spectra, e_cof ):
	from numpy import sum
	
	area = sum(don_spectra)
	print  "Area: "+str(area)	
	don_spectra = don_spectra/area
	newintegral =sum(don_spectra)
	print "Now area is: "+str(newintegral)
	e_spectra =acc_spectra*e_cof
	return sum( wl**4*don_spectra*e_spectra)
	
def sphericalVectors( latdeg =0.087266462599716474  , longdeg = 0.087266462599716474):
	from emath import rotation
	from math import pi
	from numpy import array,matrix
	
	vectors = []
	others = []
	alpha = -pi/2.
	while alpha <pi/2.:
		beta = -pi/2.
		while beta < pi/2.:
			v = array([1.,0,0])*matrix(rotation(beta,alpha,0.))
			vectors.append((v[0,0], v[0,1],v[0,2]))
			others.append((-v[0,0], v[0,1],v[0,2]))
			beta+=latdeg
		alpha +=longdeg
		others.reverse()
		cycle = len(others)*2
		for v in others:
			vectors.append(v)
		others = []
	return vectors, cycle

def create3DFRETEfficiencySphere(  f = None,distance=(1,0,0),donor = (1,0,0), acceptors = []):
	from emath import norm
	
	r = norm(distance)
	results = []
	
	for a in acceptors:
		f.calcK2( donor ,a, distance)
		results.append(( a[0],a[1],a[2],f.energyTransferEfficiency(r)))
	
	return results

def plot3DFRETEfficiencySphere( data = [] ,cycle = 0,filename= "",saveData = False ,more =("",),script = ("set hidden\n","set hidden3d\n","set pm3d\n","set pm3d depthorder\n","set ticslevel 0\n","splot \"spherescriptdata\" w l title \"Efficiency\"\n","pause 5\n",\
																							"set terminal png\n","set output \"myout\"\n","replot" )):
	import os
	
	# Gen. data file
	f = open("spherescriptdata",'w')
	c = 0
	for v in data:
		c +=1
		if c%cycle !=0:
			f.write(str(v[0])+" "+str(v[1])+" "+str(v[2])+" "+str(v[3])+"\n")
		else:
			f.write("\n")
	f.close()
	if not saveData:
		os.system("rm spherescriptdata")
	
	f = open("plotspherescript",'w')
	f.writelines(more + script)
	f.close()
	
	os.spawnlp(os.P_WAIT, 'gnuplot', 'gnuplot', 'plotspherescript')
	
	if filename != "":
		os.system("mv myout "+filename+".png")
	if not saveData:
		os.system("rm plotspherescript")
		os.system("rm myout")



from FRET import FRET


overlap = 959732990869504
qyD = 0.62
 
f = FRET( overlap, qyD,44e3)

acceptors , cycle = sphericalVectors( )

#~ print acceptors

results  = create3DFRETEfficiencySphere(  f, [0.,60.,0.] ,[0.,0.,1.], acceptors  )

#~ print results


plot3DFRETEfficiencySphere(results,cycle,"lol",True,("set view map\n",))