#!/usr/bin/python

#We thank you in advance for sending your feedback and/or 
#suggestions to:
#             horacio.v.g@gmail.com
#
#@@@@@@@@@@@@@@@@@@@@@2011@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#What we observe is not nature itself, but nature exposed to
#our method of questioning. 
#                                  W. Heisenberg
############################################################
#===========================================================
#
#========================Llajtamasis========================
# IMPORT MODULES
import os
import gtk
import gtk.glade
import pygtk # temp

#====================POINT MASS MODEL=======================
# IMPORT MODULES
from scipy import *
import mypypm7
from pylab import *
import numpy as np
# Define the path to the OS, for running python

ROOT = os.path.realpath(os.path.dirname(__file__))

# first we run the simulation

#Define programs constants
pi=mypypm7.mypmmoda.pi
kb=mypypm7.mypmmoda.kb
mu0=mypypm7.mypmmoda.mu0
nmax=mypypm7.mypmmoda.nmax

#Initial problem conditions
x=0
t=0
rk4order=2
#Initialize output vectors
y=zeros(nmax)
yout=zeros(nmax)

# Define class Dast

class dforce:

# Define modules that will be part of the Dast class

# + Module init initializes the 	
    def __init__(self):    # defines the module with the input parameter self
        self.glade_ui = gtk.glade.XML(os.path.join(ROOT, 'superharm3c.glade'))    # defines the var self.glade_ui as the same in the est module gtk.glade
        self.glade_ui.signal_autoconnect(self)    # defines the var self.glade_ui as the same in the est module gtk.glade class
        self.window1 = self.glade_ui.get_widget('window1')    # defines the var self.glade_ui as the same in the est module gtk.glade class	
        self.window1.show()
        
#Buttons 	        
    def gtk_main_quit(self, widget):
	gtk.main_quit()
    def on_run_ex_clicked(self, widget):
	# Previous vars and parameters loading
	omega0a = self.glade_ui.get_widget("f01_sp").get_value()
	omegaa = self.glade_ui.get_widget("fd1_sp").get_value()
	q=mypypm7.mypmmoda.q = self.glade_ui.get_widget("q1_sp").get_value()
	kc=mypypm7.mypmmoda.kc = self.glade_ui.get_widget("kc1_sp").get_value()
	a0a = self.glade_ui.get_widget("a01_sp").get_value()	 	
	a0 = mypypm7.mypmmoda.a0 = a0a*1e-9
	emuestraa = self.glade_ui.get_widget("saym_sp").get_value()	
	emuestra=mypypm7.mypmmoda.emuestra=emuestraa*1000000.	# in Mpa now
	epuntaa = self.glade_ui.get_widget("ymt_sp").get_value()	
	epunta=mypypm7.mypmmoda.epunta=epuntaa*1000000000.	# in Gpa now
	rada = self.glade_ui.get_widget("tr_sp").get_value()
	rad = mypypm7.mypmmoda.rad = rada*1e-9
	a00a = self.glade_ui.get_widget("a00_sp").get_value()	
	a00 = mypypm7.mypmmoda.a00 = a00a*1e-9
	hama = self.glade_ui.get_widget("ham_sp").get_value()
	ham = mypypm7.mypmmoda.ham = hama*1.60217646e-19	
	eta = mypypm7.mypmmoda.eta = self.glade_ui.get_widget("mvis_sp").get_value()	
	epsilon1 = mypypm7.mypmmoda.epsilon1 = self.glade_ui.get_widget("lrd_sp").get_value()	
	epsilon2 = mypypm7.mypmmoda.epsilon2 = self.glade_ui.get_widget("cdis_sp").get_value()	
	delta = mypypm7.mypmmoda.delta = 0.0
	ljmina = self.glade_ui.get_widget("ljdep_sp").get_value()		
	ljmin = mypypm7.mypmmoda.ljmin = ljmina*1e-9
	lengtha = self.glade_ui.get_widget("ljlen_sp").get_value()	
	length = mypypm7.mypmmoda.length = lengtha*1e-9
	mtip = mypypm7.mypmmoda.mtip = self.glade_ui.get_widget("mmtip_sp").get_value()
	msample = mypypm7.mypmmoda.msample = self.glade_ui.get_widget("mmsam_sp").get_value()
	nper=mypypm7.mypmmoda.nper = self.glade_ui.get_widget("nper_sp").get_value()
	npp=mypypm7.mypmmoda.npp = self.glade_ui.get_widget("npp_sp").get_value()
	naux=int(nper*npp) 				# need C
	nperfin = mypypm7.mypmmoda.nperfin = self.glade_ui.get_widget("nperfin_sp").get_value()
	mypypm7.mypmmoda.n = naux
	dzca = self.glade_ui.get_widget("deltazc_sp").get_value()
	dzc=mypypm7.mypmmoda.dzc = dzca*1.e-9
	zcmina = self.glade_ui.get_widget("zcmin_sp").get_value()
	zcmin=mypypm7.mypmmoda.zcmin = zcmina*1.e-9
	zcmaxa = self.glade_ui.get_widget("zcmax_sp").get_value()
	zcmax=mypypm7.mypmmoda.zcmax = zcmaxa*1.e-9
	if emuestra > 0 and epunta > 0:		# NEW
		ebarra=1./(0.25/emuestra+0.25/epunta)		
	else:	
		ebarra=0.
	mypypm7.mypmmoda.ebarra=ebarra		# need epunta & ebarra
	f0a=kc*a0/q 				# need C
	mypypm7.mypmmoda.f0=f0a
	omega0=omega0a*2.*pi*1000.			# need C
	mypypm7.mypmmoda.omega0=omega0
	omega=omegaa*2.*pi*1000.			# need C
	mypypm7.mypmmoda.omega=omega
	periodo=2.*pi/omega 			# need C
	mypypm7.mypmmoda.periodo=periodo
	tmax=nper*periodo 			# need C
	mypypm7.mypmmoda.tmax=tmax
	tmin=mypypm7.mypmmoda.tmin=0.
	saltot=(tmax-tmin)/naux			# need C
	mypypm7.mypmmoda.saltot=saltot
	m=kc/omega0**2 			# need C
	mypypm7.mypmmoda.m=m
	nzc=int(1.99999+(zcmax-zcmin)/dzc)   # need C
	mypypm7.mypmmoda.nzc=nzc
	posini=mypypm7.mypmmoda.posini=0
	veloini=mypypm7.mypmmoda.veloini=0
	#=========================== Saving Sim Inputs ===============================
	tx00=[ 'f0 [Hz]', 'fd [Hz]', 'Q [adim]', 'Kc [N/m]', 'Amp [m]', 'Esample [Pa]', 'Etip [Pa]', 'Radtip [m]', 'IntDis [m]', 'HamK [J]', 'Eta [Pa.s]', 'Epsi1 [adim]', 'Epsi2 [adim]', 'LJmin [J]', 		'LJlength [nm]', 'Mtip [A.m2]', 'Msample [A.m2]', 'Nper [adim]', 'Npper [adim]', 'Nperss [adim]', 'DeltaZc [m]', 'Zcmin [m]', 'Zcmax [m]', 'FerritinRad [m]', 'Ebarra [Pa]', 'Tperiods [s]', 'Tmax 		[s]', 'PosIni [m]', 'VelIni [m]' ]
	tx00=np.array(tx00)
	tx01=[ omega0, omega, q, kc, a0, emuestra, epunta, rad, a00, ham, eta, epsilon1, epsilon2, ljmin, length, mtip, msample, nper, npp, nperfin, dzc, zcmin, zcmax, delta, ebarra, periodo, tmax, 	posini, veloini ]	# Converted to an unique array
	tx01=np.array(tx01)
	f1 = file("params.in","w")
	headings='Names, Values\n'
	f1.write(headings)
	for i in np.ndindex(tx00.shape):
		filerow='%s,%g\n' % (tx00[i],tx01[i])
		f1.write(filerow)
	f1.close()
 	#=========================Call fortran functions ===============================
	print naux
	mypypm7.mypmmoda.inivec(naux,nzc) # Initialize the fortran dynamic vectors (IN & OUT)
	print "llega1"
	mypypm7.mypmmoda.pmapd(t,saltot,posini,veloini,5) # Runs the point-mass Model
	print "llega2"
	#========================Obtaining outputs time domain ==========================
	# "t.s","z.nm","zp.deg"
	tarray=mypypm7.mypmmoda.tarray	
	xarray=mypypm7.mypmmoda.xarray
	varray=mypypm7.mypmmoda.varray
	farray=mypypm7.mypmmoda.farray
	tx1=array([tarray, xarray, varray, farray])	# Converted to an unique array
	tx1=tx1.T				#Array transpose
	np.savetxt('timedom.out', tx1)	# With the order 1 column (tarray)
	#========================Obtaining outputs zc domain ==========================
	# "zc.nm","amp1.nm","fase1.deg","dmin.nm","Fmax.nN","defl.nm","Fmedia.nN","ets.1e20J","vir.1e20J"
	zcaux=mypypm7.mypmmoda.zcaux
	amps=mypypm7.mypmmoda.amps
	fases=mypypm7.mypmmoda.fases
	dminimo=mypypm7.mypmmoda.dminimo
	interaux=mypypm7.mypmmoda.interaux
	defl=mypypm7.mypmmoda.defl
	fmediaux=mypypm7.mypmmoda.fmediaux
	ets1=mypypm7.mypmmoda.ets1
	virialaux=mypypm7.mypmmoda.virialaux
	tx2=array([zcaux, amps, fases, dminimo, interaux, defl, fmediaux, ets1, virialaux])	# Converted to an unique array
	tx2=tx2.T				#Array transpose
	np.savetxt('zcdom.out', tx2)		# With the order 1 column (zc)
	# DN ++++++
	mypypm7.mypmmoda.killvecs() # release the fortran dynamic vectors (IN & OUT) it is not necessarily related     

    def on_plot_ex_clicked(self, widget):
	timedom = np.genfromtxt('timedom.out')
	tarray = timedom[:,0]
	zcdom = np.genfromtxt('zcdom.out')
	zcaux = zcdom[:,0]
	if self.glade_ui.get_widget("single_cb").get_active() and self.glade_ui.get_widget("pm_cb").get_active():
		if self.glade_ui.get_widget("zt_cb").get_active():
			figure(1)			
			grid(True)
			xarray = timedom[:,1]
			#xarray=mypypm7.mypmmoda.xarray
			title('Trajectory (nm) vs. Time (s)')
			plot(tarray,xarray,'g.--')
			#show()
			savefig('fig1.png')
			if self.glade_ui.get_widget("ft_cb").get_active():
				figure(2)			
				grid(True)
				farray = timedom[:,3]
				title('Force (nN) vs. Time (s)')
				plot(tarray,farray,'b.--')
				#show()
				savefig('fig2.png')
				if self.glade_ui.get_widget("ampst_cb").get_active():
					figure(3)
					grid(True)			
					#amps=mypypm7.mypmmoda.amps
					amps = zcdom[:,1]					
					title('Amplitude (nm) vs. Average distance (nm)')
					plot(zcaux,amps,'r.--')
					#show()
					savefig('fig3.png')
					if self.glade_ui.get_widget("fasest_cb").get_active():
						figure(4)
						grid(True)			
						#fases=mypypm7.mypmmoda.fases
						fases = zcdom[:,2]						
						title('Phase (deg) vs. Average distance (nm)')
						plot(zcaux,fases,'g.--')
						#show()
						savefig('fig4.png')
						if self.glade_ui.get_widget("mdef_cb").get_active():
							figure(5)
							grid(True)			
							#defl=mypypm7.mypmmoda.defl
							defl = zcdom[:,5]
							title('Mean deflection (nm) vs. Average distance (nm)')
							plot(zcaux,defl,'b.--')
							#show()
						   	savefig('fig5.png')
							if self.glade_ui.get_widget("mdis_cb").get_active():
								figure(6)
								grid(True)			
								#dminimo=mypypm7.mypmmoda.dminimo
								dminimo = zcdom[:,3]
								title('Minimum distance (nm) vs. Average distance (nm)')
								plot(zcaux,dminimo,'r.--')
								#show()
							   	savefig('fig6.png')
								if self.glade_ui.get_widget("maforce_cb").get_active():
									figure(7)
									grid(True)			
									#fmediaux=mypypm7.mypmmoda.fmediaux
									interaux = zcdom[:,4]
									title('Maximum force (nm) vs. Average distance (nm)')
									plot(zcaux, interaux,'g.--')
									#show()
								   	savefig('fig7.png')
									if self.glade_ui.get_widget("etst_cb").get_active():
										figure(8)
										grid(True)			
										#ets1=mypypm7.mypmmoda.ets1
										ets1 = zcdom[:,7]
										title('Total dissipated energy (J) vs. Average distance (nm)')
										plot(zcaux,ets1,'b.--')
										#show()
									   	savefig('fig8.png')
										if self.glade_ui.get_widget("virialt_cb").get_active():
											figure(9)
											grid(True)			
											#virialaux=mypypm7.mypmmoda.virialaux
											virialaux = zcdom[:,8]
											title('Virial (J) S vs. Average distance (nm)')
											plot(zcaux,virialaux,'r.--')
										   	savefig('fig9.png')
											show()
										else:
										        show()
									else:
									        show()
								else:
								        show()
							else:
							        show()
						else:
						        show()
					else:
					        show()
				else:
				        show()
			else:
			        show()

		else:
		        show()
	else:
		print "Esto es SUPERHARM "
	#show()			# show the plottings
	#open('fig9.png')

    def on_save_ex_clicked(self, widget):
        print "Your simulation data and curves have been saved"
    def on_inputs_ex_clicked(self, widget):
        print "please upload the file in the following format:"

if __name__ == '__main__':
    dforce()
    gtk.main()
