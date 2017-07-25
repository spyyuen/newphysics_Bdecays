'''

description:

'''
__author__    = "Stephanie Yuen"
__email__     = "stephanie.yuen@cern.ch"
__created__   = "2017-07-25"
__copyright__ = "Copyright 2017 Stephanie Yuen"
__license__   = "GPL http://www.gnu.org/licenses/gpl.html"


###################### Modules
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
import pylab
import math

###################### Global constants from Particle Data Group (http://pdg.lbl.gov/)
G_F		= 1.16637*(10**-5) 	# Fermi Coupling Constant G_F/(hbar*c)^3 [GeV^-2]
hbar		= 6.58212*(10**-25)	# Planck's constant/(2pi)
Vcb_mag		= 0.04			# Magnitude of Vcb CKM matrix element
tau_B		= 1.52*(10**-12)	# Lifetime of B meson [s]
m_B		= 5.28 			# Mass of B0 meson [GeV]
m_D		= 1.87 			# Mass of D+/- meson [GeV]
m_tau		= 1.776 		# Mass of tau lepton [GeV]
m_ratio_cb	= 1./3.			# Approximate ratio of masses of c and b quarks 


###################### Global functions
def plotFunctions(x,y=[],xlabel=None,ylabel=None,title=None):
	colors = ['i','g','r','c','m','b']
	# compose plot
	fig, ax = plt.subplots()
	if result.debug:
		for y_array in y:
			for xval, yval in zip(x,y_array):
				print 'xval ', xval, '\t yval ', yval
	for y_i in y:
		ax.plot(x,y_i,'o',markersize=1)#, c=colors[i])
	if title: plt.title(title)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	#plt.legend()
	plt.show() # show the plot

class decayRate(object):
	def __init__( self ):
		# Initialize for numeric integration of differential decay rate
		self.x 		= np.array([ 0.0, math.sqrt(3./5.), -math.sqrt(3./5.) ])
		self.w		= np.array([ 8./9., 5./9., 5./9. ])
		self.limit_hi	= ( m_B-m_D )**2
		self.limit_lo	= 0
		self.rho	= 0.2  #Constant defined for form factor
		self.DeltaS 	= 1.21 #Constant defined in second form factor to study diff between decay rates of m_l approx 0 and m_l=m_tau
		self.ratiotanb_higgspl = 0.45 #Technically tan(beta)/m(Higgs+)

	#Magnitude of D-meson momentum in the B-meson rest frame
	def pD_tilda( self, q2_i ):
		valInSqrt = ((m_B**2+m_D**2-q2_i)**2)/(4*m_B**2) - m_D**2
		if valInSqrt<0: 
			valInSqrt = 0
			print 'warning: negative in sqrt of pD_tilda'
		return math.sqrt( valInSqrt )

	def w_informfactor(self, q2_i):
		return (m_B**2 + m_D**2 - q2_i) / (2*m_B*m_D)

	#First form factor
	def fplus( self, w, q2_i ):
		#return 1-self.rho*(((m_B**2+m_D**2-q2_i)/(2*m_B*m_D))**2-1)		
		return 1 - (w**2-1)*self.rho

	#Second form factor
	def f_0( self, w, q2_i ):
		#print 'self.DeltaS ', self.DeltaS, '\t self.rho ', self.rho
		#return self.DeltaS * ( 1 - (((m_B**2+m_D**2-q2_i)/(2*m_B*m_D))**2 -1 )*self.rho)
		#Study difference between decay rates with m_l approx 0 and m_l=m_tau
		return self.DeltaS * self.fplus( w, q2_i )

	#Form factor introducing extended Higgs sector needed in supersymmetry theories
	def f_0_newphysics( self, q2_i ):
		return 1 - (self.ratiotanb_higgspl**2) *(q2_i/(1-m_ratio_cb))

	#Differential decay rate with full mass effects
	#for m_l!=0, evaluates at one q2_i=q^2 val and returns one y val
	def ddr_masseffects(self, q2_i, newPhysics=False): 
		m_l 		= m_tau
		self.limit_lo 	= m_l**2	
		w 		= self.w_informfactor(q2_i)
		prefactor 	= (G_F**2)*(Vcb_mag**2)/(24*(math.pi**3))*((self.pD_tilda(q2_i))**3)
		secondFactor	= (1-((m_l**2)/q2_i))**2
		firstTerm 	= (self.fplus(w,q2_i))**2 * (1+(m_l**2)/(2*(q2_i))) 
		secondTerm	= 3./2.*((m_l**2)/(q2_i)) * (m_B**2-m_D**2)/(2*m_B*self.pD_tilda(q2_i))
		if newPhysics:
			f_0 = self.f_0(w,q2_i) * self.f_0_newphysics(q2_i)
		else:
			f_0 = self.f_0(w, q2_i)
			if result.debug:
				print '============ ', q2_i
				print 'f+ =\t ', self.fplus(w,q2_i)
				print 'f0 =\t ', f_0
				print 'pD =\t ', self.pD_tilda(q2_i)
				print 3./2.*((m_l**2))*(m_B**2-m_D**2)/(2*m_B)
				print '============'
		return prefactor * secondFactor * (firstTerm + secondTerm*(f_0**2) )

	#Differential decay rate for m_l=0
	#for m_l==0, evaluates at one q2_i=q^2 val and returns one y val
	def ddr( self, q2_i ): 
		w 	= self.w_informfactor(q2_i)
		rho	= 0.2
		fplus 	= self.fplus( w, q2_i )
		pD_tilda= self.pD_tilda(q2_i)
		f  	= (G_F**2)*(Vcb_mag**2)/(24*(math.pi**3))*(pD_tilda**3)*(fplus**2)
		return f
	
	def numericIntegration( self, m_l=0, newPhysics=False ):
		#Compute decay rate via numeric integration
		#Solve int^b_a f(x)dx and evaluate at x_i and w_i
		if result.debug: print 'Computing decay rate for m_l = ', m_l, ' case.'
		if m_l == 0: 	self.limit_lo = 0
		else:		self.limit_lo = m_tau**2
		coeff = (self.limit_hi-self.limit_lo)/2
		summation = 0
		for i in range(0,self.x.size):
			w_i 	= self.w[i]
			x_i 	= self.x[i]
			if m_l == 0: 
				q2_i 	= (self.limit_hi-self.limit_lo)/2*x_i + (self.limit_hi+self.limit_lo)/2
				
				if result.debug: print 'x_i ', x_i, '\tq2_i ', q2_i
				summation += w_i * self.ddr(q2_i)
				if result.debug: print 'integrating ', self.limit_lo, ' to ', self.limit_hi
				if result.debug: print 'w_i ', w_i , '\t ddr ', self.ddr(q2_i)
			else: 
				q2_i 	= (self.limit_hi-self.limit_lo)/2*x_i + (self.limit_hi+self.limit_lo)/2
				if result.debug: print 'x_i ', x_i, '\tq2_i ', q2_i
				if result.debug: print 'integrating ', self.limit_lo, ' to ', self.limit_hi
				summation += w_i * self.ddr_masseffects(q2_i, newPhysics=newPhysics)
		return coeff*summation

	def plotDecayRates( self ):
		#plot decay rate for m_l approx 0
		x = np.linspace(0,11,111) # 110 linearly spaced numbers
		y = np.array([])
		self.limit_lo	= 0
		for x_i in x:
			y_i = self.ddr(x_i)
			y = np.append(y, y_i)
		
		#plot decay rate with full mass effects
		#Plot starting from min(q^2) approx 3.16.  Below this, plot 0
		qsq_min = m_tau**2 #[GeV^2]
		y_f = np.array([])
		x_f_newphysics = np.linspace(0,11,111) # 110 linearly spaced numbers
		y_f_newPhysics = np.array([])
		self.limit_lo	= m_tau**2
		self.DeltaS 	= 1.21 
		self.ratiotanb_higgspl = 0.45 

		for x_i in x:
			if x_i < qsq_min: 
				y_i = 0
				y_i_newPhysics = 0
			else:	
				self.limit_lo=m_tau**2
				y_i = self.ddr_masseffects( x_i, newPhysics=False )
				y_i_newPhysics= self.ddr_masseffects( x_i, newPhysics=True)
			y_f = np.append(y_f, y_i)
			y_f_newPhysics = np.append(y_f_newPhysics, y_i_newPhysics)
		
		if result.plot:
			plotFunctions(x,y=[y,y_f,y_f_newPhysics], xlabel='$q^2$ [GeV$^2$]', ylabel="$\\frac{d\Gamma}{dq^2}$ [GeV$^{-1}$]")
		if result.plot0:
			plotFunctions(x,y=[y], xlabel='$q^2$ [GeV$^2$]', ylabel="$\\frac{d\Gamma}{dq^2}$ [GeV$^{-1}$]")

		if result.plottau:
			plotFunctions(x,y=[y,y_f], xlabel='$q^2$ [GeV$^2$]', ylabel="$\\frac{d\Gamma}{dq^2}$ [GeV$^{-1}$]")
			
		if result.plotnp:
			plotFunctions(x,y=[y,y_f,y_f_newPhysics], xlabel='$q^2$ [GeV$^2$]', ylabel="$\\frac{d\Gamma}{dq^2}$ [GeV$^{-1}$]")

	def plotR_D( self ):
		x_dS 	= np.linspace(0,2,100)
		#Calculate R_D as function of Delta S
		y_R_D	= np.array([])
		DeltaS_init = self.DeltaS
		for x_i in x_dS:
			self.DeltaS	= x_i
			Gamma 		= self.numericIntegration()
			Gamma_mtau 	= self.numericIntegration(m_l=m_tau)
			R_D_i		= Gamma_mtau/Gamma
			y_R_D 		= np.append( y_R_D, R_D_i )
		plotFunctions(x_dS,y=[y_R_D], xlabel='$\Delta S$', ylabel='R(D)')
		#Reset DeltaS
		self.DeltaS=DeltaS_init

		#Calculate R_D as function of tan(beta)/m_H+
		x_ratiobeta_chargedHiggs= np.linspace(0.,.5,100)
		y_R_D_newphysics	= np.array([])
		ratiotanb_higgspl_init	=self.ratiotanb_higgspl
		for x_i in x_ratiobeta_chargedHiggs:
			self.ratiotanb_higgspl = x_i
			Gamma 		= self.numericIntegration(m_l=0, newPhysics=False)
			Gamma_mtau 	= self.numericIntegration(m_l=m_tau, newPhysics=True)
			R_D_i		= Gamma_mtau/Gamma
			y_R_D_newphysics 		= np.append( y_R_D_newphysics, R_D_i )
		plotFunctions(x_ratiobeta_chargedHiggs,y=[y_R_D_newphysics], xlabel="$\\frac{tan\\beta}{m_{H^+}}$", ylabel='R(D)')
		#Reset tanbeta/mH+S
		self.ratiotanb_higgspl=ratiotanb_higgspl_init


	def execute( self ):
		#Calculate decay rate for m_l=0 case
		self.limit_lo	= 0
		Gamma = self.numericIntegration()
		#Total decay width of the B meson = hbar/lifetime_B
		Gamma_B = hbar/tau_B
		print 'm_l=0 case: Gamma/Gamma_B = ', Gamma, ' / ', Gamma_B, ' = ', Gamma/Gamma_B

		print ''	
	
		self.limit_lo	= m_tau**2
		#Calculate decay rate for m_lm_tau0 case
		Gamma_mtau = self.numericIntegration(m_l=m_tau, newPhysics=False)
		R_D = Gamma_mtau/Gamma
		print 'R_D = ', Gamma_mtau, ' / ', Gamma, ' = ', R_D
		

###################### main
if __name__ == '__main__':

	parser = ArgumentParser()
	# General switches
	parser.add_argument( '-p', '--plot', help='Draw and show all plots', action='store_true' )
	parser.add_argument( '-p0', '--plot0', help='Draw and show plots for ml=0', action='store_true' )
	parser.add_argument( '-ptau', '--plottau', help='Draw and show plots ml=tau', action='store_true' )
	parser.add_argument( '-pnp', '--plotnp', help='Draw and show plots ml=tau with new physics', action='store_true' )
	parser.add_argument( '-prd', '--plotrd', help='Draw and show plots R(D)', action='store_true' )
	parser.add_argument( '-d', '--debug', help='Run with debug option', action='store_true' )
	result = parser.parse_args()

	dr 	= decayRate()
	if result.plot or result.plot0 or result.plottau or result.plotnp: 
		#Plot overlay of decay rates for different m_l
		dr.plotDecayRates()
	if result.plotrd:
		dr.plotR_D()
	dr.execute()

