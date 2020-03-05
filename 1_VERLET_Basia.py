# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
class atom:
	def __init__(self, x0, v0, mass, mobile=True):
		self.x=x0
		self.v=v0
		self.xnext=x0
		self.vnext=v0
		self.mass=mass
		self.mobile=mobile
	def count_force(self, G, atoms):	# na razie tylko dla 2 ciał
		for at in atoms:
			if at!=self:
				r2=sum((self.x-at.x)**2)
				F=G*at.mass*self.mass/r2
				return -F*(self.x-at.x)/sqrt(r2)
# xn,xn-1,vn,vn-1 - wektory dwuwymiarowe
# m,Fn wypadkowa,Fn-1, mobile, (Ekin,Epot)	

class simulation:
	def __init__(self, atoms, integration_method, dt, G):
		self.atoms=atoms
		self.integration_method=integration_method
		#self.metody={'E':  , 'V':  , 'VV':  , 'LF':  }
		self.dt=dt
		self.G=G
	def do_step(self):
		X=np.empty([2,len(self.atoms)])
		V=np.empty([2,len(self.atoms)])
		for i in range(len(self.atoms)):
			at=self.atoms[i]
			if at.mobile:
				if self.integration_method=='E':
					F=at.count_force(self.G,self.atoms)
					at.xnext = at.x+ at.v*self.dt+ F/(2*at.mass)*self.dt**2
					at.vnext = at.v+ F/at.mass*self.dt
					at.x,at.v=at.xnext,at.vnext
			X[:,i]=at.x
			V[:,i]=at.v
		return X,V
# rysuj
#M=500, m=0.1
at1=atom(np.array([0,0]),np.array([0,0]),500,False)
at2=atom(np.array([-0.7,0]),np.array([-1,2]),0.1)
S=simulation([at1,at2],'E',0.001,0.01)

plt.axis([-1, 1, -1, 1.5])

N=2000
for i in range(N):
    X,V=S.do_step()
    plt.scatter(x=X[0],y=X[1],c='black')
    #plt.pause(1e-10/(V[0,1]*V[1,1]))

plt.show()
# for i in range(10):
# 	X,V=S.wykonaj_krok()
# 	print(X)
