
class atom:
	def __init__(self, x0, v0, masa, czy_nieruchomy=False):
		self.x=x0
		self.v=v0
		self.xnext=x0
		self.vnext=v0
		self.masa=masa
		self.czy_nieruchomy=czy_nieruchomy
	def licz_sile(self, G, atomy)
		for at in atomy:
			if at!=self:
				r2=(self.x[0]-at.x[0])**2+(self.x[1]-at.x[1])**2
				return G*at.masa*self.masa/r2
# xn,xn-1,vn,vn-1 - wektory dwuwymiarowe
# m,Fn wypadkowa,Fn-1, czy_nieruchomy, (Ekin,Epot)	

class symulacja:
	def __init__(self, atomy, metoda_integracji, dt, G):
		self.atomy=atomy
		self.metoda_integracji=metoda_integracji
		#self.metody={'E':  , 'V':  , 'VV':  , 'LF':  }
		self.dt=dt
		self.G=G
	def dodaj_atom(self, atom):
		self.atomy.append(atom)
	def wykonaj_krok(self):
		for at in self.atomy:
			if metoda_integracji=='E':
				F=at.licz_sile(self.G,self.atomy)
				at.xnext = at.x+ at.v*self.dt+ F/(2*at.masa)*self.dt**2
				at.vnext = at.v+ F/at.masa*self.dt
				print(at.xnext,at.vnext)
# rysuj
#M=500, m=0.1
at1=atom([0,0],[0,0],500,True)
at1=atom([2,1],[1,1],0.1)
S=symulacja([at1,at2],'E',0.001,0.01)
