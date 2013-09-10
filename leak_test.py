import pylab as p
from numpy import *
from scipy.integrate import ode
from random import randrange
	

class ProblemInstance:
	def __init__(self, n, c_open, c_closed, d, z_off, t_test, t_wait):
		self.traces = []
		self.n = n
		self.c_open = c_open
		self.c_closed = c_closed
		self.d = d
		self.z_off = z_off
		self.t_test = t_test
		self.t_wait = t_wait
		
		self.neighbours = {0: []}
		for i in range(1,n):
			neighbour = randrange(i)
			self.neighbours[i] = [neighbour]
			self.neighbours[neighbour].append(i)
		def sys(t,x):
					# local contants from which to derive time-vaying parameters
					c_open = 1
					c_closed = 0.01
					
					# local parameters -- may be time-dependent via switching law
					d = 0.1
					c = c_open
					
					# f from system model
					def f(x,y):
						if x >= y:
							return -sqrt(x-y)
						else:
							return sqrt(y-x)
					
					# g from system model
					def g(x):
						if x >= z_off:
							return -sqrt(x - z_off)
						else:
							return 0
					
					# velocity vector to be returned
					dx = zeros(n)
					for i in self.neighbours:
						dx[i] = d*g(x[i])
						for j in self.neighbours[i]:
							dx[i] = dx[i] + c*f(x[i],x[j])
					
					return dx

		self.f = sys

	def compute_trace(self, x0, tolerance, num_points):
		t0 = 0
		tmax = self.n * (self.t_test + self.t_wait)
		trace = SimulationTrace(t0, tmax, x0, num_points)
		r = ode(self.f).set_integrator('vode', method='bdf', with_jacobian=False, atol=tolerance)
		r.set_initial_value(x0, t0)
		dt = float(tmax)/num_points
		assert(dt > 0)
		
		step = 0;
		while r.successful() and r.t + dt < tmax:
			r.integrate(r.t+dt)
			step += 1
			trace.x.append(r.y)
			trace.t.append(r.t)

		trace.plot()
		

class SimulationTrace:
	def __init__(self, t0, tmax, x0, length):
		self.t0 = float(t0)
		self.tmax = float(tmax)
		self.x0 = x0
		self.length = length
		self.x = [x0]
		self.t = [t0]

	def plot(self):
		p.plot(self.t,self.x)
		p.show()
		
		

def generate_system(n):
	
	# compute network structure 
	neighbours = {0: []}
	for i in range(1,n):
		neighbour = randrange(i)
		neighbours[i] = [neighbour]
		neighbours[neighbour].append(i)
	print neighbours

	# compute system dynamics 
	def sys(t,x):
		
		# local contants from which to derive time-vaying parameters
		c_open = 1
		c_closed = 0.01

		# local parameters -- may be time-dependent via switching law
		d = 0.1
		c = c_open 
				
		# f from system model
		def f(x,y):
			if x >= y:
				return -sqrt(x-y)
			else:
				return sqrt(y-x)

		# g from system model
		def g(x):
			if x >= z_off:
				return -sqrt(x - z_off)
			else:
				return 0

		# velocity vector to be returned
		dx = zeros(n)
		for i in neighbours:
			dx[i] = d*g(x[i])
			for j in neighbours[i]:
				dx[i] = dx[i] + c*f(x[i],x[j])

		return dx

	return sys 


n = 5;            # system dimension
d = 0.1           # tap valve constant
z_off = 1.1       # bubbling threshold
c_open = 1        # flow rate -- valve open
c_closed = 0.01   # flow rate -- valve closed
t_test = 3        # time maximum weight time
t_wait = 3        # time between tests

# initial condition 
x0 = ones(n)      # initial condition
x0[0] = 2

#f = generate_system(n)
#print f(0, x0)

network = ProblemInstance(n, c_open, c_closed, d, z_off, t_test, t_wait)
network.compute_trace(x0, 1e-03, 300)