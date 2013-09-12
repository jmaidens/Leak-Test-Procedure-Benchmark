import pylab as p
from numpy import *
from scipy.integrate import ode
from random import randrange
	

class ProblemInstance:
	def __init__(self, n, c_open, c_closed, d, z_off, t_init, t_test, t_wait):
		self.traces = []
		self.n = n
		self.c_open = c_open
		self.c_closed = c_closed
		self.d = d
		self.z_off = z_off
		self.t_init = t_init
		self.t_test = t_test
		self.t_wait = t_wait
		
		# ProblemInstance has an attribute neighbours that contains a dictionary of the neighbours of each segment 
		self.neighbours = {0: [1], 1:[0]}
		for i in range(2,n):
			neighbour = randrange(1,i)
			self.neighbours[i] = [neighbour]
			self.neighbours[neighbour].append(i)
		
		# there is a dictionary containing the distance from the root of each segment in the network
		self.distance = {0: 0}
		for i in range(1,n):
			list = [x for x in self.neighbours[i] if x<i]
			assert(len(list) == 1)
			j = list[0]
			self.distance[i] = self.distance[j] + 1

		self.depth = 0
		for i in self.neighbours:
			self.depth = max((self.depth, self.distance[i]))

		
	def compute_dynamics(self):
		# The system dynamics are defined by sys and stored in self.f
		def sys(t,x):
			
			# local contants from which to derive time-vaying parameters
			c_open = 1
			c_closed = 0.01
			d_open = 0.2
			d_closed = 0
			
			
			# dictionaries of local parameters -- may be time-dependent via switching law
			c = {}       # flow rate through internal valves
			cplus = {}   # flow rate through leaf valves
			d = {}       # flow rate through bubbler valves
			
			# case t < t_init
			if 0 <= t <= self.t_init:
				for i in self.neighbours:
					if isleaf(i, self.neighbours):
						cplus[i] = c_closed
					c[i] = c_open
					d[i] = d_closed
			
			# case t > t_init
			else:
				# threshold over which bubbler valves should be open
				distance_threshold = self.depth - (t - t_init - t_wait)/(t_test + t_wait)
				for i in self.neighbours:
					if isleaf(i, self.neighbours):
						cplus[i] = c_closed
					c[i] = c_closed
				#####  CHECK IF THE FOLLOWING LINE IS CORRECT #######
					if self.distance[i] > distance_threshold:
						d[i] = d_open
					else:
						d[i] = d_closed
			
						
			'''			
			######## TEMPORARY SO THAT THIS STILL WORKS !!!!!!!!!!
			d = d_open
			c = c_open
			######## END TEMPORARY
			'''
		
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
			'''
			# velocity vector to be returned
			dx = zeros(n+1)
			for i in self.neighbours:
				if i == 0 or i == self.n:
					dx[i] = 0
				else:
					dx[i] = d*g(x[i])
					for j in self.neighbours[i]:
						dx[i] = dx[i] + c*f(x[i],x[j])
					# add link to exterior if i is a leaf
					if isleaf(i, self.neighbours):
						dx[i] = dx[i] + c_closed*f(x[i], x[self.n])
			dx[self.n] = 0
			'''
			
			# velocity vector to be returned
			dx = zeros(self.n+1)
			for i in self.neighbours:
				assert(i < self.n)
				if i == 0:
					dx[0] = 0    # the constant state x[0] represents the input to the system
					dx[self.n] = 0  # the constant state x[n] represents the external pressure
				else:
					dx[i] = d[i]*g(x[i])
					for j in self.neighbours[i]:
						if j > i:
							dx[i] = dx[i] + c[j]*f(x[i],x[j])
						else:
							dx[i] = dx[i] + c[i]*f(x[i],x[j])
					# add link to exterior if i is a leaf
					if isleaf(i, self.neighbours):
						dx[i] = dx[i] + cplus[i]*f(x[i], x[self.n])
			
			

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

		self.traces.append(trace)
		trace.plot()

	def print_graph(self):
		
		print "neighbours = ", self.neighbours
		print "distances =", self.distance
		for i in self.neighbours:
			print "isleaf(", i, ") =", isleaf(i, self.neighbours)



class SimulationTrace:
	def __init__(self, t0, tmax, x0, length):
		self.t0 = float(t0)
		self.tmax = float(tmax)
		self.x0 = x0
		self.length = length
		self.x = [x0]
		self.t = [t0]

	def plot(self):
		ax = p.plot(self.t,self.x)
		p.legend(ax)
		p.axis([0, 30, 0, 3])
		p.show()


# The helper function isleaf checks whether node i is a leaf of the tree neighbours
def isleaf(i, neighbours):
	if i <= 1:
		return False
	elif len(neighbours[i]) == 1:
		return True
	else:
		return False


n = 5;            # system dimension
d = 0.1           # tap valve constant
z_off = 1.1       # bubbling threshold
c_open = 1        # flow rate -- valve open
c_closed = 0.01   # flow rate -- valve closed
t_init = 10       # time to pressurize network initially 
t_test = 3        # time to test 
t_wait = 3        # time between tests

# initial condition 
x0 = ones(n+1)      # initial condition
x0[0] = 2           # pressure of gas source 
x0[n] = 1		   # ambient pressure 

#f = generate_system(n)
#print f(0, x0)

network = ProblemInstance(n, c_open, c_closed, d, z_off, t_init, t_test, t_wait)
# network.print_graph()
network.compute_dynamics()
network.compute_trace(x0, 1e-03, 300)


