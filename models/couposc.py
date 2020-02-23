# couposc.py
# Simulate a discrete collection of oscillators, as in problem 5-13.
# We will use this as a model of a vibrating string.
# $Id: couposc.py,v 1.8 2008/03/24 17:23:04 doughera Exp $

from visual.graph import *

# Physical Constants.
N = 100			# Number of coupled oscillators.
A = 1.0			# Overall amplitude of motion
FT = 1.0		# Force of Tension
ell = 1.0               # length of string between masses.  Don't change this!
L = (N + 1) * ell       # Total length of the system
m = 1.0			# mass of each oscillator
M = N * m               # Total mass
mu = M/L	        # Mass per unit length

# Derived physical constants -- FIX THESE TO BE THE CORRECT VALUES!
# Use formulas so they will stay sensible even as you change N above.
omega0 = 1.0		# Natural angular frequency of undamped oscillator
omega1 = 1.0		# Angular frequency of 1st mode.

T0 = 2 * pi / omega0    # Convenient time scale for simulations.
T1 = 2 * pi / omega1	# Period of first normal mode.

# Time step.  Adjust as needed to achieve the desired precision.
dt = T0/80.0

# Graph scale
ymax = 1.0*A

# To model 'N' oscillators, we will use N+2 points, numbered
# 0, 1, 2, 3, ... N+1.  Points 0 and N+1 are actually the boundaries.
# We will keep them fixed, but adding them in as if they were
# masses makes the programming easier.  (In particular, it makes the
# acceleration calculation easier.)
# To iterate over the particles, we'll use the range(1,N+1) command, which
# prints out the integers from 1 to N.

# Set up some global arrays to use as placeholders.
# The zeros() command creates an array full of zeros.  The "Float"
# argument explicitly informs Python that these arrays will contain
# floating point numbers.  (The default would be to assume integers.)

y = zeros(N+2, Float)           # y[p] is the position of particle p.
v = zeros(N+2, Float)		# v[p] is the velocity of particle p.
a = zeros(N+2, Float)   	# a[p] is the acceleration of particle p.
ynew = zeros(N+2, Float)	# New position after integration.
vnew = zeros(N+2, Float)	# New velocity after integration.

# Acceleration function.  Adjust as necessary.
def accel(y, v, t):
    "accel(y, v,t): acceleration for the simple harmonic oscillator."
    # This next line informs Python that we're using the global array
    # a[], and that it shouldn't make its own local private copy.
    global a
    # Fix the left and right-hand sides.
    a[0] = 0.0
    a[N+1] = 0.0

    # This version loops explicitly over all the particles.
    for p in range(1,N+1):
        a[p] = -omega0**2 * (y[p] - y[p-1]) - omega0**2 * (y[p] - y[p+1])
    #
    # This version has Python do the loop implicitly.  It's faster.
    # Note how the subscripts are adjusted to take the particles to
    # the left or right as appropriate.
    # a[1:N+1] = -omega0**2 *(y[1:N+1]-y[0:N]) - omega0**2 *(y[1:N+1]-y[2:N+2])
    #
    return a

# 4th Order Runge-Kutta
# This is essentially unchanged from the previous exercises.
# Even though we are passing in arrays giving the position and velocity of
# the entire string, Python automatically operates on the entire
# arrays all at once.  This is one of the nice things about VPython.
def rk4(y, v, t, dt):
    "rk4(y,v,t,dt): Advance one time step using 4th-order Runge Kutta."
    yk1 = dt * v
    vk1 = dt * accel(y, v, t)
    yk2 = dt * (v + vk1/2.0)
    vk2 = dt * accel(y + yk1/2.0, v + vk1/2.0, t + dt/2)
    yk3 = dt * (v + vk2/2.0)
    vk3 = dt * accel(y + yk2/2.0, v + vk2/2.0, t + dt/2)
    yk4 = dt * (v + vk3)
    vk4 = dt * accel(y + yk3, v + vk3, t + dt)
    ytmp = y + (yk1 + 2.0 * yk2 + 2.0 * yk3 + yk4) / 6.0
    vtmp = v + (vk1 + 2.0 * vk2 + 2.0 * vk3 + vk4) / 6.0
    return [ytmp, vtmp]


####################################################
# Initial conditions.
# Edit these as appropriate.
####################################################

####################################################
# Initial positions.
####################################################
y[0] = 0 		# Fix the left-hand side at zero.
y[N+1] = 0 		# Fix the right-hand side at zero.
# The range(1,N+1) command only prints out [1,2,3, ... N].
for p in range(1, N+1):    # p is particle number.
    # There are several different ways to set y[p].  Here are a few:
    # Set all y's to zero and change initial velocities below.
    # y[p] = 0

    # Or, use equation 5-26 in the text.  This example uses mode 3.
    y[p] = A * sin(3 * p * pi /(N+1.0))

    # Or, explicitly set only certain particle positions (e.g. 2).
    #if p == 2:
    #	y[p] = 1
    #else:
    #	y[p] = 0

    # Or, start with a flat string.
    # y[p] = 0

####################################################
# Initial velocities.
####################################################
v[0]   = 0               # The left and right boundaries are
v[N+1] = 0               # clamped and don't move.
# This version sets them all the particle velocities to zero.
for p in range(1, N+1):
    v[p] = 0

# Or, you can give one specific particle (e.g. p = 2) a kick.
# v[2] = 1

####################################################
# Set up for the plot
####################################################
mygraph = gdisplay(xtitle='x', ytitle='y', title='Coupled Oscillators',
		    ymin = -ymax, ymax = ymax)
string = gcurve(color=color.cyan)

# Create the curve using our initial positions:
for p in range(0,N+2):
    string.plot(pos=(p, y[p]))

# Put a time ticker in the upper-right-hand corner of the graph.
clock = label(pos=(N,A), text='0', display=mygraph.display)

####################################################
# Integrate forward for 100 periods.
####################################################
t = 0
while t < 100.0* T1:
    rate(1000)		# Adjust this as needed for a nice display.
    [ynew, vnew] = rk4(y, v, t, dt)
    y = ynew
    v = vnew
    # Adjust the string positions and velocities to the new values.
    # (Rather than adding new points to the curve, this command just
    # updates the positions of the old points.)
    string.gcurve.pos[0:N+2,1] = y[0:N+2]
    clock.text = '%7.3f' % t		# Update the clock.
    t = t + dt

print "Done."

