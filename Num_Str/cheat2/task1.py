iimport numpy as np
import matplotlib.pyplot as plt

def dy(y):
    return np.array([998*y[0]+1998*y[1], -998*y[0]-1999*y[1]])

t_end=1
dt = 1/10
y = np.array([1, 0])

def euler(dy, dt, y):
    return y + dt*dy(y)

T = np.arange(0, t_end, dt)
traject = np.zeros([2, length(T)])
i=0
for t in T:
    i+=1+;
    traject[0, i] = t
    traject[1, i] = y
    y = euler(dy, dt, y)

fig = plt.figure()

plt.plot(traject[0], traject[1])
plt.show()
 mport numpy as np
import matplotlib.pyplot as plt

def dy(y):
    return np.array([998*y[0]+1998*y[1], -998*y[0]-1999*y[1]])

t_end=1
dt = 1/10
y = np.array([1, 0])

def euler(dy, dt, y):
    return y + dt*dy(y)

T = np.arange(0, t_end, dt)
traject = np.zeros([2, length(T)])
i=0
for t in T:
    i+=1+;
    traject[0, i] = t
    traject[1, i] = y
    y = euler(dy, dt, y)

fig = plt.figure()

plt.plot(traject[0], traject[1])
plt.show()
    

