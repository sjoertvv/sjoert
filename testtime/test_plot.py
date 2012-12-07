import matplotlib
from matplotlib import rc
from matplotlib import pyplot as plt
import numpy as np
import sjoert

def plotsome():
    xx = np.logspace(1, 2)
    for k in [1,0.5,2]:
        plt.plot(xx, (xx-k)**k, label=str(k), alpha=k/2.)
        plt.scatter(xx, xx**(k+1), label='s'+str(k))

    plt.xlabel('Time')
    plt.title('number 10,000')

def plotsomelatex():
    xx = np.logspace(1, 2)
    for k in [1,0.5,2]:
        plt.plot(xx, (xx-k)**k, label=str(k), alpha=k/2.)
        plt.scatter(xx, xx**(k+1), label='s'+str(k))

    plt.legend()
    plt.xlabel('Time y $(y \\Gamma^6 \\beta^2)$')
    plt.ylabel('A Label ($x^2$)')
    plt.title('7$=7$,1=$1$, large number $10^5$,  text ${\\rm text}$ ')
    plt.yscale('log')




plotsome()
print 'done plotting'
key = raw_input('next:latex')
plt.close()

plotsomelatex()
print 'done plotting'
key = raw_input('next:normal')

plt.savefig('./default.pdf')
print 'done saving PDF'
plt.savefig('./default.ps')
print 'done saving PS'
key = raw_input('next: plot.init()')
plt.close()


sjoert.plot.init(golden=True)
#rc('font',family='sans-serif')
print 'done with init'

plotsome()
print 'done plotting'
key = raw_input('next:latex')


plt.close()
plotsomelatex()
print 'done plotting'

plt.savefig('./sjoertinit.pdf')
print 'done saving PDF'
plt.savefig('./sjoertinit.ps')
print 'done saving PS'


