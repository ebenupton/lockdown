from numpy import arange, sin, pi
from numpy.random import lognormal
from math import log, sqrt
from matplotlib import pyplot
from matplotlib.ticker import MaxNLocator
import sys

# compute params of underlying normal distribution from lognormal parameters

def ln_params(mean, sd):
    return log(mean*mean/sqrt(mean*mean+sd*sd)), sqrt(log(1+sd*sd/(mean*mean)))
    
# params   
    
INCUBATION_MU, INCUBATION_SIGMA = ln_params(5, 3.9)
ONSETDEATH_MU, ONSETDEATH_SIGMA = ln_params(20.2, 11.6)

# run a simulation with specified pre- and post-lockdown R values
# and either instant (Cheianov) or delayed infectiousness

def simulate(R_pre, R_post, instant):
    print ".",
    sys.stdout.flush()

    # empirical tau values to get ~3 day doubling

    tau = 10 if instant else 5

    # counts of infectious people and deaths per day
    
    infect = [0 for i in range(1000)]
    deaths = [0 for i in range(1000)]

    # remember how many new cases we created each day
    
    cases = [0 for i in range(50)]

    # dirty startup with 1000 cases at March 1
    # ringing appears to dissipate by day 10
    # TODO: check nothing horrid lurking here

    infect[0] = 1000

    # run for 50 days

    for t in range(50):
        # number of new cases to create: spread R evenly across tau days
        
        cases[t] = int(infect[t] * (R_pre if t < 23 else R_post) / tau)

        # create the cases
        
        for i in range(cases[t]):
            # pick infection-to-onset and onset-to-death from lognormal distributions
            # TODO: surely these aren't independent variables
            
            onset = lognormal(INCUBATION_MU, INCUBATION_SIGMA)
            death = lognormal(ONSETDEATH_MU, ONSETDEATH_SIGMA)

            # increment tau infectious people counters starting either at
            # t+1 (Cheianov model), or at t+onset-1
            
            for j in range(tau):
                infect[t+max(1, int(0 if instant else onset)-1)+j] += 1

            # increment appropriate death counter
            
            deaths[t+int(onset+death)] += 1

    # return day of peak, and per-day case counts normalised to maximum

    return max([(deaths[i], i) for i in range(1000)])[1], [float(n)/max(cases) for n in cases]

# plot a cases line chart

def plot_cases(cases, filename):
    # invent plausible value to cover dirty startup

    cases[0] = cases[3]/2

    fig, ax = pyplot.subplots()
    ax.plot(arange(0, 50, 1), cases)

    ax.set(xlabel='Days from March 1', ylabel='Daily infections normalized to maximum')
    ax.grid()

    fig.savefig(filename+".png")

# plot a "day of peak" histogram

def plot_histo(peaks, filename):
    fig, ax = pyplot.subplots()

    ax.hist(peaks, bins=arange(40, 55)-0.5)
    
    ax.set(xlabel='Days from March 1', ylabel='Frequency')    
    
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_ticks(arange(40, 55, 1))

    fig.savefig(filename+".png")

# plot cases with delayed infectiousness
    
peak, cases = simulate(2.75, 0.65, True)
plot_cases(cases, "cases_instant")

print

# plot cases with instant infectiousness

peak, cases = simulate(2.75, 0.65, False)
plot_cases(cases, "cases_delayed")

print

# plot histogram for R_post = 0.65

peaks = [simulate(2.75, 0.65, True)[0] for i in range(50)]    
plot_histo(peaks, "histo_p65")

print

# plot histogram for R_post = 0.80

peaks = [simulate(2.75, 0.80, True)[0] for i in range(50)]    
plot_histo(peaks, "histo_p80")

print
