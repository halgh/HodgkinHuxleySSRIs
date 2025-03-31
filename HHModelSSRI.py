'''
    This was an old project from a module in my second year undergrad degree.
    
    It simulates the Hodgkin-Huxley neuron potential model, with the addition
    of visualising change in potential under the effect of SSRIs such as citalopram etc.
'''


from __future__ import division
from curses import window
from matplotlib import pylab as pl
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import numpy as np
from tkinter import *
import scipy.signal

class test:
    def __init__(self):
        self.flouroxMg = 1.0
        self.setralineMg = 1.0
        self.citalopramMg = 1.0
    def start(self):
        print("l")

obj = test()

def HHmodel(inhibition):
    # defining main functions
    alpha_n = lambda v: 0.01*(-v+10)/(np.exp((-v+10)/10) - 1) if v != 10 else 0.1
    beta_n = lambda v: 0.125*np.exp(-v/80)
    n_inf = lambda v: alpha_n(-v)/(alpha_n(-v)+beta_n(-v))

    alpha_m=lambda v: 0.1*(-v+25)/(np.exp((-v+25)/10) - 1 ) if v!=25 else 1
    beta_m= lambda v: 4*np.exp(-v/18)
    m_inf = lambda v: alpha_m(-v)/(alpha_m(-v)+beta_m(-v))

    alpha_h=lambda v: 0.07*np.exp(-v/20)
    beta_h= lambda v: 1/(np.exp((-v+30)/10)+1)
    h_inf = lambda v: alpha_h(-v)/(alpha_h(-v)+beta_h(-v))

    # Here we use the original Hogkin-Huxley (1950's) paper for these values
    #g_n = 120 #OG = 120
    g_k = 36
    g_l = 0.3
    v_n = 115 * (inhibition/100)
    v_k = -12 * (inhibition/100)
    v_l = 10.613
    c = 1

    # Calculate the length of the simulation with a 'step' (dt) of 0.025
    dt = 0.025
    # Return evenly spaced values between x and n with a 'step' (x, n, dt)
    # in this case, 0 to 500 with a step of 0.025 should give us a sim length of 2000ms
    time = np.arange(0,500,dt)
    # Create an array of given length Time, fill with 0s
    v = np.zeros(len(time))

    # Calculate rest values, and place into v[0]
    v_rest = 0
    m = m_inf(v_rest)
    h = h_inf(v_rest)
    n = n_inf(v_rest)
    v[0] = v_rest

    # I_s (current) between times 1000 and 5000, 7000 and 11 000, and 13 000 and 17 000 is
    # changed to 4.4, 25, and 40 respectivley (mA)
    I_s = np.zeros(len(time))
    I_s[1000:5000] = [6]*4000
    I_s[7000:11000] = [6]*4000        #25
    I_s[13000:17000] = [6]*4000        #40

    #g_n (resistance of n)
    g_n = np.zeros(len(time))
    g_n[1000:5000] = [120]*4000     # Normal
    g_n[7000:11000] = [80]*4000     # Low k+ resistance
    g_n[13000:17000] = [200]*4000   # high k+ resistance

    # from i to the time specified, calculate m(t), n(t), and h(t)
    # append to v[i] array
    for i in range(1,len(time)):
        m += (alpha_m(v[i-1])*(1-m) - beta_m(v[i-1])*m)*dt
        n += (alpha_n(v[i-1])*(1-n) - beta_n(v[i-1])*n)*dt
        h += (alpha_h(v[i-1])*(1-h) - beta_h(v[i-1])*h)*dt
        dv = (1./c)*(I_s[i-1] - g_n[i-1]*m**3*h*(v[i-1]-v_n) - g_k*n**4*(v[i-1]-v_k) - g_l*(v[i-1]-v_l))*dt
        v[i] = v[i-1]+dv

    returnOutput(v)

    # plot that shiz
    p = pl.plot(time,v,time,I_s-50)
    pl.xlabel("time (ms)")
    pl.ylabel("potential difference,V (mV)")
    #pl.show()

    return v[13000:17000]

def plot(flouroxitine, test):
    fig = Figure(figsize = (6.5,6.5), dpi = 100)
    plot1 = fig.add_subplot(111)

    #plot1.set_ylim(bottom=0, top=90)

    if test == 1:
        plot1.plot(HHmodel(100))
    else:
        plot1.plot(HHmodel(flouroxitine))

    canvas = FigureCanvasTkAgg(fig, master = window)
    canvas.get_tk_widget().grid(row = 3, column = 2)
    #canvas.draw()

def perCentInhibition(miligrams, icfifty, hillCoeff):
    inhibition = 100 / (1 + (pow((icfifty / miligrams), hillCoeff)))
    return inhibition

window = Tk()
window.title('SSRI effect on HHModel')
window.geometry("665x900")
window.resizable(False, False)

FluoxLabel = Label(master=window, text="[Fluoxitine] (ug): " + str(obj.flouroxMg))
FluoxLabel.place(x = 350, y = 695)
SetLabel = Label(master=window, text="[Sertraline] (ug): " + str(obj.setralineMg))
SetLabel.place(x = 350, y = 760)
citLabel = Label(master=window, text="[Citalopram] (ug): " + str(obj.citalopramMg))
citLabel.place(x = 350, y = 820)

def drugIncreaseFluox(index, iseefifty1, hillcoeff1):
    if index == "zerozeroone":
        obj.flouroxMg = obj.flouroxMg + 0.01
    elif index == "zeroone":
        obj.flouroxMg = obj.flouroxMg + 0.1
    elif index == "one":
        obj.flouroxMg = obj.flouroxMg + 1
    elif index == "ten":
        obj.flouroxMg = obj.flouroxMg + 10
    else:
        obj.flouroxMg = obj.flouroxMg
    plot(perCentInhibition(obj.flouroxMg, iseefifty1, hillcoeff1), 0)
    FluoxLabel.config(text = "[Fluoxitine] (ug): " + str(obj.flouroxMg))

def drugDecreaseFluox(index, iseefifty1, hillcoeff1):
    if index == "zerozeroone":
        obj.flouroxMg = obj.flouroxMg - 0.01
    elif index == "zeroone":
        obj.flouroxMg = obj.flouroxMg - 0.1
    elif index == "one":
        obj.flouroxMg = obj.flouroxMg - 1
    elif index == "ten":
        obj.flouroxMg = obj.flouroxMg - 10
    else:
        obj.flouroxMg = obj.flouroxMg
    plot(perCentInhibition(obj.flouroxMg, iseefifty1, hillcoeff1), 0)
    FluoxLabel.config(text = "[Fluoxitine] (ug): " + str(obj.flouroxMg))

def drugIncreaseSet(index, iseefifty1, hillcoeff1):
    if index == "zerozeroone":
        obj.setralineMg = obj.setralineMg + 0.01
    elif index == "zeroone":
        obj.setralineMg = obj.setralineMg + 0.1
    elif index == "one":
        obj.setralineMg = obj.setralineMg + 1
    elif index == "ten":
        obj.setralineMg = obj.setralineMg + 10
    else:
        obj.setralineMg = obj.setralineMg
    plot(perCentInhibition(obj.setralineMg, iseefifty1, hillcoeff1), 0)
    SetLabel.config(text = "[Sertraline] (ug): " + str(obj.setralineMg))

def drugDecreaseSet(index, iseefifty1, hillcoeff1):
    if index == "zerozeroone":
        obj.setralineMg = obj.setralineMg - 0.01
    elif index == "zeroone":
        obj.setralineMg = obj.setralineMg - 0.1
    elif index == "one":
        obj.setralineMg = obj.setralineMg - 1
    elif index == "ten":
        obj.setralineMg = obj.setralineMg - 10
    else:
        obj.setralineMg = obj.setralineMg
    plot(perCentInhibition(obj.setralineMg, iseefifty1, hillcoeff1), 0)
    SetLabel.config(text = "[Sertraline] (ug): " + str(obj.setralineMg))

def drugIncreaseCit(index, iseefifty1, hillcoeff1):
    if index == "zerozeroone":
        obj.citalopramMg = obj.citalopramMg + 0.01
    elif index == "zeroone":
        obj.citalopramMg = obj.citalopramMg + 0.1
    elif index == "one":
        obj.citalopramMg = obj.citalopramMg + 1
    elif index == "ten":
        obj.citalopramMg = obj.citalopramMg + 10
    else:
        obj.citalopramMg = obj.citalopramMg
    plot(perCentInhibition(obj.citalopramMg, iseefifty1, hillcoeff1), 0)
    citLabel.config(text = "[Citalopram] (ug): " + str(obj.citalopramMg))

def drugDecreaseCit(index, iseefifty1, hillcoeff1):
    if index == "zerozeroone":
        obj.citalopramMg = obj.citalopramMg - 0.01
    elif index == "zeroone":
        obj.citalopramMg = obj.citalopramMg - 0.1
    elif index == "one":
        obj.citalopramMg = obj.citalopramMg - 1
    elif index == "ten":
        obj.citalopramMg = obj.citalopramMg - 10
    else:
        obj.citalopramMg = obj.citalopramMg
    plot(perCentInhibition(obj.citalopramMg, iseefifty1, hillcoeff1), 0)
    citLabel.config(text = "[Citalopram] (ug): " + str(obj.citalopramMg))

def returnOutput(v):
    #NumOfFires = 0
    testArray = v[13000:17000]
    fileOUTPUT = open('output.txt', 'a')

    peaks,y = scipy.signal.find_peaks(testArray, height=30)

    fileOUTPUT.write('\n')
    fileOUTPUT.write(str(round(obj.flouroxMg, 5)))
    fileOUTPUT.write(',')
    fileOUTPUT.write(str(peaks.size))

def testReturn():
    for p in range(400):
        print(p)
        obj.flouroxMg += 0.1
        HHmodel(perCentInhibition(obj.flouroxMg, 27.7, 1.3))


#testReturn()

plot(perCentInhibition(obj.flouroxMg, 6, 1), 0)
buttonFluoxIncrease001 = Button(master = window, text="+0.01", width = 5, command=lambda: drugIncreaseFluox("zerozeroone", 6 , 1)).place(x = 50, y = 675)
buttonFluoxIncrease01 = Button(master = window, text="+0.1", width = 5, command=lambda: drugIncreaseFluox("zeroone",  6 , 1)).place(x = 120, y = 675)
buttonFluoxIncrease1 = Button(master = window, text="+1", width = 5, command=lambda: drugIncreaseFluox("one",  6 , 1)).place(x = 190, y = 675)
buttonFluoxIncrease10 = Button(master = window, text="+10", width = 5, command=lambda: drugIncreaseFluox("ten",  6 , 1)).place(x = 260, y = 675)
buttonFluoxDecrease001 = Button(master = window, text ="-0.01", width = 5, command = lambda: drugDecreaseFluox("zerozeroone", 6 , 1)).place(x = 50, y = 705)
buttonFluoxDecrease01 = Button(master = window, text ="-0.1", width = 5, command = lambda: drugDecreaseFluox("zeroone", 6 , 1)).place(x = 120, y = 705)
buttonFluoxDecrease1 = Button(master = window, text ="-1", width = 5, command = lambda: drugDecreaseFluox("one", 6 , 1)).place(x = 190, y = 705)
buttonFluoxDecrease10 = Button(master = window, text ="-10", width = 5, command = lambda: drugDecreaseFluox("ten", 6 , 1)).place(x = 260, y = 705)

buttonSetIncrease001 = Button(master = window, text="+0.01", width = 5, command=lambda: drugIncreaseSet("zerozeroone", 2.1 , 1.3)).place(x = 50, y = 735)
buttonSetIncrease01 = Button(master = window, text="+0.1", width = 5, command=lambda: drugIncreaseSet("zeroone",  2.1 , 1.3)).place(x = 120, y = 735)
buttonSetIncrease1 = Button(master = window, text="+1", width = 5, command=lambda: drugIncreaseSet("one",  2.1 , 1.3)).place(x = 190, y = 735)
buttonSetIncrease10 = Button(master = window, text="+10", width = 5, command=lambda: drugIncreaseSet("ten",  2.1 , 1.3)).place(x = 260, y = 735)
buttonSetIncrease001 = Button(master = window, text="-0.01", width = 5, command=lambda: drugDecreaseSet("zerozeroone", 2.1 , 1.3)).place(x = 50, y = 765)
buttonSetIncrease01 = Button(master = window, text="-0.1", width = 5, command=lambda: drugDecreaseSet("zeroone",  2.1 , 1.3)).place(x = 120, y = 765)
buttonSetIncrease1 = Button(master = window, text="-1", width = 5, command=lambda: drugDecreaseSet("one",  2.1 , 1.3)).place(x = 190, y = 765)
buttonSetIncrease10 = Button(master = window, text="-10", width = 5, command=lambda: drugDecreaseSet("ten",  2.1 , 1.3)).place(x = 260, y = 765)

buttonCitIncrease001 = Button(master = window, text="+0.01", width = 5, command=lambda: drugIncreaseCit("zerozeroone", 27.7 , 1.3)).place(x = 50, y = 795)
buttonCitIncrease01 = Button(master = window, text="+0.1", width = 5, command=lambda: drugIncreaseCit("zeroone", 27.7 , 1.3)).place(x = 120, y = 795)
buttonCitIncrease1 = Button(master = window, text="+1", width = 5, command=lambda: drugIncreaseCit("one", 27.7 , 1.3)).place(x = 190, y = 795)
buttonCitIncrease10 = Button(master = window, text="+10", width = 5, command=lambda: drugIncreaseCit("ten", 27.7 , 1.3)).place(x = 260, y = 795)
buttonCitDecrease001 = Button(master = window, text="-0.01", width = 5, command=lambda: drugDecreaseCit("zerozeroone", 27.7 , 1.3)).place(x = 50, y = 825)
buttonCitDecrease01 = Button(master = window, text="-0.1", width = 5, command=lambda: drugDecreaseCit("zeroone", 27.7 , 1.3)).place(x = 120, y = 825)
buttonCitDecrease1 = Button(master = window, text="-1", width = 5, command=lambda: drugDecreaseCit("one", 27.7 , 1.3)).place(x = 190, y = 825)
buttonCitDecrease10 = Button(master = window, text="-10", width = 5, command=lambda: drugDecreaseCit("ten",27.7 , 1.3)).place(x = 260, y = 825)

window.mainloop()


# Many thanks to the unknown creator at http://funpythonprojects.blogspot.com/2013/07/hodgkin-huxley-model-part-1.html 
# for the amazingly detailed and easy to understand tutorial