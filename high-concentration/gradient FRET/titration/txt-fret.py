import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import PhotonSeriesCalculate as psc
from matplotlib.widgets import RadioButtons
from matplotlib.ticker import FuncFormatter
plt.style.use('ggplot')

data = '0.7uM/'
rep = 1
save_pic = 0
vary = 0

xhist = pd.read_table(data+'x_hist.txt', names=['A']).A.values
yhist = pd.read_table(data+'bootstrap_'+str(rep)+'.txt', names=['A']).A.values

fret = psc.FretDist(xhist, yhist)
fig, ax = plt.subplots()
fret.showFret(ax)
fig.canvas.mpl_connect('button_press_event', fret.fret_click)
plt.show()

fig2, ax2 = plt.subplots()
fret.showFret(ax2)
fret.fret_fit(ax2, vary=vary)
print fret.fit_notes
plt.show()


if save_pic:
    fig2.savefig(data+'fret_fit_'+str(rep)+'.png', dpi=300)
    with open(data+'fit_note_'+str(rep)+'.txt', 'w') as f:
        f.write(fret.fit_notes)
