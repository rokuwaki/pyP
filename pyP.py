"""
pyP
=======
`pyP` is a Python GUI tool to pick arrival time of P-phase
Basis usage: execute `pyP.py` where the SAC data is stored
>>> python pyP.py 10
"""

class Index(object):
    '''
    Switch panel for the waveform traces by clcking the buttons
    '''
    ind = 0
    def next(self, event):
        self.ind += 1
        if self.ind > totalnumtrace-1:
            self.ind = 0
        for ax in axs:
            ax.set_visible(False)
        axinit = self.ind
        gsind = 0
        for j in np.arange(axinit, axinit+worknumtrace, 1):
            if j > totalnumtrace-1:
                j = j - totalnumtrace
            axs[j].set_visible(True)
            axs[j].set_position(gs[gsind].get_position(fig))
            gsind += 1
        fig.canvas.draw()

    def nextFurther(self, event):
        self.ind += worknumtrace
        if self.ind > totalnumtrace-1:
            self.ind = self.ind - totalnumtrace
        for ax in axs:
            ax.set_visible(False)
        axinit = self.ind
        gsind = 0
        for j in np.arange(axinit, axinit+worknumtrace, 1):
            if j > totalnumtrace-1:
                j = j - totalnumtrace
            axs[j].set_visible(True)
            axs[j].set_position(gs[gsind].get_position(fig))
            gsind += 1
        fig.canvas.draw()

    def prev(self, event):
        self.ind -= 1
        if self.ind < 0:
            self.ind = totalnumtrace-1
        for ax in axs:
            ax.set_visible(False)
        axinit = self.ind
        gsind = 0
        for j in np.arange(axinit, axinit+worknumtrace, 1):
            if j > totalnumtrace-1:
                j = j - totalnumtrace
            axs[j].set_visible(True)
            axs[j].set_position(gs[gsind].get_position(fig))
            gsind += 1
        fig.canvas.draw()

    def prevFurther(self, event):
        self.ind -= worknumtrace
        if self.ind < 0:
            self.ind = self.ind + totalnumtrace
        for ax in axs:
            ax.set_visible(False)
        axinit = self.ind
        gsind = 0
        for j in np.arange(axinit, axinit+worknumtrace, 1):
            if j > totalnumtrace-1:
                j = j - totalnumtrace
            axs[j].set_visible(True)
            axs[j].set_position(gs[gsind].get_position(fig))
            gsind += 1
        fig.canvas.draw()

    def nextKey(self, event):
        if event.key == 'down':
            self.ind += 1
            if self.ind > totalnumtrace-1:
                self.ind = 0
            for ax in axs:
                ax.set_visible(False)
            axinit = self.ind
            gsind = 0
            for j in np.arange(axinit, axinit+worknumtrace, 1):
                if j > totalnumtrace-1:
                    j = j - totalnumtrace
                axs[j].set_visible(True)
                axs[j].set_position(gs[gsind].get_position(fig))
                gsind += 1
            fig.canvas.draw()

    def nextFurtherKey(self, event):
        if event.key == 'right':
            self.ind += worknumtrace
            if self.ind > totalnumtrace-1:
                self.ind = self.ind - totalnumtrace
            for ax in axs:
                ax.set_visible(False)
            axinit = self.ind
            gsind = 0
            for j in np.arange(axinit, axinit+worknumtrace, 1):
                if j > totalnumtrace-1:
                    j = j - totalnumtrace
                axs[j].set_visible(True)
                axs[j].set_position(gs[gsind].get_position(fig))
                gsind += 1
            fig.canvas.draw()

    def prevKey(self, event):
        if event.key == 'up':
            self.ind -= 1
            if self.ind < 0:
                self.ind = totalnumtrace-1
            for ax in axs:
                ax.set_visible(False)
            axinit = self.ind
            gsind = 0
            for j in np.arange(axinit, axinit+worknumtrace, 1):
                if j > totalnumtrace-1:
                    j = j - totalnumtrace
                axs[j].set_visible(True)
                axs[j].set_position(gs[gsind].get_position(fig))
                gsind += 1
            fig.canvas.draw()

    def prevFurtherKey(self, event):
        if event.key == 'left':
            self.ind -= worknumtrace
            if self.ind < 0:
                self.ind = self.ind + totalnumtrace
            for ax in axs:
                ax.set_visible(False)
            axinit = self.ind
            gsind = 0
            for j in np.arange(axinit, axinit+worknumtrace, 1):
                if j > totalnumtrace-1:
                    j = j - totalnumtrace
                axs[j].set_visible(True)
                axs[j].set_position(gs[gsind].get_position(fig))
                gsind += 1
            fig.canvas.draw()

def oncpick(event):
    '''
    Left click: Pick arrival time
    Click "Save" button: Keep picked time as `amarker` and save to csv file
    '''
    if event.inaxes and event.inaxes in axs:
        if event.button == 1: # Left click only
            ax = event.inaxes
            x = ax.lines[0].get_xdata()
            ax.lines[0].set_xdata(x-event.xdata)
            axind = fig.axes.index(ax)
            amarkerlist[axind] = ax.get_xlim()[0]
    elif event.inaxes and event.inaxes == axbutton:
        if event.button == 1: # Left click only
            ax = event.inaxes
            print(amarkerlist)
            tmp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            fig.texts[-1].set_text('Saved! '+tmp)

    fig.canvas.draw()

def keypress(event):
    '''
    press ",": Enlarge x range
    press ".": Shrink x range
    press "a": Enlarge y range
    press "z": Shrink y range
    press "x": Reset xrange
    '''
    for ax in axs:
        x = ax.lines[0].get_xdata()
        y = ax.lines[0].get_ydata()
        xmin = ax.get_xlim()[0]
        xmax = ax.get_xlim()[1]
        ymin = ax.get_ylim()[0]
        ymax = ax.get_ylim()[1]
        if event.key == '.':
            ax.set_ylim(ymin*0.8, ymax*0.8)
        if event.key == ',':
            ax.set_ylim(ymin*1.2, ymax*1.2)
        if event.key == 'a':
            xmin = xmin * 0.8
            xmax = xmax * 0.8
            if xmax-xmin < 1:
                exit
            else:
                ax.set_xlim(xmin, xmax)
        if event.key == 'z':
            ax.set_xlim(xmin*1.2, xmax*1.2)
        if event.key == 'x':
            ax.set_xlim(xmin0, xmax0)
        xmin = ax.get_xlim()[0]
        xmax = ax.get_xlim()[1]
        ymin = ax.get_ylim()[0]
        ymax = ax.get_ylim()[1]
        ax.texts[0].set_x(xmin+(xmax-xmin)*0.005)
        ax.texts[0].set_y(ymax-(ymax-ymin)*0.1)
    fig.canvas.draw()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, MultiCursor
from matplotlib.gridspec import GridSpec
import sys
import glob
import argparse
from obspy import read
import mplcursors
import matplotlib.patheffects as path_effects
from datetime import datetime
import matplotlib as mpl
from cycler import cycler
cmap = plt.get_cmap('Set3', 12)
custom_color_cycle = [ str(mpl.colors.rgb2hex(cmap(i)[:3])) for i in range(cmap.N) ]
plt.rc('axes', prop_cycle=(cycler(color=custom_color_cycle)))

textpe = [path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()]
xmin0 = -100
xmax0 = 400

parser = argparse.ArgumentParser(description='pyP: Python based P-arrival picking tool')
parser.add_argument('arg1', help='Number of traces shown in display (e.g., 7)')
parser.add_argument('arg2', help='SAC files you want to pick P arrival (e.g., "./*.SAC") *Do not forget quotation marks!')
args = parser.parse_args()
try:
    worknumtrace = int(args.arg1)
except ValueError:
    worknumtrace = 7
    sacfiles = args.arg1
    print('')
    print('  Caution: You did not input the number of traces shown at once.')
    print('           --> I set 7 as the number of traces on display at once: `worknumtrace`.')
    print('  Next time, input the number of traces you like as:')
    print('  > python pyP.py 7 "./*.SAC"')
    print('')
sacfiles = args.arg2
if len(glob.glob(sacfiles)) == 0:
    print('')
    print('  Error: No SAC files identified. Pleae specify the correct SAC files. Stopped.')
    print('')
    sys.exit(1)


st = read(sacfiles)
totalnumtrace = len(glob.glob(sacfiles))
#worknumtrace = 8
amarkerlist = np.zeros(totalnumtrace)

fig, axs = plt.subplots(totalnumtrace)

for j,ax in enumerate(axs.flat):
    trace = st[j].copy()
    amarker = trace.stats.sac.a
    amarkerlist[j] =  amarker
    df = trace.stats.sampling_rate
    x = np.arange(0, trace.stats.npts, 1) / df - amarker
    y = trace.data - np.mean(trace.data[ int((amarker+xmin0)*df):int((amarker)*df) ])
    ax.plot(x,y, lw=1, color='k')
    ax.axvline(0, lw=0.75, color='k')

    ax.set_xlim(xmin0, xmax0)
    ymin0 = -max( np.abs(y[ int((amarker+xmin0)*df):int((amarker+xmax0)*df) ]) ) * 1.2
    ymax0 = max( np.abs(y[ int((amarker+xmin0)*df):int((amarker+xmax0)*df) ]) ) * 1.2
    ax.set_ylim(ymin0, ymax0)
    textlabel = str(j+1)+'/'+str(totalnumtrace)+'\n'+str(trace.stats.network)+'.'+str(trace.stats.station)+'.'+str(trace.stats.location)+'.'+str(trace.stats.channel)
    #axp = ax.get_position()
    text = ax.text(xmin0+(xmax0-xmin0)*0.005, ymax0-(ymax0-ymin0)*0.1, textlabel, fontsize=8, ha='left', va='top').set_path_effects(textpe)

    ax.patch.set_facecolor('C'+str(j))
    ax.patch.set_alpha(0.2)
    ax.tick_params(labelleft=False)
    ax.tick_params(labelbottom=False)

gs = GridSpec(worknumtrace,1, hspace=0)



callback = Index()

axp = axs[-1].get_position()

#axp.x1-axp.width*0.1, axp.y1+0.01, axp.width*0.1, axp.height*0.75
tmpw = axp.width*0.08
axprevFurther = plt.axes([axp.x0, axp.y0-0.05, tmpw, 0.035])
axprev = plt.axes([axp.x0+tmpw+0.02, axp.y0-0.05, tmpw, 0.035])
axnext = plt.axes([axp.x0+tmpw*2+0.02*2, axp.y0-0.05, tmpw, 0.035])
axnextFurther = plt.axes([axp.x0+tmpw*3+0.02*3, axp.y0-0.05, tmpw, 0.035])

bprevFurther = Button(axprevFurther, '<<')
bprevFurther.on_clicked(callback.prevFurther)

bprev = Button(axprev, '<')
bprev.on_clicked(callback.prev)

bnext = Button(axnext, '>')
bnext.on_clicked(callback.next)

bnextFurther = Button(axnextFurther, '>>')
bnextFurther.on_clicked(callback.nextFurther)

kpnext = fig.canvas.mpl_connect('key_press_event', callback.nextKey)
kpnextFurther = fig.canvas.mpl_connect('key_press_event', callback.nextFurtherKey)
kpprev = fig.canvas.mpl_connect('key_press_event', callback.prevKey)
kpprevFurther = fig.canvas.mpl_connect('key_press_event', callback.prevFurtherKey)

curser = MultiCursor(fig.canvas, axs, color='#B12763', lw=1, alpha=0.75)

fig.canvas.mpl_connect('button_press_event', oncpick)
fig.canvas.mpl_connect('key_press_event', keypress)


for ax in axs:
    ax.set_visible(False)

for j in np.arange(0, worknumtrace, 1):
    axs[j].set_visible(True)
    axs[j].set_position(gs[j].get_position(fig))

axp = axs[-1].get_position()
axbutton = fig.add_axes([axp.x1-0.1, axp.y0-0.05, 0.1, 0.035])

axp = axbutton.get_position()
fig.text(axp.x0-0.01, axp.y0+axp.height/2, '', ha='right', va='center')
bsave = Button(axbutton, 'Save')

plt.show()
