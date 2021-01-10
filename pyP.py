"""
pyP
=======
`pyP` is a Python GUI tool to pick arrival time of P-phase
"""
class Index(object):
    '''
    Switch panel for the waveform traces by clcking the buttons
    '''
    ind = 0
    def nextKey(self, event):
        if event.key == 'down':
            self.ind += 1
            if self.ind > totalnumtrace-1:
                self.ind = 0
            for ax in axs:
                ax.set_visible(False)
                ax.tick_params(labelbottom=False)
            axinit = self.ind
            gsind = 0
            for j in np.arange(axinit, axinit+worknumtrace, 1):
                if j > totalnumtrace-1:
                    j = j - totalnumtrace
                axs[j].set_visible(True)
                axs[j].set_position(gs[gsind].get_position(fig))
                gsind += 1

            axs[j].tick_params(labelbottom=True)
            fig.canvas.draw()

    def nextFurtherKey(self, event):
        if event.key == 'right':
            self.ind += worknumtrace - 1
            if self.ind > totalnumtrace-1:
                self.ind = self.ind - totalnumtrace
            for ax in axs:
                ax.set_visible(False)
                ax.tick_params(labelbottom=False)
            axinit = self.ind
            gsind = 0
            for j in np.arange(axinit, axinit+worknumtrace, 1):
                if j > totalnumtrace-1:
                    j = j - totalnumtrace
                axs[j].set_visible(True)
                axs[j].set_position(gs[gsind].get_position(fig))
                gsind += 1
            axs[j].tick_params(labelbottom=True)
            fig.canvas.draw()

    def prevKey(self, event):
        if event.key == 'up':
            self.ind -= 1
            if self.ind < 0:
                self.ind = totalnumtrace-1
            for ax in axs:
                ax.set_visible(False)
                ax.tick_params(labelbottom=False)
            axinit = self.ind
            gsind = 0
            for j in np.arange(axinit, axinit+worknumtrace, 1):
                if j > totalnumtrace-1:
                    j = j - totalnumtrace
                axs[j].set_visible(True)
                axs[j].set_position(gs[gsind].get_position(fig))
                gsind += 1
            axs[j].tick_params(labelbottom=True)
            fig.canvas.draw()

    def prevFurtherKey(self, event):
        if event.key == 'left':
            self.ind -= worknumtrace - 1
            if self.ind < 0:
                self.ind = self.ind + totalnumtrace
            for ax in axs:
                ax.set_visible(False)
                ax.tick_params(labelbottom=False)
            axinit = self.ind
            gsind = 0
            for j in np.arange(axinit, axinit+worknumtrace, 1):
                if j > totalnumtrace-1:
                    j = j - totalnumtrace
                axs[j].set_visible(True)
                axs[j].set_position(gs[gsind].get_position(fig))
                gsind += 1
            axs[j].tick_params(labelbottom=True)
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
            tmpx = ax.lines[0].get_xdata()
            amarkerlist[axind] = -tmpx[0] # arrival time relative to data starttime
            fig.texts[-2].set_text('Last pick: '+str(stanamelist[axind])+' '+str('{:.2f}'.format(-tmpx[0]))+' s')
            fig.texts[-2].set_color('C'+str(axind))

    elif event.inaxes and event.inaxes == axbutton:
        if event.button == 1: # Left click only
            ax = event.inaxes
            tmptime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            fig.texts[-1].set_text('Saved! '+tmptime)

            #columns = ['amarker', 'Station']
            #df = pd.DataFrame(np.array([amarkerlist, stanamelist]).T, columns = columns)
            #df.to_csv('log.csv')

            tmp = glob.glob(args.sacfiles)[0].split('/')[:-1]
            dirname = '/'.join(tmp)
            if len(dirname) == 0:
                dirname = '.'
            tmplist = np.zeros(len(st))
            for j in range(len(st)):
                trace = st[j].copy()
                tmplist[j] = trace.stats.sac.a
                if trace.stats.sac.a != amarkerlist[j]:
                    filename = glob.glob(dirname+'/'+trace.stats.network+'.'+trace.stats.station+'.'+trace.stats.location+'.'+trace.stats.channel+'*.SAC')
                    if len(filename) > 1:
                        tmptime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                        print('['+tmptime+'] Error!', trace.stats.network+'.'+trace.stats.station+'.'+trace.stats.location+'.'+trace.stats.channel,
                               'is duplicated. Stopped.')
                        print('['+tmptime+'] Error!', trace.stats.network+'.'+trace.stats.station+'.'+trace.stats.location+'.'+trace.stats.channel,
                               'is duplicated. Stopped.', file=logfile)
                        print('['+str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+'] pyP ended, unexpectedly. See above error.', file=logfile)
                        logfile.close()
                        exit(1)

                    abspath = os.path.abspath(filename[0])
                    print('['+tmptime+'] Saved:', '{:.3f}'.format(trace.stats.sac.a), '-->', '{:.3f}'.format(amarkerlist[j]), '|', abspath, file=logfile)
                    print('['+tmptime+'] Saved:', '{:.3f}'.format(trace.stats.sac.a), '-->', '{:.3f}'.format(amarkerlist[j]), '|', abspath)
                    trace.stats.sac.a = amarkerlist[j]
                    st[j] = trace
                    trace.write(abspath, format='SAC')
            #if (tmplist == amarkerlist).all():
            #    print('['+tmptime+'] Save log: Nothing has been changed since the last pick.')

    fig.canvas.draw()

def keypress(event):
    '''
    press ",": zoom in ylim
    press ".": zoom out ylim
    press "a": zoom in xlim
    press "z": zoom out xlim
    press "x": Reset xlim
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

def aziequi(ax, azilist, dellist):
    d, a = np.array(dellist), 90-np.array(azilist)
    x, y=(d*np.cos(a*np.pi/180.0), d*np.sin(a*np.pi/180.0))
    sc=ax.scatter(x, y, s=15, marker='^', edgecolor='none', facecolor='gray', alpha=0.85, zorder=10)
    sc=ax.scatter(x, y, s=15, marker='^', edgecolor='k', facecolor='none', alpha=1, lw=0.5, zorder=10)
    theta=np.linspace(0, 360, 360)
    for i in [30, 60, 90]:
        x, y=(i*np.cos(theta*np.pi/180.0), i*np.sin(theta*np.pi/180.0))
        ax.plot(x, y, color='k', zorder=0, solid_capstyle='round', lw=0.1)
        x, y=(i*np.cos(-90*np.pi/180.0), i*np.sin(-90*np.pi/180.0))
        text = ax.text(x, y, str(i)+'$\degree$', size=8, va='center', ha='center')
        text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()])
    delrange = np.linspace(0, 100, 10)
    for i in np.arange(0, 360, 30):
        x, y=(delrange*np.cos(i*np.pi/180.0), delrange*np.sin(i*np.pi/180.0))
        ax.plot(x, y, color='k', zorder=0, solid_capstyle='round', lw=0.1)

    x, y=(100*np.cos(theta*np.pi/180.0), 100*np.sin(theta*np.pi/180.0))
    ax.plot(x, y, color='k', solid_capstyle='round', lw=1)
    ax.fill(x, y, edgecolor='none', facecolor='w', zorder=0)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_yticks([])
    return sc

def updateColorCycle(cmapkey):
    cmap = plt.get_cmap(cmapkey, 8)
    custom_color_cycle = [ str(mpl.colors.rgb2hex(cmap(i)[:3])) for i in range(cmap.N) ]
    plt.rc('axes', prop_cycle=(cycler(color=custom_color_cycle)))

def checkArgparse():
    # Argument check
    parser = argparse.ArgumentParser(description='pyP: Python based P-arrival picking tool',
                                    usage='pyP.py [-h] displayNum sacfiles elat elon'+\
                                    '\n\nExample:'+\
                                    '\n>>> python pyP.py 7 "../II*.SAC, ../IU.*.SAC" 36.11 140.10\n ')
    parser.add_argument('displayNum', help='Number of traces shown in display (e.g., 7)', type=int)
    parser.add_argument('sacfiles', help='SAC files you want to pick P arrival (e.g., "./*.SAC"). comma-separated list is available. *Do not forget quotation marks!', type=str)
    parser.add_argument('elat', help='Latitude of epicentre (for station azimuth)', type=float)
    parser.add_argument('elon', help='Longitude of epicentre (for station azimuth)', type=float)
    args = parser.parse_args()
    return args

def loadStream(args):
    tmplist = [x.strip() for x in args.sacfiles.split(',')]
    st = read(tmplist[0])
    for sacfiles in tmplist[1:]:
        st += read(sacfiles)
    elat, elon = float(args.elat), float(args.elon) #38.392, 39.085
    for j in range(len(st)):
        trace = st[j].copy()
        tmp = geod.Inverse(elat, elon, trace.stats.sac.stla, trace.stats.sac.stlo)
        azi = tmp['azi1']
        if azi < 0:
            azi = 360 + azi
        trace.stats.sac.user0 = azi
        trace.stats.sac.user1 = tmp['a12']
        try:
            tmpa = trace.stats.sac.a
        except AttributeError:
            trace.stats.sac.a = 300.0
        st[j] = trace
    st.traces.sort(key=lambda x: x.stats.sac.user0) # sort by azimuth
    return st

def outputlog(st, args):
    for j in range(len(st)):
        trace = st[j].copy()
        filename = trace.stats.network+'.'+trace.stats.station+'.'+trace.stats.location+'.'+trace.stats.channel
        tmptime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print('['+tmptime+']', str(j+1).rjust(3)+'/'+str(totalnumtrace), '|', '{:.3f}'.format(trace.stats.sac.a), '|', filename, file=logfile)




import os
import sys
import argparse
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from obspy import read
import matplotlib as mpl
from cycler import cycler
updateColorCycle('Set2')
from datetime import datetime
from matplotlib.gridspec import GridSpec
from geographiclib.geodesic import Geodesic
geod = Geodesic.WGS84
from matplotlib.widgets import Button, MultiCursor
import matplotlib.patheffects as path_effects

# Check parsed arguments, and load data
args = checkArgparse()
logfile = open('pyP.log', 'a')
print('['+str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+'] pyP started '+os.getcwd(), file=logfile)
st = loadStream(args)

# Base subplots (not shown in display)
totalnumtrace = len(st)
fig, axs = plt.subplots(totalnumtrace, figsize=(15, 10))
fig.subplots_adjust(left=0.3)
fig.canvas.set_window_title('pyP')

textpe = [path_effects.Stroke(linewidth=2, foreground='w', alpha=1), path_effects.Normal()]
xmin0, xmax0 = -100, 400 # default xlim
amarkerlist = np.zeros(totalnumtrace)
stanamelist, azilist, dellist = [], [], []
for j,ax in enumerate(axs.flat):
    trace = st[j].copy()
    amarker = trace.stats.sac.a
    amarkerlist[j] =  amarker
    stanamelist.append(trace.stats.station)
    azilist.append(trace.stats.sac.user0)
    dellist.append(trace.stats.sac.user1)
    df = trace.stats.sampling_rate
    x = np.arange(0, trace.stats.npts, 1) / df - amarker
    y = trace.data - np.mean(trace.data[ int((amarker+xmin0)*df):int((amarker)*df) ])
    ax.plot(x,y, lw=1, color='k')
    ax.axvline(0, lw=0.75, color='k')

    ax.set_xlim(xmin0, xmax0)
    ymin0 = -max( np.abs(y[ int((amarker+xmin0)*df):int((amarker+xmax0)*df) ]) ) * 1.2
    ymax0 = max( np.abs(y[ int((amarker+xmin0)*df):int((amarker+xmax0)*df) ]) ) * 1.2
    ax.set_ylim(ymin0, ymax0)
    textlabel0 = str(j+1)+'/'+str(totalnumtrace)+'\n'+str(trace.stats.network)+'.'+str(trace.stats.station)+'.'+str(trace.stats.location)+'.'+str(trace.stats.channel)
    textlabel1 = '\nAzi: '+str('{:.2f}'.format(trace.stats.sac.user0)) + '\nDel: ' + str('{:.2f}'.format(trace.stats.sac.user1))
    text = ax.text(xmin0+(xmax0-xmin0)*0.005, ymax0-(ymax0-ymin0)*0.1, textlabel0+textlabel1, fontsize=8, ha='left', va='top').set_path_effects(textpe)

    ax.patch.set_facecolor('C'+str(j))
    ax.patch.set_alpha(0.1)
    ax.tick_params(labelleft=False)
    ax.tick_params(labelbottom=False)

# Display interactive canvus, define click actions
worknumtrace = args.displayNum
gs = GridSpec(worknumtrace,1,hspace=0)
callback = Index()
kpnext = fig.canvas.mpl_connect('key_press_event', callback.nextKey)
kpnextFurther = fig.canvas.mpl_connect('key_press_event', callback.nextFurtherKey)
kpprev = fig.canvas.mpl_connect('key_press_event', callback.prevKey)
kpprevFurther = fig.canvas.mpl_connect('key_press_event', callback.prevFurtherKey)

curser = MultiCursor(fig.canvas, axs, color='#B12763', lw=1, alpha=0.75)
fig.canvas.mpl_connect('button_press_event', oncpick)
fig.canvas.mpl_connect('key_press_event', keypress)

# Initial setting: hide base subplots and display first `worknumtrace` subplots
for ax in axs:
    ax.set_visible(False)

for j in np.arange(0, worknumtrace, 1):
    axs[j].set_visible(True)
    axs[j].set_position(gs[j].get_position(fig))
axs[j].tick_params(labelbottom=True)

# Usage at top
axp = axs[0].get_position()
usagetext = 'Usage: press key\n[A] Zoom-in xlim      [Z] Zoom-out xlim      [X] Reset xlim\n[.] Zoom-in ylim        [,] Zoom-out ylim\n'+\
r'[$\leftarrow$] Prev more      [$\uparrow$] Prev      [$\downarrow$] Next      [$\rightarrow$] Next more'
fig.text(axp.x0, axp.y1+0.01, usagetext, size=8)

# Station map
axp = axs[-1].get_position()
axstamap = fig.add_axes([axp.x0-0.255, axp.y0, 0.25, 0.25])
axstamap.set_aspect(1)
sc = aziequi(axstamap, azilist, dellist)

# Save function (save amarker list in log.csv) and some information
axbutton = fig.add_axes([axp.x1-0.1, axp.y0-0.075, 0.1, 0.035])
axpb = axbutton.get_position()
fig.text(axp.x0+(axp.width/2), axpb.y0+axpb.height/2, 'Time (s)', ha='center', va='center')
fig.text(axp.x0, axpb.y0+axpb.height/2, 'Last pick: None', ha='left', va='center')
fig.text(axpb.x0-0.01, axpb.y0+axpb.height/2, '', ha='right', va='center')
bsave = Button(axbutton, 'Save')

plt.show()
outputlog(st, args)
print('['+str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+'] pyP ended', file=logfile)
print('['+str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+'] =========', file=logfile)
logfile.close()
