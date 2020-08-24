import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import curve_fit
from scipy.stats import poisson
from scipy.stats import expon
from scipy.stats import chisquare
from matplotlib.widgets import Button
from lmfit import Minimizer, Parameters, report_fit

def bintime(trace, dt, dtBin):
    if dt == dtBin:
        return trace, trace
    binNum = int(dtBin/dt)
    trace = trace[:trace.size/binNum*binNum]
    return trace, np.sum(trace.reshape((binNum, -1), order = 'F'), axis = 0)

def readSdt(fileName):
    fid = open(fileName,'rb')

    #File Header
    revision = np.fromfile(fid,'i2',1)
    info_offs = np.fromfile(fid,'i4',1)
    info_length = np.fromfile(fid,'i2',1)
    setup_offs = np.fromfile(fid,'i4',1)
    setup_length = np.fromfile(fid,'i2',1)
    data_block_offs = np.fromfile(fid,'i4',1)
    no_of_data_blocks = np.fromfile(fid,'i2',1)
    data_block_length = np.fromfile(fid,'i4',1)
    meas_desc_block_offs = np.fromfile(fid,'i4',1)
    no_of_meas_desc_blocks = np.fromfile(fid,'i2',1)
    meas_desc_block_length = np.fromfile(fid,'i2',1)

    #Meas Info
    fid.seek(meas_desc_block_offs,0)
    Meas_time = np.fromfile(fid,'S9',1)
    Meas_date = np.fromfile(fid,'S11',1)

    #Find the Data Block A
    fid.seek(data_block_offs,0)
    block_A_No = np.fromfile(fid,'i2',1)
    data_A_offset = np.fromfile(fid,'i4',1)
    next_block_A_offset = np.fromfile(fid,'i4',1)
    block_type_A = np.fromfile(fid,'i2',1)
    meas_desc_block_no_A = np.fromfile(fid,'i2',1)
    lblock_no_A = np.fromfile(fid,'i4',1)
    block_length_A = np.fromfile(fid,'i4',1)

    #Find the Data Block B
    fid.seek(next_block_A_offset,0)
    block_B_No = np.fromfile(fid,'i2',1)
    data_B_offset = np.fromfile(fid,'i4',1)
    next_block_B_offset = np.fromfile(fid,'i4',1)
    block_type_B = np.fromfile(fid,'i2',1)
    meas_desc_block_no_B = np.fromfile(fid,'i2',1)
    lblock_no_B = np.fromfile(fid,'i4',1)
    block_length_B = np.fromfile(fid,'i4',1)

    #Find the Data Block C
    fid.seek(next_block_B_offset,0)
    block_C_No = np.fromfile(fid,'i2',1)
    data_C_offset = np.fromfile(fid,'i4',1)
    next_block_C_offset = np.fromfile(fid,'i4',1)
    block_type_C = np.fromfile(fid,'i2',1)
    meas_desc_block_no_C = np.fromfile(fid,'i2',1)
    lblock_no_C = np.fromfile(fid,'i4',1)
    block_length_C = np.fromfile(fid,'i4',1)

    #Find the Data Block D
    fid.seek(next_block_C_offset,0)
    block_D_No = np.fromfile(fid,'i2',1)
    data_D_offset = np.fromfile(fid,'i4',1)
    next_block_D_offset = np.fromfile(fid,'i4',1)
    block_type_D = np.fromfile(fid,'i2',1)
    meas_desc_block_no_D = np.fromfile(fid,'i2',1)
    lblock_no_D = np.fromfile(fid,'i4',1)
    block_length_D = np.fromfile(fid,'i4',1)

    #Get trajectories
    fid.seek(data_A_offset,0)
    pointA = block_length_A/4
    if pointA==0:
        pointA = 65536
    TraceA = np.fromfile(fid,'i4',pointA) # Donor_s

    fid.seek(data_B_offset,0)
    pointB = block_length_B/4
    if pointB==0:
        pointB = 65536
    TraceB = np.fromfile(fid,'i4',pointB) # Donor_p

    fid.seek(data_C_offset,0)
    pointC = block_length_C/4
    if pointC==0:
        pointC = 65536
    TraceC = np.fromfile(fid,'i4',pointC) # Acceptor_s

    fid.seek(data_D_offset,0)
    pointD = block_length_D/4
    if pointD==0:
        pointD = 65536
    TraceD = np.fromfile(fid,'i4',pointD) # Acceptor_p

    points = min(pointA,pointB,pointC,pointD)
    TraceA = TraceA[:points]
    TraceB = TraceB[:points]
    TraceC = TraceC[:points]
    TraceD = TraceD[:points]

    return TraceA,TraceB,TraceC,TraceD #Donor_s,Donor_p,Acceptor_s,Acceptor_p

    fid.close()

def boundaryPosition(index):
    index = 1*index
    boundary = index[1:] - index[:-1]
    boundary = np.append(index[0], boundary)
    boundary = np.append(boundary, -index[-1])
    #np.where always return a tuple even in a 1D array
    leftSeq = np.where(boundary == 1)[0]
    rightSeq = np.where(boundary == -1)[0]
    return leftSeq, rightSeq

def pois(x, lamb, offset):
    return poisson.pmf(x, lamb)+offset

def one_gau(params, x, y):
    h = params['h']
    c = params['c']
    w = params['w']
    #o = params['o']

    return h*np.exp(-(x-c)**2/(2*w**2))-y#+o

def two_gau(params, x, y):
    h1 = params['h1']
    c1 = params['c1']
    w1 = params['w1']

    h2 = params['h2']
    c2 = params['c2']
    w2 = params['w2']
    
    o  = params['o']

    return h1*np.exp(-(x-c1)**2/(2*w1**2))+h2*np.exp(-(x-c2)**2/(2*w2**2))-y+o

def three_gau(params, x, y):
    h1 = params['h1']
    c1 = params['c1']
    w1 = params['w1']

    h2 = params['h2']
    c2 = params['c2']
    w2 = params['w2']

    h3 = params['h3']
    c3 = params['c3']
    w3 = params['w3']

    o  = params['o']

    return h1*np.exp(-(x-c1)**2/(2*w1**2))+h2*np.exp(-(x-c2)**2/(2*w2**2))+h3*np.exp(-(x-c3)**2/(2*w3**2))-y+o

def subBurstChunk(Dr, Ac, begInd, endInd, chunkPhotons, dt=10, dtBin=1000):
    Tr = Dr + Ac
    subAc = np.array([])
    start = begInd*dtBin/dt #index correct since bintime: dt, dtBin
    stop = start+1
    while (stop <= (endInd*dtBin/dt+1)): #index correct since bintime: dt, dtBin
        if Tr[start:stop].sum() < chunkPhotons:
            stop = stop+1
        elif Tr[start:stop].sum() > chunkPhotons:
            start = start+1
        else:
            subAc = np.append(subAc, Ac[start:stop].sum())
            start = stop
            stop = start+1
    return np.std(subAc)/chunkPhotons

def natural_sort(l):
    prog = re.compile(r'(\d+)(\.txt$)')
    alphanum_key = lambda key: int(re.search(prog, key).group(1))
    return sorted(l, key = alphanum_key)

#correction
class Correct:
    def __init__(self, chl_a_bin, chl_b_bin, dtBin, eta, cross):
        self.chl_a_bin = chl_a_bin
        self.chl_b_bin = chl_b_bin
        self.chl_t_bin = chl_a_bin+chl_b_bin
        self.traceNum = self.chl_a_bin.size
        self.dtBin = dtBin
        self.eta = eta # Anisotropy setup
        self.cross= cross # FRET setup
        self.time = (np.arange(self.traceNum)+1.)*self.dtBin/1e6 # end time point, s, float format in case of overflow

    def fitPoisson(self, trace, cut_offset):
        hist, bins = np.histogram(trace, bins=np.arange(15), density=True)
        bins = bins[:-1]
        cut = bins[hist.cumsum()<.85][-1]+cut_offset+1

        #partial fitting
        bins_fit = bins[:cut]
        hist_fit = hist[:cut]
        parameter, cov_matrix = curve_fit(pois, bins_fit, hist_fit)
        chi2, p_val = chisquare(pois(bins_fit, *parameter), hist_fit)
        return parameter, bins, hist, bins_fit, hist_fit, chi2

    def getNoise(self, cut_offset=0):
        noise_bin_fold = 1

        noise_temp_A,noise_temp_BinA = bintime(self.chl_a_bin, 1, noise_bin_fold)
        noise_temp_B,noise_temp_BinB = bintime(self.chl_b_bin, 1, noise_bin_fold)

        self.a_noise_par, self.a_bins, self.a_hist, self.a_bins_fit, self.a_hist_fit, self.a_chi2 = self.fitPoisson(noise_temp_BinA, cut_offset=0) #bintime = 2*dtBin for better noise fitting
        self.b_noise_par, self.b_bins, self.b_hist, self.b_bins_fit, self.b_hist_fit, self.b_chi2 = self.fitPoisson(noise_temp_BinB, cut_offset=0)
        # self.t_noise_par, self.t_bins, self.t_hist, self.t_bins_fit, self.t_hist_fit, self.t_p_val = self.fitPoisson(bintime(self.chl_t_bin, 1, noise_bin_fold), cut_offset=0)
        self.a_noise_rate = self.a_noise_par[0]/noise_bin_fold/self.dtBin*1e3 #ms-1
        self.b_noise_rate = self.b_noise_par[0]/noise_bin_fold/self.dtBin*1e3 #ms-1
        self.t_noise_rate = self.a_noise_rate+self.b_noise_rate
        # self.t_noise_rate = self.t_noise_par[0]/noise_bin_fold/self.dtBin*1e3 #ms-1

    def showNoise_a(self, ax):
        noise_info = \
            'noise: '+ \
            str(round(self.a_noise_rate, 2))+' KHz' \
            '\noffset: '+ \
            str(round(self.a_noise_par[1], 2))+ \
            '\nchi2: '+ \
            str(round(self.a_chi2, 4))
    
        ax.plot(self.a_bins, self.a_hist, 'o', color='grey')
        ax.plot(self.a_bins_fit, self.a_hist_fit, 'ro')
        ax.plot(self.a_bins, pois(self.a_bins, *self.a_noise_par), '-')
        ax.annotate(noise_info, (10, .1))
        ax.annotate(str(int(self.a_hist_fit.sum()*100))+'%', (self.a_bins_fit[-1], self.a_hist_fit[-1]*1.1))
        ax.set_xlabel('A channel trace distribuion')
        ax.set_ylabel('frequency density')
        ax.set_title('photon counts distribution')

    def showNoise_b(self, ax):
        noise_info = \
            'noise: '+ \
            str(round(self.b_noise_rate, 2))+' KHz' \
            '\noffset: '+ \
            str(round(self.b_noise_par[1], 2))+ \
            '\nchi2: '+ \
            str(round(self.b_chi2, 4))
    
        ax.plot(self.b_bins, self.b_hist, 'o', color='grey')
        ax.plot(self.b_bins_fit, self.b_hist_fit, 'ro')
        ax.plot(self.b_bins, pois(self.b_bins, *self.b_noise_par), '-')
        ax.annotate(noise_info, (10, .1))
        ax.annotate(str(int(self.b_hist_fit.sum()*100))+'%', (self.b_bins_fit[-1], self.b_hist_fit[-1]*1.1))
        ax.set_xlabel('B channel trace distribuion')
        ax.set_ylabel('frequency density')
        ax.set_title('photon counts distribution')

    def getTrace(self):
        self.chl_a_align = (self.chl_a_bin-self.a_noise_rate*self.dtBin/1e3)
        self.chl_b_align = (self.chl_b_bin-self.b_noise_rate*self.dtBin/1e3)*self.eta-self.chl_a_align*self.cross
        self.chl_t_align = self.chl_a_align+self.chl_b_align

    def showTrace(self, ax,  beg=0, end=-1):
        beg = int(beg*1e6/self.dtBin)
        if end == -1:
            end = self.traceNum
        else:
            end = int(end*1e6/self.dtBin)

        ax.plot(self.time[beg:end],  self.chl_a_bin[beg:end], alpha=.7, color='c')
        ax.plot(self.time[beg:end], -self.chl_b_bin[beg:end], alpha=.7, color='m')
        ax.set_xlabel('time / s')
        ax.set_ylabel('corrected fluorescence counts')

#calculation
class Burst(Correct):
    'FRET calculation using summed-points in a burst'
    def __init__(self, chl_a_bin, chl_b_bin, dtBin, eta, cross):
        Correct.__init__(self, chl_a_bin, chl_b_bin, dtBin, eta, cross)

    def findLimPosition(self):
        #assign self.peak_bool and self.back_bool outside
        #find boundary of peakIndex & backIndex
        leftPeakSeq, rightPeakSeq = boundaryPosition(self.peak_bool)
        leftBackSeq, rightBackSeq = boundaryPosition(self.back_bool)

        #unique position
        leftLimSeq = map(lambda x: leftBackSeq[np.where(leftBackSeq <= x)[0][-1]], leftPeakSeq)
        rightLimSeq = map(lambda x: rightBackSeq[np.where(rightBackSeq >= x)[0][0]], rightPeakSeq)
        leftLimSeq = np.unique(np.array(leftLimSeq))
        rightLimSeq = np.unique(np.array(rightLimSeq))

        #remove potential incomplete bursts
        leftLimSeq = np.delete(leftLimSeq, 0)
        rightLimSeq = np.delete(rightLimSeq, 0)
        leftLimSeq = np.delete(leftLimSeq, -1)
        rightLimSeq = np.delete(rightLimSeq, -1)

        self.leftLimSeq = leftLimSeq
        self.rightLimSeq = rightLimSeq
        self.LimSeqNum = self.leftLimSeq.size

    def findSelectIndex(self):
        self.selectIndBin = np.zeros(len(self.chl_t_align), dtype=bool)
        for i in np.arange(self.LimSeqNum):
            self.selectIndBin[self.leftLimSeq[i]:self.rightLimSeq[i]] = True

    def sumPeak(self):
        self.a_peak_sum = np.zeros(self.LimSeqNum)
        self.b_peak_sum = np.zeros(self.LimSeqNum)
        self.t_peak_sum = np.zeros(self.LimSeqNum)
        for i in np.arange(self.LimSeqNum):
            l, r = self.leftLimSeq[i], self.rightLimSeq[i]
            self.a_peak_sum[i] = np.sum(self.chl_a_align[l:r+1])
            self.b_peak_sum[i] = np.sum(self.chl_b_align[l:r+1])
            self.t_peak_sum[i] = np.sum(self.chl_t_align[l:r+1])

    def maxPeak(self):
        self.a_peak_max = np.zeros(self.LimSeqNum)
        self.b_peak_max = np.zeros(self.LimSeqNum)
        self.t_peak_max = np.zeros(self.LimSeqNum)
        for i in np.arange(self.LimSeqNum):
            l, r = self.leftLimSeq[i], self.rightLimSeq[i]
            self.a_peak_max[i] = np.max(self.chl_a_align[l:r+1])
            self.b_peak_max[i] = np.max(self.chl_b_align[l:r+1])
            self.t_peak_max[i] = np.max(self.chl_t_align[l:r+1])

    def sumPoint(self):
        self.a_point_sum = np.ones(self.traceNum)
        self.b_point_sum = np.ones(self.traceNum)
        self.t_point_sum = np.ones(self.traceNum)
        for i in np.arange(self.LimSeqNum):
            l, r = self.leftLimSeq[i], self.rightLimSeq[i]
            self.a_point_sum[l:r+1] = np.repeat(np.sum(self.chl_a_align[l:r+1]), r-l+1)
            self.b_point_sum[l:r+1] = np.repeat(np.sum(self.chl_b_align[l:r+1]), r-l+1)
            self.t_point_sum[l:r+1] = np.repeat(np.sum(self.chl_t_align[l:r+1]), r-l+1)

    def maxPoint(self):
        self.a_point_max = np.ones(self.traceNum)
        self.b_point_max = np.ones(self.traceNum)
        self.t_point_max = np.ones(self.traceNum)
        for i in np.arange(self.LimSeqNum):
            l, r = self.leftLimSeq[i], self.rightLimSeq[i]
            self.a_point_max[l:r+1] = np.repeat(np.max(self.chl_a_align[l:r+1]), r-l+1)
            self.b_point_max[l:r+1] = np.repeat(np.max(self.chl_b_align[l:r+1]), r-l+1)
            self.t_point_max[l:r+1] = np.repeat(np.max(self.chl_t_align[l:r+1]), r-l+1)

    def getBurst(self):
        self.findLimPosition()
        self.findSelectIndex()    
        self.sumPeak()
        self.maxPeak()
        self.sumPoint()
        self.maxPoint()

    def showBurst(self, ax, beg=0, end=-1):
        beg = int(beg*1e6/self.dtBin)
        if end == -1:
            end = self.traceNum
        else:
            end = int(end*1e6/self.dtBin)

        chl_a_select = self.chl_a_bin.copy().astype('float')
        chl_b_select = self.chl_b_bin.copy().astype('float')
        chl_a_select[~self.selectIndBin] = np.nan
        chl_b_select[~self.selectIndBin] = np.nan

        ax.plot(self.time[beg:end],  self.chl_a_bin[beg:end], color='grey')
        ax.plot(self.time[beg:end], -self.chl_b_bin[beg:end], color='grey')
        ax.plot(self.time[beg:end],  chl_a_select[beg:end], alpha=.7, color='c')
        ax.plot(self.time[beg:end], -chl_b_select[beg:end], alpha=.7, color='m')
        ax.set_xlabel('time / s')
        ax.set_ylabel('corrected fluorescence counts')

class Fret(Burst):
    'FRET of Dr and Ac channels'
    def __init__(self, chl_dr_bin, chl_ac_bin, dtBin, eta, cross):
        Burst.__init__(self, chl_dr_bin, chl_ac_bin, dtBin, eta, cross)
        self.x_ini = []
        self.y_ini = []

    def getFret(self):
        Burst.getBurst(self)
        self.fret = 1.*self.b_peak_sum/(self.a_peak_sum+self.b_peak_sum)
        self.fret_point = 1.*self.b_point_sum/(self.a_point_sum+self.b_point_sum)
        self.fret_fake = 1.*self.chl_b_bin[self.peak_bool]/(self.chl_a_bin[self.peak_bool]+self.chl_b_bin[self.peak_bool])

    def showFret(self, ax, norm=0):
        self.yhist, self.xhist, patches = ax.hist(self.fret, bins=np.arange(-.2, 1.2, .02), normed=norm, ec='black')
        self.xhist = (self.xhist[:-1]+self.xhist[1:])*.5
        ax.text(.2, max(self.yhist*.8), '#molecule: '+str(len(self.fret)))
        ax.set_xlim((-.2, 1.2))
        if norm:
            ax.set_ylim((0, 10))
        ax.set_xlabel('FRET efficiency')
        ax.set_ylabel('# bursts')

    def showFret_fake(self, ax):
        self.yhist_fake, self.xhist_fake, patches = ax.hist(self.fret_fake, bins=np.arange(-.2, 1.2, .02))
        self.xhist_fake = (self.xhist_fake[:-1]+self.xhist_fake[1:])*.5
        ax.text(.3, max(self.yhist*1), '#molecule: '+str(len(self.fret_fake)))
        ax.set_xlim((-.2, 1.2))
        ax.set_xlabel('fake FRET efficiency')
        ax.set_ylabel('events')

    def fret_click(self, event):
        event.inaxes.plot(event.xdata, event.ydata, 'v', color='black')
        self.x_ini.append(event.xdata)
        self.y_ini.append(event.ydata)
        event.canvas.draw()

    def fret_fit(self, ax, vary):
        dis_peak_n = len(self.x_ini)
        params = Parameters()

        if dis_peak_n == 1:
            params.add('h', value = self.y_ini[0]) # height
            params.add('c', value = self.x_ini[0]) # center
            params.add('w', value = .05)           # width
            #params.add('o', value = 0)             # offset
            minner = Minimizer(one_gau, params, fcn_args=(self.xhist, self.yhist))

        if dis_peak_n == 2:
            params.add('h1', value = self.y_ini[0])
            params.add('h2', value = self.y_ini[1])
            params.add('c1', value = self.x_ini[0])
            params.add('c2', value = self.x_ini[1], vary=vary)
            params.add('w1', value = .05)
            params.add('w2', value = .05)
            params.add('o', value = 0)
            minner = Minimizer(two_gau, params, fcn_args=(self.xhist, self.yhist))

        if dis_peak_n == 3:
            params.add('h1', value = self.y_ini[0])
            params.add('h2', value = self.y_ini[1])
            params.add('h3', value = self.y_ini[2])
            params.add('c1', value = self.x_ini[0])
            params.add('c2', value = self.x_ini[1], vary=vary)
            params.add('c3', value = self.x_ini[2])
            params.add('w1', value = .05)
            params.add('w2', value = .05)
            params.add('w3', value = .05)
            params.add('o', value = 0)
            minner = Minimizer(three_gau, params, fcn_args=(self.xhist, self.yhist))

        self.fret_fit_result = minner.minimize()

        self.yhist_fit = self.yhist+self.fret_fit_result.residual
        ax.plot(self.xhist, self.yhist_fit, 'black')
        ax.set_ylabel('Frequency')

        if dis_peak_n == 1:
            params1 = Parameters()
            params1.add('h', value = self.fret_fit_result.params['h'].value)
            params1.add('c', value = self.fret_fit_result.params['c'].value)
            params1.add('w', value = self.fret_fit_result.params['w'].value)
            #params1.add('o', value = self.fret_fit_result.params['o'].value)

            params1.add('h_std', value = self.fret_fit_result.params['h'].stderr)
            params1.add('c_std', value = self.fret_fit_result.params['c'].stderr)
            params1.add('w_std', value = self.fret_fit_result.params['w'].stderr)
            #params1.add('o_std', value = self.fret_fit_result.params['o'].stderr)

            self.fit_notes = '\npeak 1:\n' \
                             +'center: '+str(round(params1['c'].value, 2))+'+/-'+str(round(params1['c_std'].value, 2))+' ('+str(round(params1['c_std'].value/params1['c'].value, 2)*100)+'%)\n' \
                             +'height: '+str(round(params1['h'].value, 2))+'+/-'+str(round(params1['h_std'].value, 2))+' ('+str(round(params1['h_std'].value/params1['h'].value, 2)*100)+'%)\n' \
                             +'sigma: '+str(round(params1['w'].value, 2))+'+/-'+str(round(params1['w_std'].value, 2))+' ('+str(round(params1['w_std'].value/params1['w'].value, 2)*100)+'%)\n' \
                             +'area: '+str(round(params1['h'].value*params1['w'].value/0.3989, 2))+'\n\n' \
                             +'reduced chisqr: '+str(self.fret_fit_result.redchi)
                             #+'offset: '+str(round(params1['o'].value, 2))+'+/-'+str(round(params1['o_std'].value, 2))+' ('+str(round(params1['o_std'].value/params1['o'].value, 2)*100)+'%)\n' \

        if dis_peak_n == 2:
            params1 = Parameters()
            params2 = Parameters()

            params1.add('h', value = self.fret_fit_result.params['h1'].value)
            params1.add('c', value = self.fret_fit_result.params['c1'].value)
            params1.add('w', value = self.fret_fit_result.params['w1'].value)

            params1.add('h_std', value = self.fret_fit_result.params['h1'].stderr)
            params1.add('c_std', value = self.fret_fit_result.params['c1'].stderr)
            params1.add('w_std', value = self.fret_fit_result.params['w1'].stderr)

            params2.add('h', value = self.fret_fit_result.params['h2'].value)
            params2.add('c', value = self.fret_fit_result.params['c2'].value)
            params2.add('w', value = self.fret_fit_result.params['w2'].value)

            params2.add('h_std', value = self.fret_fit_result.params['h2'].stderr)
            params2.add('c_std', value = self.fret_fit_result.params['c2'].stderr)
            params2.add('w_std', value = self.fret_fit_result.params['w2'].stderr)

            params1.add('o', value = self.fret_fit_result.params['o'].value)
            params2.add('o', value = self.fret_fit_result.params['o'].value)
            params1.add('o_std', value = self.fret_fit_result.params['o'].stderr)
            params2.add('o_std', value = self.fret_fit_result.params['o'].stderr)

            self.fit_notes = '\npeak 1:\n' \
                             +'center: '+str(round(params1['c'].value, 2))+'+/-'+str(round(params1['c_std'].value, 2))+' ('+str(round(params1['c_std'].value/params1['c'].value, 2)*100)+'%)\n' \
                             +'height: '+str(round(params1['h'].value, 2))+'+/-'+str(round(params1['h_std'].value, 2))+' ('+str(round(params1['h_std'].value/params1['h'].value, 2)*100)+'%)\n' \
                             +'sigma: '+str(round(params1['w'].value, 2))+'+/-'+str(round(params1['w_std'].value, 2))+' ('+str(round(params1['w_std'].value/params1['w'].value, 2)*100)+'%)\n' \
                             +'area: '+str(round(params1['h'].value*params1['w'].value/0.3989, 2)) \
                             +'\n\npeak 2:\n' \
                             +'center: '+str(round(params2['c'].value, 2))+'+/-'+str(round(params2['c_std'].value, 2))+' ('+str(round(params2['c_std'].value/params2['c'].value, 2)*100)+'%)\n' \
                             +'height: '+str(round(params2['h'].value, 2))+'+/-'+str(round(params2['h_std'].value, 2))+' ('+str(round(params2['h_std'].value/params2['h'].value, 2)*100)+'%)\n' \
                             +'sigma: '+str(round(params2['w'].value, 2))+'+/-'+str(round(params2['w_std'].value, 2))+' ('+str(round(params2['w_std'].value/params2['w'].value, 2)*100)+'%)\n' \
                             +'area: '+str(round(params2['h'].value*params2['w'].value/0.3989, 2))+'\n\n' \
                             +'reduced chisqr: '+str(self.fret_fit_result.redchi) +'\n\n' \
                             +'offset: '+str(round(params1['o'].value, 2))+'+/-'+str(round(params1['o_std'].value, 2))+' ('+str(round(params1['o_std'].value/params1['o'].value, 2)*100)+'%)\n' \
    


            mask1 = (self.xhist>params1['c']-params1['w']*3) & (self.xhist<params1['c']+params1['w']*3)
            mask2 = (self.xhist>params2['c']-params2['w']*3) & (self.xhist<params2['c']+params2['w']*3)

            ax.plot(self.xhist[mask1], one_gau(params1, self.xhist[mask1], self.yhist[mask1])+self.yhist[mask1], color='blue')
            ax.plot(self.xhist[mask2], one_gau(params2, self.xhist[mask2], self.yhist[mask2])+self.yhist[mask2], color='blue')

        if dis_peak_n == 3:
            params1 = Parameters()
            params2 = Parameters()
            params3 = Parameters()

            params1.add('h', value = self.fret_fit_result.params['h1'].value)
            params1.add('c', value = self.fret_fit_result.params['c1'].value)
            params1.add('w', value = self.fret_fit_result.params['w1'].value)

            params1.add('h_std', value = self.fret_fit_result.params['h1'].stderr)
            params1.add('c_std', value = self.fret_fit_result.params['c1'].stderr)
            params1.add('w_std', value = self.fret_fit_result.params['w1'].stderr)

            params2.add('h', value = self.fret_fit_result.params['h2'].value)
            params2.add('c', value = self.fret_fit_result.params['c2'].value)
            params2.add('w', value = self.fret_fit_result.params['w2'].value)

            params2.add('h_std', value = self.fret_fit_result.params['h2'].stderr)
            params2.add('c_std', value = self.fret_fit_result.params['c2'].stderr)
            params2.add('w_std', value = self.fret_fit_result.params['w2'].stderr)

            params3.add('h', value = self.fret_fit_result.params['h3'].value)
            params3.add('c', value = self.fret_fit_result.params['c3'].value)
            params3.add('w', value = self.fret_fit_result.params['w3'].value)

            params3.add('h_std', value = self.fret_fit_result.params['h3'].stderr)
            params3.add('c_std', value = self.fret_fit_result.params['c3'].stderr)
            params3.add('w_std', value = self.fret_fit_result.params['w3'].stderr)

            params1.add('o', value = self.fret_fit_result.params['o'].value)
            params2.add('o', value = self.fret_fit_result.params['o'].value)
            params3.add('o', value = self.fret_fit_result.params['o'].value)
            params1.add('o_std', value = self.fret_fit_result.params['o'].stderr)
            params2.add('o_std', value = self.fret_fit_result.params['o'].stderr)
            params3.add('o_std', value = self.fret_fit_result.params['o'].stderr)

            self.fit_notes = '\npeak 1:\n' \
                             +'center: '+str(round(params1['c'].value, 2))+'+/-'+str(round(params1['c_std'].value, 2))+' ('+str(round(params1['c_std'].value/params1['c'].value, 2)*100)+'%)\n' \
                             +'height: '+str(round(params1['h'].value, 2))+'+/-'+str(round(params1['h_std'].value, 2))+' ('+str(round(params1['h_std'].value/params1['h'].value, 2)*100)+'%)\n' \
                             +'sigma: '+str(round(params1['w'].value, 2))+'+/-'+str(round(params1['w_std'].value, 2))+' ('+str(round(params1['w_std'].value/params1['w'].value, 2)*100)+'%)\n' \
                             +'area: '+str(round(params1['h'].value*params1['w'].value/0.3989, 2)) \
                             +'\n\npeak 2:\n' \
                             +'center: '+str(round(params2['c'].value, 2))+'+/-'+str(round(params2['c_std'].value, 2))+' ('+str(round(params2['c_std'].value/params2['c'].value, 2)*100)+'%)\n' \
                             +'height: '+str(round(params2['h'].value, 2))+'+/-'+str(round(params2['h_std'].value, 2))+' ('+str(round(params2['h_std'].value/params2['h'].value, 2)*100)+'%)\n' \
                             +'sigma: '+str(round(params2['w'].value, 2))+'+/-'+str(round(params2['w_std'].value, 2))+' ('+str(round(params2['w_std'].value/params2['w'].value, 2)*100)+'%)\n' \
                             +'area: '+str(round(params2['h'].value*params2['w'].value/0.3989, 2)) \
                             +'\n\npeak 3:\n' \
                             +'center: '+str(round(params3['c'].value, 2))+'+/-'+str(round(params3['c_std'].value, 2))+' ('+str(round(params3['c_std'].value/params3['c'].value, 2)*100)+'%)\n' \
                             +'height: '+str(round(params3['h'].value, 2))+'+/-'+str(round(params3['h_std'].value, 2))+' ('+str(round(params3['h_std'].value/params3['h'].value, 2)*100)+'%)\n' \
                             +'sigma: '+str(round(params3['w'].value, 2))+'+/-'+str(round(params3['w_std'].value, 2))+' ('+str(round(params3['w_std'].value/params3['w'].value, 2)*100)+'%)\n' \
                             +'area: '+str(round(params3['h'].value*params3['w'].value/0.3989, 2))+'\n\n' \
                             +'reduced chisqr: '+str(self.fret_fit_result.redchi)+'\n\n' \
                             +'offset: '+str(round(params1['o'].value, 2))+'+/-'+str(round(params1['o_std'].value, 2))+' ('+str(round(params1['o_std'].value/params1['o'].value, 2)*100)+'%)\n' \
            
            mask1 = (self.xhist>params1['c']-params1['w']*3) & (self.xhist<params1['c']+params1['w']*3)
            mask2 = (self.xhist>params2['c']-params2['w']*3) & (self.xhist<params2['c']+params2['w']*3)
            mask3 = (self.xhist>params3['c']-params3['w']*3) & (self.xhist<params3['c']+params3['w']*3)

            ax.plot(self.xhist[mask1], one_gau(params1, self.xhist[mask1], self.yhist[mask1])+self.yhist[mask1], color='blue')
            ax.plot(self.xhist[mask2], one_gau(params2, self.xhist[mask2], self.yhist[mask2])+self.yhist[mask2], color='blue')
            ax.plot(self.xhist[mask3], one_gau(params3, self.xhist[mask3], self.yhist[mask3])+self.yhist[mask3], color='blue')

class FretDist(Fret):
    'FRET distribution from txt hitogram files'
    def __init__(self, xhist, yhist):
        self.xhist = xhist
        self.yhist = yhist
        self.x_ini = []
        self.y_ini = []

    def showFret(self, ax):
        ax.bar(self.xhist, self.yhist, align='center', width=.02, color='#E24A33', ec='black')
        ax.set_ylim((0, 10))
        ax.set_xlabel('FRET efficiency')
        ax.set_ylabel('frequency')

class Aniso(Burst):
    'Anisotropy of parallel and senkrecht channels'
    def __init__(self, chl_s_bin, chl_p_bin, dtBin, eta, cross):
        Burst.__init__(self, chl_s_bin, chl_p_bin, dtBin, eta, cross)

    def getAniso(self):
        Burst.getBurst(self)
        self.aniso = 1.*(self.a_peak_sum-self.b_peak_sum)/(self.a_peak_sum+2*self.b_peak_sum)

class Fcs():
    def __init__(self, trace_h5_path, dt, lag_end = 20):
        # lag_end 20: 2**20 ~= 1 s, after quasi-exponential ~0.6 s
        self.trace_h5_path = trace_h5_path
        self.dt = dt
        self.lag_end = lag_end
        self.n_power = np.int(lag_end-5-np.log2(dt))
        self.channel = 2**(self.n_power+6)

    def corrFFT(self, trace1, trace2):
        # trace1 cor trace2
        # trace1: correlated trace n(t)
        # trace2: lagged trace n(t+tau)
        # channel: number of points of correlated function
        # Max delay time before quasi-exponential bin: channel*dt/1e6 s

        FtraceCor = np.fft.rfft(trace1[:self.channel/2], n=self.channel) #n(t)
        FtraceCor = np.conj(FtraceCor)
        FtraceLag = np.fft.rfft(trace2, n=self.channel) #n(t+tau)
        Fcorr = FtraceCor*FtraceLag
        del FtraceCor, FtraceLag
        corr = np.fft.irfft(Fcorr)[1:self.channel/2]/(self.channel/2)
        del Fcorr
        return corr

    def chunkFcs(self, accc, single_resolve):
        # single_resolve: select specific peaks and fill noises elsewhere if single_resolve is True.
        if single_resolve:
            self.G = np.zeros(self.channel/2-1)
            read_i = 0
            fcs_i = 0
            trace_append_flag = 0
            mean_cor = 0
            mean_lag = 0

            while (read_i+1)*self.channel/2 <= len(self.selectInd):
                start = read_i*self.channel/2
                stop = (read_i+1)*self.channel/2
                temp_select = self.selectInd[start:stop]

                if np.count_nonzero(temp_select):
                    if trace_append_flag:
                        trace_temp = pd.read_hdf(self.trace_h5_path, key='trace', start=start, stop=stop)
                        noise_temp = pd.read_hdf(self.noise_h5_path, key='noise', start=start, stop=stop)

                        tempA_tail = trace_temp.A.values
                        tempB_tail = trace_temp.B.values
                        tempA_tail[~temp_select] = noise_temp.A.values[~temp_select]
                        tempB_tail[~temp_select] = noise_temp.B.values[~temp_select]
                        # tempA_tail[~temp_select] = np.random.poisson(self.a_noise_rate/1e3*self.dt, np.count_nonzero(~temp_select))
                        # tempB_tail[~temp_select] = np.random.poisson(self.b_noise_rate/1e3*self.dt, np.count_nonzero(~temp_select))

                        tempA = np.append(tempA_head, tempA_tail)
                        tempB = np.append(tempB_head, tempB_tail)

                        if accc == 'ab':
                            temp_cor = tempA
                            temp_lag = tempB
                        if accc == 'ba':
                            temp_cor = tempB
                            temp_lag = tempA
                        if accc == 'aa':
                            temp_cor = tempA
                            temp_lag = tempA
                        if accc == 'bb':
                            temp_cor = tempB
                            temp_lag = tempB
                        if accc == 'tt':
                            temp_cor = tempA+tempB
                            temp_lag = tempA+tempB

                        self.G = self.G+self.corrFFT(temp_cor, temp_lag)
                        mean_cor = mean_cor+np.mean(temp_cor)
                        mean_lag = mean_lag+np.mean(temp_lag)
                        fcs_i = fcs_i+1
                        print fcs_i

                        tempA_head = tempA_tail
                        tempB_head = tempB_tail
                        read_i = read_i+1

                    else:
                        trace_temp = pd.read_hdf(self.trace_h5_path, key='trace', start=start, stop=stop)
                        noise_temp = pd.read_hdf(self.noise_h5_path, key='noise', start=start, stop=stop)

                        tempA_head = trace_temp.A.values
                        tempB_head = trace_temp.B.values

                        tempA_head[~temp_select] = noise_temp.A.values[~temp_select]
                        tempB_head[~temp_select] = noise_temp.B.values[~temp_select]
                        # tempA_head[~temp_select] = np.random.poisson(self.a_noise_rate/1e3*self.dt, np.count_nonzero(~temp_select))
                        # tempB_head[~temp_select] = np.random.poisson(self.b_noise_rate/1e3*self.dt, np.count_nonzero(~temp_select))

                        trace_append_flag = 1
                        read_i = read_i+1
                else:
                    read_i = read_i+1

            self.G = self.G*fcs_i/mean_cor/mean_lag

        else:
            self.G = np.zeros(self.channel/2-1)
            i = 0
            mean_cor = 0
            mean_lag = 0

            temp_trace = pd.read_hdf(self.trace_h5_path, key='trace', start=i*self.channel/2, stop=(i+2)*self.channel/2)
            while len(temp_trace) == self.channel:
                print 'fcs interval: '+str(i)
                if accc=='ab':
                    temp_cor = temp_trace.A.values
                    temp_lag = temp_trace.B.values

                if accc=='ba':
                    temp_cor = temp_trace.B.values
                    temp_lag = temp_trace.A.values

                if accc=='aa':
                    temp_cor = temp_trace.A.values
                    temp_lag = temp_trace.A.values

                if accc=='bb':
                    temp_cor = temp_trace.B.values
                    temp_lag = temp_trace.B.values

                if accc=='tt':
                    temp_cor = temp_trace.A.values+temp_trace.B.values
                    temp_lag = temp_trace.A.values+temp_trace.B.values

                self.G = self.G+self.corrFFT(temp_cor, temp_lag)
                mean_cor = mean_cor+np.mean(temp_cor)
                mean_lag = mean_lag+np.mean(temp_lag)

                i = i+1
                temp_trace = pd.read_hdf(self.trace_h5_path, key='trace', start=i*self.channel/2, stop=(i+2)*self.channel/2)

            self.G = self.G*i/meanA/meanB
        self.T = (np.arange(self.G.size)+1)/1e6*self.dt

    def quasiExponent(self):
        #first 32 (2**5) points: dt
        #next 16 points: 2*dt
        #...
        #final 16 points: 2**n*dt
        #the end time = 2**(n+5)*dt

        self.BinG = np.zeros(self.n_power*16+32)
        self.BinT = np.zeros(self.n_power*16+32)
        self.BinG[:32] = self.G[:32]
        self.BinT[:32] = self.T[:32]
        for binFold in np.arange(self.n_power)+1:
            BegIndex = 2**(4+binFold)
            EndIndex = 2**(5+binFold)
            if EndIndex <= self.G.size:
                self.BinG[32+16*(binFold-1):32+16*binFold] = np.mean(self.G[BegIndex:EndIndex].reshape((-1, 16), order='F'), axis=0)
                self.BinT[32+16*(binFold-1):32+16*binFold] = np.mean(self.T[BegIndex:EndIndex].reshape((-1, 16), order='F'), axis=0)
            else:
                EndIndex = EndIndex-2**binFold
                self.BinG[32+16*(binFold-1):31+16*binFold] = np.mean(self.G[BegIndex:EndIndex].reshape((-1, 15), order='F'), axis=0)
                self.BinT[32+16*(binFold-1):31+16*binFold] = np.mean(self.T[BegIndex:EndIndex].reshape((-1, 15), order='F'), axis=0)
                self.BinG = np.delete(self.BinG, -1)
                self.BinT = np.delete(self.BinT, -1)

        print 'FCS calculated!\n'

    def getFcs(self, accc='tt', single_resolve=False):
        self.chunkFcs(accc, single_resolve)
        self.quasiExponent()

    def showFcs(self):
        plt.semilogx(self.BinT, self.BinG)
        plt.xlabel('lag time / s')
        plt.ylabel('correlation')
        plt.show()

# single-molecule resolve FCS
class srFcs(Fret, Fcs):
    def __init__(self, trace_h5_path, noise_h5_path, chl_a_bin, chl_b_bin, dt, dtBin, eta, cross, lag_end = 20):
        Fret.__init__(self, chl_a_bin, chl_b_bin, dtBin, eta, cross)
        Fcs.__init__(self, trace_h5_path, dt, lag_end)
        self.noise_h5_path = noise_h5_path
        # lag_end 20: 2**20 ~= 1 s, after quasi-exponential ~0.6 s

    def getFcs(self, accc='tt', single_resolve=True):
        self.chunkFcs(accc, single_resolve)
        self.quasiExponent()

    def showFcs(self):
        plt.semilogx(self.BinT, self.BinG)
        plt.xlabel('lag time / s')
        plt.ylabel('correlation')
        plt.show()
