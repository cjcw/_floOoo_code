import numpy as np
import matplotlib.pyplot as plt
import json
import os
import pandas as pd
from scipy import stats


def help_figplot():
    '''plt.figure(figsize=(wTot, hTot))
    grid = plt.GridSpec(hTot, wTot, hspace=hspace, wspace=wspace)
    plt.subplot(grid[:, currW:(currW+w)], facecolor=facecolor)
    '''
    print('see docstring')

def help_tsne():
    '''
    emb = TSNE(n_components=2, perplexity=30, early_exaggeration=5,
    learning_rate=10, n_iter=1000).fit_transform(obs)
    '''
    print('see docstring')

### --- DEFAULTS --- ###
print('WARNING: sample_rate to check with Ore')
sample_rate = 1000. # fix - check this
channels = ['CP3', 'C3', 'F5', 'PO3', 'PO4', 'F6', 'C4', 'CP4']
analysisDir_temp = '/Users/charl/OneDrive/Documents/__water/floOoo/code/analysis/neurofusion/_temp/'



### --- LOAD --- ###
def load_defaults():
    '''
    returns:
    sample_rate, channels, analysisDir_temp
    '''
    return sample_rate, channels, analysisDir_temp

def load_events(rootDir):
    events = pd.read_csv(rootDir+'events'+'.csv')
    return events




### --- CLASSES --- ###

class Filename_RawBrainwave:
    def __init__(self, folder, filename):
        self.folder = folder
        self.filename = filename
        self.fileID = int(filename.split('_')[1].split('.')[0])
        self.channels = channels
        self.sample_rate = sample_rate

    def load_chData(self):
        data = pd.read_csv(self.folder + self.filename)
        chData = np.row_stack([data[channel] for channel in self.channels])
        return chData




### --- GET --- ###
def get_filenames(rootDir, subDirStr='rawBrainwaves'):
    '''
    returns:
    folder, filenames
    '''
    folder = rootDir+subDirStr+'/'
    filenames = os.listdir(folder)
    if subDirStr == 'rawBrainwaves':
        def check_file(filename, subDirStr=subDirStr):
            '''
            checks for expected filetype
            fix - add filesize threshold?
            '''
            checks = []
            if subDirStr == 'rawBrainwaves':
                checks.append(filename[-4:] == '.csv')
                checks.append('rawBrainwaves' in filename)
                checks.append(len(filename.split('_')) == 2)
            else:
                raise ValueError('subDirStr not recognised')
            return all(checks)
    else:
        raise ValueError('subDirStr not recognised')
    filenames = [filename for filename in filenames if check_file(filename)]
    return folder, filenames








### --- UTILS --- ###

# --- get
def get_psd(X, sample_rate):
    """
    !! Warning: This function may need to be modified so that the returned variable freqAx_psd covers a suitable frequency range
    Computes and return the Power Spectral Density (PSD) estimate of a signal.

    Parameters
    ----------
    X : ndarray
        1D array signal
    sample_rate : float
        The sampling rate of X

    Returns
    ----------
    freqAx_psd, psd

    freqAx_psd : ndarray
        1D array frequency axis
    psd : ndarray
        1D array PSD estimate
    """

    window = 'hann' # use a Hanning window for Fourier transform

    def get_psd_(X, sample_rate, maxFreq, pointsPerHz):
        from scipy.signal import welch
        psdIndMax = int(maxFreq*pointsPerHz)
        psd = welch(X, fs=sample_rate, window=window, nperseg=int(sample_rate)*pointsPerHz)
        freqAx_psd, psd = psd[0][0:psdIndMax], psd[1][0:psdIndMax]
        return freqAx_psd, psd

    freqAx_psd, psd = [], []
    maxFreqs = [1, 10, 100, 200, 500]
    pointsPerHzs = [20, 4, 1, 0.04, 0.02]
    for maxFreq, pointsPerHz in zip(maxFreqs, pointsPerHzs):
        freqAx_psd_, psd_ = get_psd_(X, sample_rate, maxFreq, pointsPerHz)
        if not len(freqAx_psd):
            i = 0
        else:
            i = np.flatnonzero(freqAx_psd_ > freqAx_psd[-1][-1])[0]
        freqAx_psd.append(freqAx_psd_[i:])
        psd.append(psd_[i:])
    freqAx_psd = np.concatenate(freqAx_psd)
    psd = np.concatenate(psd)

    return freqAx_psd, psd



def edges2cens(edges):
    cens = np.convolve(edges,[.5,.5],'same')
    cens = cens[1::]
    return cens

def cens2edges(cens):
    d=np.diff(cens)/2.
    edges = np.add(cens, d[0])
    edges = np.concatenate([[cens[0]-d[0]], edges])
    return edges



def get_hist(data, bins, norm=False):

    counts, edges = np.histogram(data, bins)
    cens = edges2cens(edges)

    if norm:
        counts = counts / counts.sum(dtype=float)

    return cens, counts, edges

def get_xy_4_histFill(cens, Ys):
    edges = cens2edges(cens)

    x=[]
    y=[]
    for i in range(cens.shape[0]):
        if i:
            x.append([edges[i], edges[i]])
            y.append([Ys[i-1], Ys[i]])

    x = np.concatenate([[edges[0], edges[0], edges[1]], np.concatenate(x), [edges[-2], edges[-1], edges[-1]]])
    y = np.concatenate([[0, Ys[0], Ys[0]], np.concatenate(y), [Ys[-1], Ys[-1], 0]])

    return x, y


def get_motifStartInds(motif, arr):
    '''returns the start inds of a motif found in an array/list'''
    #
    #startInds = [ind for ind, val in enumerate(Arr) if val == motif[0]]
    #startInds = [ind for ind in startInds if ind <= len(Arr)-len(motif)]
    startInds = np.where(arr == motif[0])[0]
    startInds = startInds[np.where(startInds<=len(arr)-len(motif))]
    #
    motifStartInds = []
    #
    for candidateInd in startInds:
        #
        if arr[(candidateInd+(len(motif)-1))] == motif[-1]:
            candidateMotif = arr[candidateInd:(candidateInd+(len(motif)))]
            if np.array_equal(candidateMotif, motif):
                motifStartInds.append(candidateInd)
    #
    motifStartInds = np.array(motifStartInds)
    return motifStartInds


def get_boutWindows(boutStatus):
    boutStartInds = list(np.array(get_motifStartInds([False, True], boutStatus))+1)
    boutEndInds = list(np.array(get_motifStartInds([True, False], boutStatus)))
    if all([len(boutStartInds), len(boutEndInds)]): #techincall might want option to incliude start/end
        if boutEndInds[0] < boutStartInds[0]:
            boutEndInds = boutEndInds[1:]
        if boutStartInds[-1] > boutEndInds[-1]:
            boutEndInds.append(len(boutStatus)-1)
        bouts = np.column_stack([boutStartInds, boutEndInds])
    else:
        bouts = []

    return bouts

def get_nanMSE4obs(obs, med=False, std=False, axis=0):
    '''
    default will return:
    m, se
    '''
    def f_m(obs):
        f = [np.nanmean, np.nanmedian][med]
        return f(obs, axis=axis)

    if callable(std):
        f_err = std
    else:
        def f_se(obs):
            return stats.sem(obs, axis=axis, nan_policy='omit')

        def f_std(obs):
            return np.nanstd(obs, axis=axis)

        f_err = [f_se, f_std][std]

    m, err = [f(obs) for f in [f_m, f_err]]
    return m, err



# --- plot
def plot_title(s, fontsize=None, color=None, loc='left', fontweight='bold'):
    plt.title(s, fontsize=fontsize, color=color, loc=loc, fontweight=fontweight)

# --- save
def save_image(filename, saveDir=analysisDir_temp, ext='.png'):
    os.chdir(saveDir)
    plt.savefig(filename+ext)
    plt.close()
