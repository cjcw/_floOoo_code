import ccwvBaseFunctions as vbf
import os as os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random as rdm
import ccwspectra as spec
import pickle as pk
import json
import itertools
from math import*
import scipy
from scipy import stats
from scipy import signal as sig
import myemd
import seaborn as sb
import networkx as nx
#import emd as emd
# extra modules

import ccwdifferenceEstimation as de

from sklearn.decomposition import FastICA


emailAdress = 'charlieclarkewilliams@gmail.com'

emdStr='.opt_mEMD'
eegHz, sr = 1250., 1250.

cellTypeCols = {'glu' : '#2ebc4dff', 'da' : '#b84e5fff', 'gaba' : '#5971ecff'}

salCol = sb.color_palette('Blues', 100)[-5]
cocCol = sb.color_palette('Reds', 100)[-5]
highCol = '#fe9121ff'
lowCol = '#3fa8dbff'

gray1, gray2, gray3, gray4 = '#dadadaff', '#bbbbbbff', '#8d8d8dff', '#585858ff'
black = '#000000ff'
pastBlue = '#289ff7ff'
pastBlue2 = '#286df7f0'
pastRed = '#eb4346ff'
pastPurp = '#a466c9ff'
turq = '#66c9c8ff'

pulseCol_chR = '#33c5e9ff'
pulseCol_archT, pulseCol_archT2 = '#eeee85ff', '#eed285ff'


cocCol2='#b400a6ff'
salCol2='#d3d3d3ff'
regions = np.array(['pfc', 'nac', 'bla', '1', 'vta'])
regions_old = np.array(['1', 'bla', 'nac', 'pfc', 'vta'])
zones = np.array(['sal', 'coc'])
zones3 = np.array(['sal', 'coc', 'out'])
nRegs = len(regions)
led1Suffix = '.led1_pulse'
led2Suffix = '.led2_pulse'

salStateCol = '#b9b9b9ff'

salCol_LED = '#5fd35fff'
cocCol_LED = '#ffdf28ff'
 
salCol2_LED = '#008000ff' # darker
cocCol2_LED = '#ff9800ff' #

condCols_coc = sb.color_palette('Reds', 4)[1:]
condCols_sal = sb.color_palette('Blues', 4)[1:]


imfTits_suraTheta = ['beta', 'slow-gamma', 'mid-gamma', 'fast-gamma']


aniIDs_orig = np.array(['mccw03', 'mccw05', 'mhb17', 'mme18', 'mrr02', 'mrr03', 'mrr04', 'mrr06', 'mrr05'])
aniIDs_orig5site = np.array(['mhb17', 'mrr02', 'mrr04', 'mccw03', 'mccw05', 'mrr03', 'mrr05', 'mrr06'])
aniIDs_cl = np.array(['mccw06', 'mccw06_2', 'mccw07', 'mccw07_2', 'mccw09', 
                      'mccw10', 'mccw11', 'mccw11_2', 'mccw12', 'mccw13', 'mccw14', 
                      'mccw16s', 'mccw16', 'mccw17s', 'mccw17', 'mccw16_2', 'mccw18s', 'mccw18', 'mccw17_2', 'mccw20', 
                      # daGlu
                      'mccw21', 'mccw21_2', 'mccw21_3', 'mccw22', 'mccw22_2', 'mccw22s'])
aniIDs_cpp = np.array(['mccw18', 'mccw20', 'mccw23'])
aniIDs_mor = np.array(['mccw16', 'mccw17', 'mccw18', 'mccw23'])
aniIDs_phase = ['mccw11']
aniIDs_chR_zone = np.array(['mrr11', 'mrr12', 'mccw24', 'mccw25'])




bsnms = ['mhb17-171222', 'mrr02-180503', 'mrr04-180520', 'mccw03-180219',          
         'mccw05-180525', 'mrr03-180602', 'mrr05-181125', 'mrr06-190904']

bsnms_pre = ['mhb17-171220', 'mrr02-180501', 'mrr04-180517', 'mccw03-180217', 'mccw05-180522',
             'mrr03-180530', 'mrr06-190901']
stageBsnms = {}
stageBsnms['pre'] = bsnms_pre
stageBsnms['rec'] = bsnms
stageBsnms['ext'] = bsnms
stageBsnms['ren'] = bsnms

bsnms_chR = ['mrr11a-191106', 'mrr11b-191106', 'mrr11a-191107', 'mrr11b-191107', 'mrr11a-191108', 'mrr11b-191108', 'mrr11a-191109', 'mrr11b-191109', 
             'mrr12a-191106', 'mrr12a-191108', 'mrr12b-191108','mrr12c-191108', 'mrr12a-191109', 'mrr12b-191109']

bsnms_da = ['mdm42-0112-0127', 'mdm42-0312-0121', 'mdm42-0412-0120', 'mdm42-0512-0120', 
            'mdm42-0612-0120', 'mdm42-0712-0117', 'mdm42-2611-0115', 'mdm42-2711-0115', 
            'mdm42-2811-0117', 'mdm42-3011-0116']


bsnms_gaba = np.concatenate([['mccw15-220410'], 
                             ['mccw15'+s+'-220411' for s in ['a', 'b', 'c', 'd', 'e']], 
                             ['mccw15'+s+'-220414' for s in ['a', 'b', 'c', 'd', 'e']], 
                             ['mccw15a-220419'], 
                             ['mccw15'+s+'-220420' for s in ['a', 'b']], 
                             ['mccw15'+s+'-220421' for s in ['a', 'b']]])
                             

bsnms2 = ['mam5-160307', 'mam9-160524', 'mam10-170213']

bsnms_il = ['mccw05-180521']

bsnmsall = np.concatenate([bsnms, bsnms2, bsnms_chR, bsnms_il])
bsnmsallvta = np.concatenate([bsnms, bsnms2, bsnms_chR, bsnms_da])
bsnms_db = ['mvl15-210812', 'mvl15-210813']
bsnms_leila = ['sub056-220224']

rootFolder_rr = '/Dupret_Lab/merged/rrothaermel_merged/'
rootFolder_ccw = '/Dupret_Lab/merged/ccwilliams_merged/'
rootFolder_db = '/Dupret_Lab/merged/cNOR_merged/'
rootFolder_opto = '/Dupret_Lab/merged/optotagging_merged/'
rootFolder_steph = '/Dupret_Lab/merged/strouche_merged/'
rootFolder_mc = '/Dupret_Lab/merged/mcastelli_merged/'
rootFolder_vl = '/Dupret_Lab/merged/vlopes_merged/'


bsnms_rr = np.array(['mam9-160527', 'mccw05-180509', 'mccw05-180510', 'mccw05-180513', 'mccw05-180515',
                     'mhb17-171212', 'mhb17-171213', 'mhb17-171214', 
                     'mrr02-180417', 'mrr02-180418', 'mrr02-180419',
                     'mrr05-181024',
                     'mrr06-190712', 'mrr06-190720', 'mrr06-190721', 'mrr06-190722', 'mrr06-190723', 'mrr06-190724',
                     'mrr07-190920', 'mrr07-190921', 'mrr07-190922', 'mrr07-190924', 'mrr07-190925',
                     'mrr08-190927', 'mrr08-190928', 'mrr08-190929', 'mrr08-190930', 'mrr08-191001', 'mrr08-191002',
                     'mrr09-191029', 'mrr09-191030', 'mrr09-191031', 'mrr09-191101', 'mrr09-191103', 'mrr09-191104',
                     'mrr10-191029', 'mrr10-191030', 'mrr10-191031', 'mrr10-191101', 'mrr10-191102', 'mrr10-191103', 'mrr10-191104'])


bsnms_rr_nopulse = np.array(['mam9-160527', 'mccw05-180509', 'mccw05-180510', 'mccw05-180513', 'mccw05-180515', 'mhb17-171212', 'mhb17-171213', 'mhb17-171214', 'mrr02-180417', 'mrr02-180418', 'mrr02-180419', 'mrr05-181024'])

#np.array(['mdm195-170328', 'mdm195-170321', 'mdm195-170319', 'mdm195-170310', 'dm195-170308', 'mst973-170503', 'mdm189-170510'])
bsnms_steph_sucCPP = np.array(['mdm220-170731', 'mdm227-171031', 'mdm168-161102', 'mdm227-171110', 'mdm169-161104', 
                               'mst1095-170927', 'mdm195-170308', 'mdm227-171102', 'mst1095-170928', 'mst1095-170926',
                               'mdm169-161108', 'mdm227-171101'])
bsnms_steph_sucCPP_good = np.array(['mdm227-171031', 'mdm168-161102', 'mdm195-170308', 'mdm227-171102',
                                    'mst1095-170926', 'mdm169-161108'])

testStages_sucCPP = np.array(['pre', 'suc', 'wat', 'test'])


temdFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/tmEMD/'
cocrcFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/'
spks2imfsFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/basics/spks2imfs/' #'/Dupret_Lab/analysis/ccwilliams_analysis/figure_3/spks2imfs/'
cellTypesFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/figure_5/vta_cellTypes/'
fCorr2behFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/f_corr_2_beh/'
burstFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/figure_2/burst_analysis/'
sucroseFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/sucrose/'
conditioningFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/conditioning/'
extinctionFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/extinction/'
ledTrigSpecsFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/ledTrigSpecs/'
imfZonePercTroughTrigFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/imfZonePercTroughTrigX/'
spkZonePredFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/figure_3/spkZonePred/'
cpp_zoneTroughTrigFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/cpp/zoneTroughTrig/'
morphineFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/morphine/'
daGluFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/daGlu/'
cppFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/cpp/'
landOn4HzFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/landOn4Hz/'
raw_egsFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/raw_egs/'
harmonicsFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/coc_rc/harmonics/'


def set_nProcesses(nProcesses):
    '''
    cap the max nProcesses used
    '''
    os.environ['OMP_NUM_THREADS'] = str(nProcesses)


def get_nonDiagInds(n, mode='all'):
    '''
    returns the 2D indices of a square array not for the diagonal 
    mode : str
        'all' | 'upper' | 'lower'
        - to return upper+lower, upper, or lower triangles, respectively
    '''
    a, b = np.triu_indices(n, k=1)
    upper = np.row_stack([[a_, b_] for a_, b_ in zip(a, b)])

    a, b = np.tril_indices(n, k=-1)
    lower = np.row_stack([[a_, b_] for a_, b_ in zip(a, b)])

    if mode == 'all':
        nonDigInds = np.row_stack([upper, lower])
    elif mode == 'upper':
        nonDigInds = upper
    elif mode == 'lower':
        nonDigInds = lower
    else:
        raise ValueError("mode should be one of: 'all' | 'upper' | 'lower' ")
    return nonDigInds

def get_ellipse(x,y):
    '''
    Author: Guisseppe Gava
    '''
    xm = x.mean(); ym = y.mean()
    x -= xm; y -= ym
    U, S, V = np.linalg.svd(np.stack((x, y)))
    tt = np.linspace(0, 2*np.pi, 1000)
    circle = np.stack((np.cos(tt), np.sin(tt)))    # unit circle
    transform = np.sqrt(2/len(x)) * U.dot(np.diag(S))   # transformation matrix
    fit = transform.dot(circle) + np.array([[xm], [ym]])
    return fit

def getTitColForReg(reg, new=True):
    if reg == '1' or reg == 'ca1':
        titleOut = ['CA1', 'Hpc'][new] #'CA1'
        colOut = '#00cdffff'
    if reg == 'nac':
        titleOut = 'NAc'
        colOut = '#ff6900ff'
    if reg == 'bla':
        titleOut = 'Amy'
        colOut = '#ff2a7fff'
    if reg == 'vta':
        titleOut = 'VTA'
        colOut = '#55d700ff'
    if reg == 'pfc':
        titleOut = 'PFC'
        colOut = '#9955ffff'
    #
    #
    return titleOut, colOut

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False

def getAniID(bsnm):
    end_ind = bsnm.find('-')
    aniID = bsnm[0:end_ind]

    if is_number(aniID[-1]):
        return aniID
    else:
        for i in range(end_ind):
            aniID = bsnm[0:(end_ind-i)]
            if is_number(aniID[-1]):
                break 

    return aniID


bsnms_archT = np.array([b for b in bsnms_rr if getAniID(b) in ['mrr06', 'mrr07', 'mrr08']])




### CCW COC PAPER ANIS ##
bsnms_pre = ['mhb17-171220', 'mrr02-180501', 'mrr04-180517', 'mccw03-180217', 'mccw05-180522', 
             'mrr03-180530', 'mrr06-190901', 'mme18-180423']

bsnms_post = ['mhb17-171222', 'mrr02-180503', 'mrr04-180520', 'mccw03-180219', \
              'mccw05-180525', 'mrr03-180602', 'mrr06-190904', 'mme18-180425']
aniBsnms = {}
aniIDs_pre = [getAniID(bsnm) for bsnm in bsnms_pre]
aniIDs_post = [getAniID(bsnm) for bsnm in bsnms_post]
for ai, aniID in enumerate(aniIDs_pre):
    aj = aniIDs_post.index(aniID)
    #
    aniBsnms[aniID] = {'pre' : bsnms_pre[ai], 'post' : bsnms_post[aj]}
#
#aniIDs = np.intersect1d(aniIDs_pre, aniIDs_post)
bsnms_beh = np.array([aniBsnms['mme18'][testStage] for testStage in aniBsnms['mme18']])

mse = [np.mean, stats.sem]
xyt = [plt.xticks, plt.yticks]
xyl = [plt.xlabel, plt.ylabel]
xyli = [plt.xlim, plt.ylim]


'''
---------------------- NOTES ------------------------------------------------------
# LOG BINS 4 pdist
a=5
b=10**a
nHistBins=100
histBins_ = np.logspace(0.000001, a, nHistBins)
histBins_ /= b
'''



def help_legend():
    for s in ["l = plt.legend()", "for text in l.get_texts():", "text.set_color(color)"]:
        print(s)


    

def help_load():
    for s in ["rootFolder = ccw.get_rootFolder4bsnm(bsnm)", 
              "tetRegs, tetLabs, seshLabs, tetUseInds = ccw.get_bsnmInfo(bsnm, rootFolder)"]:
        print(s)

def help_tsne():
    print("emb = TSNE(n_components=2, perplexity=30, early_exaggeration=5, learning_rate=10, n_iter=1000).fit_transform(classDat)")

def help_kmeans():
    print("kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(pca.components_.T)")
    
def help_readBinary():
    for s in ["size = os.path.getsize(datfname)", "size = int(size/np.dtype(dtype).itemsize)",
              "chTraces = np.memmap(datfname, mode='r', dtype=dtype, order='F', shape=(nCh, int(size/nCh)))"]:
        print(s)

def set_plotStyle(i=0):
    s = ['Solarize_Light2', 'dark_background', 'fivethirtyeight'][i]
    print(s)
    plt.style.use(s)


def help_resetTheme():
    print("!jt -t chesterish -T -cellw 80%")

def help_figplot():
    for s in ['plt.figure(figsize=(wTot, hTot))', 
              'grid = plt.GridSpec(hTot, wTot, hspace=3, wspace=3)', 
              'plt.subplot(grid[currH:(currH+h_1r), currW:(currW+w_1)], facecolor=facecolor)']:
        print(s)

def help_instX():
    print("IP,IF,IA = emd.spectra.frequency_transform( imfs, sr, 'hilbert')")

def help_tickCol():
    print("plt.gca().get_xticklabels()[i].set_color(col)")

def help_optMaskFreqs():
    for s in ["from ccw_tailorEMD_functions import get_optMaskFreqs", "get_optMaskFreqs('nac', [300, 150, 0, 0], 5)"]:
        print(s)

def help_emd_spectra():
    print("IP,IF,IA = emd.spectra.frequency_transform( imfs, sr, 'hilbert' )")
        
def help_anova():
    '''
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    
    if ANCOVA:
        formula = 'amp ~ zone + C(bsnm) + speed'
        formula = 'speed ~ zone + C(bsnm)' 
        formula = 'speed ~ zone + speed'
        lm = ols(formula, data).fit()
        results = sm.stats.anova_lm(lm, typ=2)
        p = results['PR(>F)']['zone']
    else:
        formula = 'amp ~ zone + C(bsnm)
        lm = ols(formula, data).fit()
        results = sm.stats.anova_lm(lm, typ=2)
        p = results['PR(>F)']['zone']
    '''
    print('see docstring')

def help_multiprocessing():
    '''
    import multiprocessing
    import multiprocessing.pool
    class NoDaemonProcess(multiprocessing.Process):
        # make 'daemon' attribute always return False
        def _get_daemon(self):
            return False
        def _set_daemon(self, value):
            pass
        daemon = property(_get_daemon, _set_daemon)

    class _pool(multiprocessing.pool.Pool):
        Process = NoDaemonProcess

    def run_subIteration(args):
        a, b, c = args
        np.random.seed() # !! important for random processes
        work
        return __
    
    pool = _pool(nprocesses)
    args = (Xs, f_freqs, imfis_4_scoring, sample_rate, mask_args, mixScore_func, consistency_func, compute_consistency, 
            f_ranges0, nprocesses)
    it_outputs = pool.map(run_subIteration, [args for i in range(n_per_it)])
    pool.close()
    pool.join()

    it_mask_freqs = np.row_stack([it_outputsi[0] for it_outputsi in it_outputs])
    it_mix_scores = np.row_stack([it_outputsi[1] for it_outputsi in it_outputs])
    it_adj_mix_scores = np.array([it_outputsi[2] for it_outputsi in it_outputs])
    '''
    print('see docstring')
        
        
        
        
        
### FUNCTIONS ##############

opto_markers = {'da' : '^', 'glu' : 's', 'gaba' : 'o'}

def load_tmEMD_outputs(region, label='new', mixScoreStr='_corr'):
    
    if label == 'new':
        saveStr = '_pre-eEMD'
        optimise_method = 'm'
        if mixScoreStr in ['_corr', '_psds', '_psd']:
            print('ccw to check:     ccw.load_tmEMD_outputs')
            outputFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/tmEMD/tmEMD-JOSS_paper/'
        else:
            outputFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/tmEMD/tmEMD-JOSS_paper/example/'
        regTit, _ = getTitColForReg(region)
        #
        with open(outputFolder+'eg_tmEMD_outputs.'+regTit+saveStr+'_'+optimise_method+mixScoreStr+'.pkl', 'rb') as h:
            outputs = pk.load(h)

        it_mask_freqs, it_mix_scores, it_adj_mix_scores, it_consistency_scores, it_is, optimised_mask_freqs, converged, time_taken_secs = \
        [outputs[k] for k in ['it_mask_freqs', 'it_mix_scores', 'it_adj_mix_scores', 'it_consistency_scores', 'it_is', 'optimised_mask_freqs', 'converged', 'time_taken_secs']]
        
        extra_outputs = {}
        for x, k in zip([it_adj_mix_scores, it_is, optimised_mask_freqs, converged, time_taken_secs], 
                        ['it_adj_mix_scores', 'it_is', 'optimised_mask_freqs', 'converged', 'time_taken_secs']):
            extra_outputs[k] = x
    
    elif label == 'old':
        from ccw_tailorEMD_functions import get_optStr
        maskFreqs2add=np.array([300, 150, 0, 0])
        nMainFreqs = 5
        optStr = te.get_optStr(maskFreqs2add, nMainFreqs)

        outputFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/figure_2/tailorEMD/'
        it_mask_freqs = np.load(outputFolder+'anis_'+optStr+region+'.it_maskFreqs'+'.npy')
        it_mix_scores = np.load(outputFolder+'anis_'+optStr+region+'.it_mixScores'+'.npy')
        it_consistency_scores = np.load(outputFolder+'anis_'+optStr+region+'.it_consistencyScores'+'.npy')
        
        extra_outputs = None
        
    return it_mask_freqs, it_mix_scores, it_consistency_scores, extra_outputs
        
    
    
def plot_emd(lfp, imfs, amps=None, ses=None, st=12500, en=15500, timeAx=None, lfpCol='k', imfCols=None, ampCol='k', cmap='gray', 
             alpha=1, ls='-', lw_lfp=1, lw_imfs=1, lw_amps=2, spaceFactor=0.2, lfp_shift=0., imfs_shift=0., flipCols=False, 
             focusImfis=None, unFocusCol='gray', unFocusAlpha=0.3, unFocusLw=1, zorder=2, sr=eegHz, alpha_se=0.5, return_imfYs=False):
    
    """imfs is  [nImfs x t] array """
    
    
    if imfCols is None:
        try:
            imfCols = sb.color_palette(cmap, imfs.shape[1])
        except:
            imfCols = [cmap]*imfs.shape[1]
    #
    if flipCols:
        imfCols = imfCols[::-1]
    
    if lfp is None:
        lfp4plot = imfs[st:en, :].sum(axis=1)
        
    else:
        lfp4plot = lfp[st:en]
    
    if amps is not None:
        amps4plot = amps[st:en, :].T
    if timeAx is None:
        timeAx = np.linspace(0, len(lfp4plot)/sr, len(lfp4plot))
    
    plt.plot(timeAx, lfp4plot+lfp_shift, color=lfpCol, lw=lw_lfp, zorder=zorder)

    lfpMin, lfpMax = [f(lfp4plot) for f in [np.min, np.max]]
    emdYSt = lfpMin - (lfpMax-lfpMin)*spaceFactor + imfs_shift
    imfSpace = (lfpMax-lfpMin)*spaceFactor
    
    imfs4plot = imfs[st:en, :].T
    
    lfpMin, lfpMax = [f(lfp4plot) for f in [np.min, np.max]]
    emdYSt = lfpMin - (lfpMax-lfpMin)*spaceFactor + imfs_shift
    imfSpace = (lfpMax-lfpMin)*spaceFactor
    imfYs = []
    for imfi, imfTrace in enumerate(imfs4plot):
        if focusImfis is None:
            col = imfCols[imfi]
            alpha=1
            lw = lw_imfs
        else:
            if imfi in focusImfis:
                col = imfCols[imfi]
                alpha=alpha
                lw = lw_imfs
                zorder=3
            else:
                col = unFocusCol
                alpha = unFocusAlpha
                lw = unFocusLw
                zorder=2
        plt.plot(timeAx, imfTrace+emdYSt-(imfSpace*imfi), color=col, alpha=alpha, ls=ls, lw=lw, zorder=zorder)
        
        imfYs.append(emdYSt-(imfSpace*imfi))
        
        if ses is not None:
            #
            plt.fill_between(timeAx, imfTrace+emdYSt-(imfSpace*imfi)+ses[imfi, :], imfTrace+emdYSt-(imfSpace*imfi)-ses[imfi, :], color=col, alpha=alpha_se, zorder=zorder)
        
        if amps is not None:
            plt.plot(timeAx, amps4plot[imfi]+emdYSt-(imfSpace*imfi), color=ampCol, alpha=alpha, ls=ls, lw=lw_amps, zorder=zorder)
    plt.xlim(timeAx[0], timeAx[-1])
    if return_imfYs:
        return imfYs

    
    
    
def figplot_tmEMD(Xs, xi, it_mask_freqs, it_X_scores, sample_rate, show_variants=True, variants=['EMD', 'eEMD', 'itEMD'], k4score='m_corr', 
                  show_egs=True, window=None, eg_percs=[80, 30, 0], fontsize=16, ms=4, ms_=6, nSecs=30, title=None, color_title='w', opt2xi=False):
    
    import ccw_tmEMD as temd
    import emd
    
    facecolor=None
    #imfCols = sb.color_palette('Spectral', it_mask_freqs.shape[1])
    
    if opt2xi:
        it_X_scores_M = it_X_scores[:, xi] 
    else:
        it_X_scores_M = it_X_scores.mean(axis=1) 
    X = Xs[xi]
    
    if window is None:
        ll = int(sample_rate*nSecs)
        st = np.random.choice(np.arange(len(X)-ll))
        en = st+ll
    else:
        st, en = window

    w_freqs = [8, 6][show_egs]
    w_imfs = 20
    w_psd = 6

    if show_egs:
        eg_inds = np.array([findNearestInd(np.nanpercentile(np.unique(it_X_scores_M), p), it_X_scores_M) for p in eg_percs])
        wTot = w_freqs + w_imfs + w_psd
        h_eg = 7
        hTot = h_eg*len(eg_percs)
    else:
        eg_inds = np.array([])
        wTot = w_freqs
        hTot = 6

    
    
    plt.figure(figsize=(wTot, hTot))

    grid = plt.GridSpec(hTot, wTot, hspace=3, wspace=3)
    
    currW = 0
    ### ---------  PLOT mask freq space --------- ###
    plt.subplot(grid[:, currW:(currW+w_freqs)], facecolor=facecolor)
    plot_title(title, color=color_title)
    plt.xticks(fontsize=fontsize-2)
    plt.yticks(fontsize=fontsize-2)
    imfCols = temd.plot_mask_freq_scores(it_mask_freqs, it_X_scores, inds=eg_inds, ms=ms, color_='w', ms_=ms_)
    
    
    if show_variants:
        
        max_imfs = it_mask_freqs.shape[-1]
        fmin, fmax = 0, np.round(np.nanmax(it_mask_freqs), -1)
        
        for variant in variants:
            kScores = get_scores4emd(Xs, sample_rate, variant, temd.get_psd2)
            
            score = kScores[k4score].mean()
            
            plt.hlines(score, fmin, fmax, color='gray', zorder=1)
            #plt.text(0, mixScore+0.02, variant, color='gray')
            
                
    
    plt.xlabel('Mask freq. (Hz)', fontsize=fontsize)
    plt.ylabel('Mode mixing score (r)', fontsize=fontsize)
    plt.xscale('log')
    currW += w_freqs

    ### --------- PLOT eg mEMDs --------- ###
    if show_egs:
        currW0 = np.copy(currW)
        currH = 0
        for egi, ind in enumerate(eg_inds):

            mask_freqs = it_mask_freqs[ind]
            sift_config = emd.sift.get_config('mask_sift')
            sift_config['mask_freqs'] = mask_freqs/sample_rate
            sift_config['max_imfs'] = len(mask_freqs)

            currW = np.copy(currW0)
            imfs = emd.sift.mask_sift(X, **sift_config)
            freqAx_psd, imfPSDs = temd.get_imfPSDs(imfs, sample_rate)

            plt.subplot(grid[currH:(currH+h_eg), currW:(currW+w_imfs)], facecolor=facecolor)
            plt.xticks(fontsize=fontsize-2)
            plt.yticks([])
            plot_emd(None, imfs, st=st, en=en, sr=sample_rate, imfCols=imfCols, lfpCol='w', lw_imfs=2)
            if egi == len(eg_inds)-1:
                plt.xlabel('Time (s)', fontsize=fontsize)
            currW += w_imfs
            plt.subplot(grid[currH:(currH+h_eg), currW:(currW+w_psd)], facecolor=facecolor)
            plt.xticks(fontsize=fontsize-2)
            plot_imfPSDs(freqAx_psd, imfPSDs, imfCols=imfCols)
            if egi == len(eg_inds)-1:
                plt.xlabel('Freq. (Hz)', fontsize=fontsize)
            currW += w_psd

            currH += h_eg

        
'''
def figplot_tmEMD(X, it_mask_freqs, it_mix_scores, sample_rate, window=None, eg_percs=[80, 30, 0], fontsize=16, ms=4, ms_=6, nSecs=30, title=None, color_title='w', 
                  show_variants=True, variants=['EMD', 'eEMD', 'itEMD']):
    
    import ccw_tmEMD as temd
    import emd
    
    facecolor=None
    #imfCols = sb.color_palette('Spectral', it_mask_freqs.shape[1])
    it_mix_scores_M = it_mix_scores.mean(axis=1) #it_mix_scores[:, 0] 
    
    if window is None:
        ll=int(sample_rate*nSecs)
        st = np.random.choice(np.arange(len(X)-ll))
        en = st+ll
    else:
        st, en = window

    w_freqs = 6
    w_imfs = 20
    w_psd = 6

    wTot = w_freqs + w_imfs + w_psd
    h_eg = 7
    hTot = h_eg*len(eg_percs)
    eg_inds = np.array([findNearestInd(np.nanpercentile(np.unique(it_mix_scores_M), p), it_mix_scores_M) for p in eg_percs])

    plt.figure(figsize=(wTot, hTot))

    grid = plt.GridSpec(hTot, wTot, hspace=3, wspace=3)
    
    currW = 0
    # PLOT mask freq space 
    plt.subplot(grid[:, currW:(currW+w_freqs)], facecolor=facecolor)
    plot_title(title, color=color_title)
    plt.xticks(fontsize=fontsize-2)
    plt.yticks(fontsize=fontsize-2)
    imfCols = temd.plot_mask_freq_scores(it_mask_freqs, it_mix_scores, inds=eg_inds, ms=ms, color_='w', ms_=ms_)
    
    if show_variants:
        max_imfs = it_mask_freqs.shape[-1]
        fmin, fmax = 0, np.round(np.nanmax(it_mask_freqs), -1)
        for variant in variants[:1]:
            if variant == 'EMD':
                imfs = emd.sift.sift(X, max_imfs=max_imfs)
            elif variant == 'eEMD':
                imfs = emd.sift.ensemble_sift(X, max_imfs=max_imfs)
            elif variant == 'itEMD':
                from ccw_it_emd import it_emd
                imfs = it_emd(X, sample_rate, N_imf=max_imfs)[0]
            else:
                print('could not find variant')
                imfs = None
            
            if imfs is None:
                continue
            
            mixScore = temd.get_modeMixScore(imfs, np.arange(imfs.shape[1]))
            
            plt.hlines(mixScore, fmin, fmax, color='gray', zorder=1)
            #plt.text(0, mixScore+0.02, variant, color='gray')
            
                
    
    plt.xlabel('Mask freq. (Hz)', fontsize=fontsize)
    plt.ylabel('Mode mixing score (r)', fontsize=fontsize)
    plt.xscale('log')
    currW += w_freqs

    currW0 = np.copy(currW)
    currH = 0
    for egi, ind in enumerate(eg_inds):

        mask_freqs = it_mask_freqs[ind]
        sift_config = emd.sift.get_config('mask_sift')
        sift_config['mask_freqs'] = mask_freqs/sample_rate
        sift_config['max_imfs'] = len(mask_freqs)

        currW = np.copy(currW0)
        imfs = emd.sift.mask_sift(X, **sift_config)
        freqAx_psd, imfPSDs = temd.get_imfPSDs(imfs, sample_rate)

        plt.subplot(grid[currH:(currH+h_eg), currW:(currW+w_imfs)], facecolor=facecolor)
        plt.xticks(fontsize=fontsize-2)
        plt.yticks([])
        plot_emd(None, imfs, st=st, en=en, sr=sample_rate, imfCols=imfCols, lfpCol='w', lw_imfs=2)
        if egi == len(eg_inds)-1:
            plt.xlabel('Time (s)', fontsize=fontsize)
        currW += w_imfs
        plt.subplot(grid[currH:(currH+h_eg), currW:(currW+w_psd)], facecolor=facecolor)
        plt.xticks(fontsize=fontsize-2)
        plot_imfPSDs(freqAx_psd, imfPSDs, imfCols=imfCols)
        if egi == len(eg_inds)-1:
            plt.xlabel('Freq. (Hz)', fontsize=fontsize)
        currW += w_psd

        currH += h_eg

'''
    
def get_slidingBins(start, end, bin_step, bin_size):
    bins = np.column_stack([np.arange(start, end, bin_step), np.arange(start, end, bin_step)+bin_size])
    bins = bins[np.flatnonzero(bins[:, 1] < end)]
    return bins

def get_sliding_ledScores(bsnm, seshLab, rootFolder, behWin_mins=10, behInt_secs=30, last20=True):
    
    behWin = int(behWin_mins*60*sr)
    behInt = int(behInt_secs*sr)

    seshLen = get_seshLen4seshLab(bsnm, seshLab, rootFolder)
    ledTimes = getLEDtimes(bsnm, seshLab, rootFolder)

    if last20:
        behBins = np.column_stack([np.arange(int(seshLen-20*60*sr), seshLen, behInt), 
                                   np.arange(int(seshLen-20*60*sr), seshLen, behInt)+behWin])
    else:
        behBins = np.column_stack([np.arange(0, seshLen, behInt), 
                                   np.arange(0, seshLen, behInt)+behWin])
    
    behBins = behBins[np.flatnonzero(behBins[:, 1] < seshLen)]

    behScores = []
    for st, en in behBins:
        nSal, nCoc = [np.sum(np.digitize(ledTimes[zone][:,0], [st, en]) == 1) for zone in zones]
        behScores.append((nCoc - nSal) / behWin_mins)
    behScores = np.array(behScores)
    
    return behBins, behScores

def print_pStr(p, nComps=1, return_str=False):
    if p < 0.001/nComps:
        pStr = '***'
    elif p < 0.01/nComps:
        pStr = '**'
    elif p < 0.05/nComps:
        pStr = '*'
    else:
        pStr = 'n.s.'
    if return_str:
        return pStr
    print(pStr)
    
def print_pStr_diffEsts(diffs, return_str=False, null=0., nComps=1, two_sided=True):
    #if len(np.unique(diffs)) == 1 and diffs.mean() == null:
    #    print('oooooo')
    #    percs = [None, None, None]
    #    def sig_test(percVal):
    #        return False
    nSides = [1, 2][two_sided]
    if np.median(diffs) > null:
        
        percs = (100-np.array([95, 99, 99.9]))/(nSides*nComps)
        def sig_test(percVal):
            return percVal > null
    else:
        percs = 100-(100-np.array([95, 99, 99.9]))/(nSides*nComps)
        def sig_test(percVal):
            return percVal < null

    pStr = np.array(['n.s.', '*', '**', '***'])[np.concatenate([[True], [sig_test(np.percentile(diffs, p)) for p in percs]])][-1]
    if return_str:
        return pStr
    print(pStr)

def get_wvShiftSpec(spec, specFreqs):
    freqShifts = [int(len(sig.morlet(int(10.*sr/f),w=5,s=1,complete=True))/2) for f in specFreqs]
    spec = np.row_stack([np.roll(specf, shift) for specf, shift in zip(spec, freqShifts)])
    return spec

def get_psd2(X, sample_rate, window='hann'):

        def get_psd_(X, sample_rate, maxFreq, pointsPerHz):
            from scipy.signal import welch
            psdIndMax = int(maxFreq*pointsPerHz)
            psd = welch(X, fs=sample_rate, window=window, nperseg=int(sample_rate)*pointsPerHz)
            freqAx_psd, psd = psd[0][0:psdIndMax], psd[1][0:psdIndMax]
            return freqAx_psd, psd

        freqAx_psd, psd = [], []
        maxFreqs = [1, 10, 100, 200, 500]
        pointsPerHzs = [20, 4, 4, 1, 0.1]
        for maxFreq, pointsPerHz in zip(maxFreqs, pointsPerHzs):
            freqAx_psd_, psd_ = get_psd_(X, sample_rate, maxFreq, pointsPerHz)
            if not len(freqAx_psd):
                i = 0
            else:
                i = np.where(freqAx_psd_ > freqAx_psd[-1][-1])[0][0]
            freqAx_psd.append(freqAx_psd_[i:])
            psd.append(psd_[i:])
        freqAx_psd = np.concatenate(freqAx_psd)
        psd = np.concatenate(psd)

        return freqAx_psd, psd
    
def get_coh(x, y, sr, maxFreqs=[12, 40, 100, 200, 500], 
            npersegs=[4, 2, 1, 0.1, 0.05]):
    '''
    returns
    f, coh
    '''
    if any([arg is None for arg in [maxFreqs, npersegs]]):
        maxFreqs = [500]
        npersegs = [4]
    
    if len(maxFreqs) != len(npersegs):
        raise ValueError('len(maxFreqs) != len(npersegs)')
    npersegs = np.array(npersegs)*sr
    f, coh = [], []
    for i, maxFreq, nperseg in zip(range(len(maxFreqs)), maxFreqs, npersegs):
        f_, coh_ = sig.coherence(x, y, fs=sr, nperseg=nperseg)

        inds = np.digitize(f_, maxFreqs) == i
        f.append(f_[inds]) 
        coh.append(coh_[inds])
    f = np.concatenate(f)
    coh = np.concatenate(coh)
    return f, coh

    
def plot_psd(x, sr, color='k', lw=2, alpha=1):
    f, p = get_psd2(x, sr)
    plt.plot(f, p, color=color, lw=lw, alpha=alpha)
    plt.xscale('log')
    
    
    
def plot_imfPSDs(freqAx, imfPSDs, space=0.7, imfCols=None, alpha=0.5, zorder_add=0):
    if imfCols is None:
        imfCols = sb.color_palette('husl', imfPSDs.shape[0])[::-1]
    for imfi, psd in enumerate(imfPSDs):
        zorder = imfi + zorder_add
        y = psd/psd.max()-imfi*space 
        plt.plot(freqAx, y, color=imfCols[imfi], zorder=zorder)
        plt.fill_between(freqAx, np.zeros_like(y)-imfi*space, y, color=imfCols[imfi], zorder=zorder, alpha=alpha)
    plt.yticks([])
    plt.xscale('log')
    

def figplot_emd(x, imfs, sr, win=None, imfCols=None, lfpCol='k', lw_imfs=1, title=None):
    w_emd = 20
    w_psd = 5
    facecolor=None
    
    if win is None:
        st = 0
        en = imfs.shape[0]-1
    else:
        st, en = win

    if imfCols is None:
        imfCols = sb.color_palette('husl', imfs.shape[1])[::-1]
        
    
    hTot = 8
    wTot = w_emd + w_psd

    plt.figure(figsize=(wTot, hTot))
    grid = plt.GridSpec(hTot, wTot) #, hspace=3, wspace=3)


    plt.subplot(grid[:, :w_emd], facecolor=facecolor)
    if title is not None:
        plot_title(title)
    from ccw_tailorEMD_functions import plot_emd
    plot_emd(x, imfs, st=st, en=en, sr=sr, imfCols=imfCols, lfpCol=lfpCol, lw_imfs=lw_imfs)

    plt.subplot(grid[:, w_emd:], facecolor=facecolor)
    space = 0.7
    for imfi, imf in enumerate(imfs.T):
        f, p = get_psd2(imf, sr)
        y = p/p.max()-imfi*space 
        plt.plot(f, y, color=imfCols[imfi])
        plt.fill_between(f, np.zeros_like(y)-imfi*space, y, color=imfCols[imfi], zorder=imfi, alpha=0.5)
    plt.yticks([])
    plt.xscale('log')

def fix_cellInfo(cellInfo, bsnm):
    cellInfo['bsnm'] = np.repeat(bsnm, cellInfo.shape[0])
    cellInfo['ci'] = np.arange(cellInfo.shape[0])
    return cellInfo
    
def get_pinkGreen():
    pg = [sb.color_palette('PiYG', 4)[i] for i in [0, -1]]
    return pg

def shuffInds(inds, maxLen=None, add_jitter=True):
    intervals = np.diff(sorted(inds))
    shuffInds = np.cumsum(np.random.choice(intervals, len(inds)-1, replace=False))
    if not add_jitter:
        return shuffInds
    
    def f_():
        return 1
    
    if maxLen is None:
        a, b = 0, inds.min()
        f = f_
    else:
        a, b = sorted([inds.min(), maxLen-inds.max()])
        if np.sign(b) == -1:
            a, b = 0, inds.min()
            f = f_
        else:
            f = randSign
            
    def jitter():
        return np.random.choice(np.arange(a, b))*f()
    
    return shuffInds + jitter()


def get_trained_neuralNetwork(X_train, y_train, solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(5, 2), random_state=1, max_iter=2000):
    from sklearn.neural_network import MLPClassifier
    model = MLPClassifier(solver=solver, alpha=alpha, hidden_layer_sizes=hidden_layer_sizes, random_state=random_state, max_iter=max_iter)
    model.fit(X_train, y_train)
    return model
    
def get_z4perc(perc=95):
    return stats.norm.ppf(perc/100)   

def load_imfs4region(region, bsnm, seshLab, rootFolder, tetUseInds, tetLabs, imfStr='imf', emdStr='.opt_mEMD'):
    
    b = rootFolder+bsnm+'/'+bsnm+'_'+seshLab
    tetLab = getFocusRegTetLab(region, tetUseInds, tetLabs)
    
    return np.load(b+emdStr+'.'+imfStr+'.'+tetLab+'.npy')


def get_most_activWin(bsnm, dataFolder, nMins=5):
    with open(dataFolder+'activWins/'+bsnm+'.activWins_'+str(nMins)+'mins'+'.pkl', 'rb') as h:
        jj = pk.load(h)

    wins, winPropActiv = [jj[k] for k in ['wins', 'winPropActiv']]

    return wins[winPropActiv.argmax()]

def get_fibbonacci(n=12):
    seq = [1, 1]
    for i in range(n-2):
        seq.append(seq[-1]+seq[-2])
    return np.array(seq)

def get_inds4propAUC(Y, prop):

    a0 = Y.sum()
    maxi = Y.argmax()
    
    if Y[maxi] / a0 > prop:
        if maxi == 0:
            st, en = [maxi, maxi+2]
        elif maxi == len(Y)-1:
            st, en = [maxi-2, maxi]
        else:
            st, en = [maxi-1, maxi+1]
        return st, en
    
    
    if maxi in [0, len(Y)-1]:
        if maxi == 0:
            st, en = 0, 1
            finished_start, finished_end = True, False
        else:
            st, en = maxi-1, maxi
            finished_start, finished_end = False, True
    else:
        st, en = maxi-1, maxi+1
        finished_start, finished_end = False, False
    
    for _ in Y:
        if st == 0:
            finished_start = True
        elif en == len(Y)-1:
            finished_end = True
        
        if finished_start:
            en += 1
        elif finished_end:
            st -= 1
        elif Y[st-1] > Y[en+1]:
            st -= 1
        else:
            en += 1
        
        if Y[st:en].sum()/a0 > prop:
            break
    
    return st, en



barcodeFolder = '/Dupret_Lab/analysis/ccwilliams_analysis/figure_2/barcodes/imfComponents/detected_barcodes/'

from ccw_imfBarcode_functions import regImfUseInds_v1 ##### DEFAULT for below function ######
def get_barcodes(barcodeFolder, emdStr='.opt_mEMD', regImfUseInds=regImfUseInds_v1, nICs=30):
    ''' 
    returns:
    barcodeStr, allStates, intDimInds, stateCols
    '''

    from ccw_imfBarcode_functions import get_regImf_saveStr

    barcodeStr = get_regImf_saveStr(regImfUseInds, emdStr)+'_'+str(nICs)
    allStates = np.load(barcodeFolder+'all'+barcodeStr+'.detStates'+'.npy')
    intDimInds = np.load(barcodeFolder+'all'+barcodeStr+'.intDimInds'+'.npy')
    stateCols = sb.color_palette('husl', nICs)

    return barcodeStr, allStates, intDimInds, stateCols

def figplot_barcodes(barcodeFolder, data=None, regImfUseInds=None, nCols=6, w=5, h=5):
    if data is None:
        barcodeStr, allStates, intDimInds, stateCols = get_barcodes(barcodeFolder)
        regImfUseInds = regImfUseInds_v1
    else:
        barcodeStr, allStates, intDimInds, stateCols = data

    nICs = allStates.shape[0]
    nRows = int(np.ceil(nICs/nCols))
    
    plt.figure(figsize=(nCols*w, nRows*h))
    for sti, stateCol in enumerate(stateCols):
        plt.subplot(nRows, nCols, sti+1)
        plot_barcodeTit(sti, stateCols, fontsize=16)
        plot_barcode(allStates[sti], regImfUseInds, intDimInds)
        
        
def load_assemblyData(bsnm, seshLab, assemblyFolder):
    ''' 
    returns:
    assemblies, assStrenghs, assThreshs, cellInfo, cellCols
    '''

    rootFolder = get_rootFolder4bsnm(bsnm)

    assemblies = np.load(assemblyFolder+bsnm+'.cellAssemblies'+'.npy')
    assStrenghs = np.load(assemblyFolder+bsnm+'_'+seshLab+'.cellAssemblies'+'.strengths'+'.npy')
    assThreshs = np.load(assemblyFolder+bsnm+'_'+seshLab+'.cellAssemblies'+'.threshs'+'.npy')

    cellInfo = load_cellInfo(bsnm)
    cellCols = [[] for i in range(cellInfo.shape[0])]
    for region in regions:
        regTit, regCol = getTitColForReg(region)
        for i in get_spkInds(region, cellInfo):
            cellCols[i] = regCol

    return assemblies, assStrenghs, assThreshs, cellInfo, cellCols

def load_cellInfo(bsnm):
    rootFolder = get_rootFolder4bsnm(bsnm)
    with open(rootFolder+bsnm+'/'+bsnm+'.cellInfo'+'.pkl', 'rb') as h:
        cellInfo = pk.load(h)
    cellInfo['bsnm'] = np.repeat(bsnm, cellInfo.shape[0])
    cellInfo['ci'] = np.arange(cellInfo.shape[0])

    return cellInfo





def plot_assemblies(assemblies, cellCols, assAlphas=None, sd_thresh=2, sep=1.3):
    if assAlphas is None:
        assAlphas = np.repeat(1, assemblies.shape[0])

    for assi, ass in enumerate(assemblies):
        inds = np.where(np.abs(ass) > np.abs(ass).std()*sd_thresh)[0]
        cellCols_ = []
        for i, cellCol in enumerate(cellCols):
            if i in inds:
                cellCols_.append(cellCol)
            else:
                cellCols_.append('gray')
        plot_vec(ass, xshift=assi*sep, cols=cellCols_, alpha=assAlphas[assi])

def get_nMemberRegs(ass, cellInfo, sd_thresh=2):

    member_inds = np.where(np.abs(ass) < np.abs(ass).std()*sd_thresh)[0]
    reg_positiveMembers = []
    for region in regions:
        inds = get_spkInds(region, cellInfo)
        inds = np.intersect1d(inds, member_inds)
        if len(inds):
            reg_positiveMembers.append(np.sum(ass[inds] > 0) > len(inds)/2.)
        else:
            reg_positiveMembers.append(0)

    nMemberRegs = np.array(reg_positiveMembers).sum()

    return nMemberRegs



def get_optoTrigSpks(bsnms2run):

    trigLen_secs = 4
    trigLen = int(trigLen_secs*ccw.sr)

    cellInfo = []
    cellSpikesTrig0 = []
    splitCellSpikesTrig0 = [[] for spliti in range(2)]
    for bsnm in bsnms2run:

        try:
            rootFolder = ccw.get_rootFolder4bsnm(bsnm)
            tetRegs, tetLabs, nTets, seshLabs, nSessions = ccw.getTetAndSeshInfo(bsnm, rootFolder)
        except:
            print(bsnm, 'fail')

        cellSpikes = []
        triggers = []
        seshTot = 0
        load_fail = False
        for seshLab in seshLabs:
            b = rootFolder+bsnm+'/'+bsnm+'_'+seshLab

            try:
                cellSpikes_, cellInfo_ = vbf.loadActivityMatrix(b)
            except:
                load_fail = True
                continue
            cellSpikes.append(cellSpikes_)

            triggers_ = np.array(vbf.LoadIntervals(b, '.light_pulse'))[:,0]
            triggers_ += seshTot
            triggers.append(triggers_)

            seshTot += len(ccw.loadLFPs(bsnm, seshLab, rootFolder, tetRegs[0]))
        if load_fail:
            print(bsnm, 'fail')
            continue
        cellSpikes = np.column_stack(cellSpikes)
        triggers = np.concatenate(triggers)

        ##
        trigAx0, cellSpikesTrig0_ = ccw.triggeredAverage(cellSpikes, triggers, trigLen)

        for spliti, triggers_ in enumerate(ccw.splitX(triggers)):
            splitCellSpikesTrig0[spliti].append(ccw.triggeredAverage(cellSpikes, triggers_, trigLen)[1])


        cellSpikesTrig0.append(cellSpikesTrig0_)
        cellInfo_['bsnm']  = np.repeat(bsnm, cellInfo_.shape[0])
        cellInfo_['ci']  = np.arange(cellInfo_.shape[0])
        cellInfo.append(cellInfo_)

    cellSpikesTrig0 = np.row_stack(cellSpikesTrig0)
    for spliti in range(2):
        splitCellSpikesTrig0[spliti] = np.row_stack(splitCellSpikesTrig0[spliti])
    splitCellSpikesTrig0 = np.array(splitCellSpikesTrig0)
    cellInfo = pd.concat(cellInfo)

    return trigAx0, cellSpikesTrig0, splitCellSpikesTrig0, cellInfo




def plot_rankOverRep(scores, labels, labelCols, nShuffles=1000, sigThresh=(2.5, 97.5), alpha=1, 
                     lw=2, lw_sig=3, lw_insig=1, alpha_sig=1, alpha_insig=0.5, hi2lo=True):

    checks_failed = []
    checks_failed.append(len(scores) != len(labels))
    checks_failed.append(not all([l in list(labelCols.keys()) for l in np.unique(labels)]))

    if any(checks_failed):
        raise ValueError("scores and labels must be the same length. Also labelCols must contain entries for all labels")

    if hi2lo:
        sortInds = np.argsort(np.argsort(-scores))
    else:
        sortInds = np.argsort(np.argsort(scores))

    sortedLabels = ['' for i in range(len(scores))]

    for i, sortInd in enumerate(sortInds):
        sortedLabels[sortInd] = labels[i]
    sortedLabels = np.array(sortedLabels)

    n = len(scores)
    x = np.linspace(0, 100, n)

    for label in np.unique(labels):
        color = labelCols[label]
        predCounts = np.array(sortedLabels==label, dtype=int)

        predCumDist = np.cumsum(predCounts) / predCounts.sum()
        shuffDist = []
        for shi in range(nShuffles): 
            jjdist = np.zeros_like(predCumDist)
            for i in np.random.choice(np.arange(len(predCumDist)), predCounts.sum(dtype=int), replace=False):
                jjdist[i] = 1
            shuffDist.append(np.cumsum(jjdist)/jjdist.sum(dtype=float))
        shuffDist = np.row_stack(shuffDist)
        m = np.mean(shuffDist, axis=0)
        threshLo, threshHi = [np.percentile(shuffDist, p, axis=0) for p in sigThresh]

        sig = np.array([any([obi<loi, obi>hii]) for obi, loi, hii in zip(predCumDist, threshLo, threshHi)])
        plt.plot(x, predCumDist, color=color, lw=lw_insig, alpha=alpha_insig)
        
        boutTimes = getBoutTimes(sig)
        if boutTimes is not None:
            for st, en in boutTimes:
                if (en-st) > 2:
                    plt.plot(x[st:en], predCumDist[st:en], color=color, alpha=alpha_sig, lw=lw_sig)

    plt.plot(x, np.linspace(0, 1, n), color='gray', lw=lw, alpha=0.5, ls='--')    
    if hi2lo:
        plt.xticks(np.arange(0, 125, 25), np.arange(0, 125, 25)[::-1])
    else:
        plt.xticks(np.arange(0, 125, 25))



    
def get_rootFolder4aniID(aniID):
    aniIDs_ccw = np.concatenate([aniIDs_cl, aniIDs_mor, [getAniID(b) for b in np.concatenate([bsnmsall, bsnms_da, bsnms_pre, bsnms_beh])]])
    aniIDs_rr = np.array([getAniID(b) for b in bsnms_rr])
    aniIDs_db = np.array([getAniID(b) for b in bsnms_db])
    aniIDs_opto = np.array([getAniID(b) for b in bsnms_gaba])
    #print(aniIDs_ccw)
    
    if aniID in aniIDs_ccw:
        rootFolder = rootFolder_ccw
    elif aniID in aniIDs_rr:
        rootFolder = rootFolder_rr
    elif aniID in aniIDs_db:
        rootFolder = rootFolder_db
    elif aniID in aniIDs_opto:
        rootFolder = rootFolder_opto

    return rootFolder

def get_rootFolder4bsnm(bsnm):
    from rrMetaData import rr_bsnmL
    aniID = getAniID(bsnm)
    if aniID in np.concatenate([aniIDs_orig, aniIDs_cl, aniIDs_cpp, aniIDs_mor, aniIDs_chR_zone]) or bsnm in np.concatenate([bsnmsall, bsnms_da, bsnms_pre, bsnms_beh]):
        rootFolder = rootFolder_ccw
    elif bsnm in bsnms_rr:
        rootFolder = rootFolder_rr
    elif bsnm in bsnms_db:
        rootFolder = rootFolder_db
    elif bsnm in bsnms_gaba:
        rootFolder = rootFolder_opto
    elif bsnm[1:] in rr_bsnmL()[1]:
        rootFolder = rootFolder_rr
    elif bsnm in bsnms_steph_sucCPP:
        rootFolder = rootFolder_steph
    elif bsnm in bsnms_leila:
        print('to change when more data from Leila')
        rootFolder = '/Dupret_Lab/raw/lreddy_data/'
    elif bsnm in os.listdir(rootFolder_mc):
        rootFolder = rootFolder_mc
    elif bsnm in os.listdir(rootFolder_vl):
        rootFolder = rootFolder_vl
    else:
        print('no rootFolder found')

    return rootFolder

def get_seshLen4seshLab(bsnm, seshLab, rootFolder, sr0=20000., sr=sr):
    seshInfo  = vbf.LoadStages(rootFolder+bsnm+'/'+bsnm)
    seshLen = int(sr*[(seshInfo.iloc[si]['end_t'] - seshInfo.iloc[si]['start_t'])/sr0 \
                      for si in range(seshInfo.shape[0]) if seshInfo.iloc[si]['filebase'][-1] == seshLab][0])
    return seshLen


def load_cellFRs(cellInfo):
    binSize_ms=10
    cell_frs = []
    for i in range(cellInfo.shape[0]):
        bsnm, ci = [cellInfo.iloc[i][k] for k in ['bsnm', 'ci']]

        rootFolder = get_rootFolder4bsnm(bsnm)
        try:
            with open(rootFolder+bsnm+'/'+bsnm+'.binnedSpikesStats_'+str(binSize_ms)+'ms'+'.pkl', 'rb') as h:
                jj = pk.load(h)
            fr = jj['cellStats'][ci, 0] / (binSize_ms/1000.)
        except:
            fr = np.nan

        cell_frs.append(fr)
    cell_frs = np.array(cell_frs)
    #
    return cell_frs


def get_pulseTrigRaster(bsnm, ci, trigLen=int(1250*1.2), ext='.light_pulse'):

    rootFolder = get_rootFolder4bsnm(bsnm)

    tetRegs, tetLabs, seshLabs, _ = get_bsnmInfo(bsnm, rootFolder, multi_site=False) #tetIDs, tetLabs, nTets, seshLabs, nSessions = getTetAndSeshInfo(bsnm, rootFolder)

    pulses = []
    spikes = []
    seshTot = 0
    for seshLab in seshLabs:
        cellSpikes_, cellInfo_ = vbf.loadActivityMatrix(rootFolder+bsnm+'/'+bsnm+'_'+seshLab)
        seshLen = cellSpikes_.shape[1]
        pulses_ = np.array(LoadIntervals(rootFolder+bsnm+'/'+bsnm+'_'+seshLab, ext))[:,0]
        pulses.append(pulses_+seshTot)
        seshTot += seshLen

        spikes.append(cellSpikes_[ci, :])


    spikes = np.concatenate(spikes)
    pulses = np.concatenate(pulses)

    trigAx, trigMat = triggeredMatrix(spikes, pulses, trigLen)


    return trigAx, trigMat

def get_sidak_p(p, nComps):
    return 1 - (1-p)**(1/nComps)

def get_trigRast_args(binSize_ms=None, xlim=None, x_ext=None):
    if binSize_ms is None:
        binSize_ms=1000/1250.
    if xlim is None:
        xlim=(-20, 40)
    if x_ext is None:
        x_ext=50

    return binSize_ms, xlim, x_ext


def plot_trigRast_hist(trigAx, trigMat, pulseDur_ms=None, color='k', alpha=0.5, xlim=None, lw=2, y_shift=0,
                       alpha_pulse=0.4, color_pulse=pulseCol_chR, hist_fill=True, zorder=2, f=np.nansum):
    
    binSize_ms, xlim, x_ext = get_trigRast_args(xlim=xlim)

    trigTot = f(trigMat, axis=0)

    trigTot_ = binSpikes(trigTot, binSize_ms, np.sum)
    trigAx_ = binSpikes(trigAx, binSize_ms, np.median)
    xi, xj = [findNearestInd(x, trigAx_) for x in [xlim[0]-x_ext, xlim[1]+x_ext]]
    x, y = get_xy_4_histFill(trigAx_[xi:xj], trigTot_[xi:xj])
    
    y += y_shift

    if hist_fill:
        plt.fill_between(x, np.zeros_like(x)+y_shift, y, color=color, alpha=alpha, zorder=zorder)
    plt.plot(x, y, lw=lw, color=color, zorder=zorder)
    ymax = y.max()*1.1

    if pulseDur_ms is not None:
        plt.fill_between([0, pulseDur_ms], [y_shift, y_shift], [ymax, ymax], color=color_pulse, alpha=alpha_pulse, zorder=1)

    #trigAx_edges_ = cens2edges(trigAx_[xi:xj])
    #plt.xlim(trigAx_edges_[0], trigAx_edges_[-1])
    
    #plt.xlim([x[(findNearestInd(t, x))] for t in xlim])
    plt.xlim(xlim)
    #return x


def plot_trigRast_rast(trigAx, trigMat, cmap='Greys', xlim=None, show_zero=False, vmin=None, vmax=None, mask=None, maskCol='w'):
    binSize_ms, xlim, x_ext = get_trigRast_args(xlim=xlim)
    x1, x2 = [findNearestInd(x, trigAx) for x in xlim]
    x2 += 1
    if mask is not None:
        mask = mask[:, x1:x2]
    m = sb.heatmap(trigMat[:, x1:x2], vmin=vmin, vmax=vmax, cbar=False, cmap=cmap, mask=mask)
    m.set_facecolor(maskCol)
    fix_sbHeatmap()
    plt.xticks([])
    plt.yticks([])
    if show_zero:
        color, lw, alpha, linestyle = 'k', 2, 1, '--'
        zInd = findNearestInd(0, trigAx[x1:x2])

        plt.vlines(zInd, 0, trigMat.shape[0], color=color, lw=lw, alpha=alpha, linestyle=linestyle, zorder=3)




def get_frexAx4noiseRange(freqAx, fmin=48, fmax=52):
    freqAx_ = np.full_like(freqAx, np.nan)
    for i, fi in enumerate(freqAx):
        if any([fi<fmin, fi>fmax]):
            freqAx_[i] = fi
    return freqAx_



def load_stageZoneSpkCoup2dat(testStages, rootFolder, zoneCoup2Folder, zones2load=zones, emdStr='.opt_mEMD'):
    stageCellInfo = {}
    stageZoneCellRegImfPhaseMod = {}
    for testStage in testStages:
        stageCellInfo[testStage] = []
        stageZoneCellRegImfPhaseMod[testStage] =  {}
        for zone in zones2load:
            stageZoneCellRegImfPhaseMod[testStage][zone] = []

        for bsnm in stageBsnms[testStage]:
            seshLab  = getSeshLab4testStage(bsnm, testStage, rootFolder)
            path = zoneCoup2Folder+bsnm+'_'+seshLab+emdStr+'.cellZoneSpkCoup2'+'.pkl'
            if not os.path.exists(path):
                continue
            with open(path, 'rb') as h:
                zMod = pk.load(h)
            cellInfo_ = vbf.loadActivityMatrix(rootFolder+bsnm+'/'+bsnm+'_'+seshLab)[1]
            cellInfo_['bsnm'] = np.repeat(bsnm, cellInfo_.shape[0])
            cellInfo_['ci'] = np.arange(cellInfo_.shape[0])
            #
            stageCellInfo[testStage].append(cellInfo_)
            for zone in zones2load:
                stageZoneCellRegImfPhaseMod[testStage][zone].append(zMod[zone])

        try:
            stageCellInfo[testStage] = pd.concat(stageCellInfo[testStage])
            for zone in stageZoneCellRegImfPhaseMod[testStage]:
                stageZoneCellRegImfPhaseMod[testStage][zone] = np.concatenate(stageZoneCellRegImfPhaseMod[testStage][zone], axis=0)
        except:
            stageCellInfo[testStage] = None
            stageZoneCellRegImfPhaseMod[testStage] = None


    nPhaseBins = zMod['sal'].shape[-1]
    phaseCens = edges2cens(np.linspace(0, 2*np.pi, nPhaseBins+1))

    return phaseCens, stageCellInfo, stageZoneCellRegImfPhaseMod



def plot_spkCoup2ImfTit4stage(spkReg, coupReg, coupImfTit, testStage, phaseCens, stageCellInfo, stageZoneCellRegImfPhaseMod, 
                              smoothSDs=1, nCy=2, alpha=0.5, lw=2, add_xlab=True, add_ylab=True):

    phaseCens4plot = get_phaseCens4plot(phaseCens, nCy)
    spkTit, spkCol = getTitColForReg(spkReg)
    coupTit, coupCol = getTitColForReg(spkReg)

    if stageCellInfo[testStage] is None:
        plotNothing()
        return

    spkInds = get_spkInds(spkReg, stageCellInfo[testStage])

    ri = np.where(coupReg==regions)[0][0]
    imfi = get_imfi4imfTit(coupImfTit, coupReg)
    if imfi is None:
        plotNothing()
    else:
        for zi, zone in enumerate(zones):
            color = ['gray', spkCol][zi]


            obs = ([smooth(np.divide(histi, histi.sum(dtype=float)), smoothSDs) for histi in stageZoneCellRegImfPhaseMod[testStage][zone][spkInds, ri, imfi, :]])

            m = np.nanmean(obs, axis=0)
            se = stats.sem(obs, axis=0, nan_policy='omit')

            m, se = [np.concatenate([x]*nCy) for x in [m, se]]

            plt.plot(phaseCens4plot, m, color=color, lw=lw)
            plt.fill_between(phaseCens4plot, m-se, m+se, color=color, alpha=alpha)

    if add_ylab:
        plt.ylabel(spkTit+' spike prob.', color=spkCol, fontweight='bold')
    if add_xlab:
        plot_coupXlab(coupReg, coupImfTit, phaseCens4plot, nCy)
    else:
        plt.xticks(np.arange(0, 360*nCy+180, xdeg), fontsize=fontsize)
        plotWave(phaseCens4plot, nCy, col=coupCol, lw=lw)


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





def get_logBins(a=5, d=0, nBins=100):
    b=10**(a-d)
    bins_ = np.logspace(0.000001, a, nBins)
    bins_ /= b
    return bins_

def get_perc4val(a, x, nBinsPer1th=100, pExt=5.):
    perc0 = stats.percentileofscore(a, x)
    pRange = [(perc0-(pExt/2.)), (perc0+(pExt/2.))]
    ps = np.linspace(pRange[0], pRange[1], int(nBinsPer1th*pExt))
    ps = ps[np.where(np.logical_and(ps>=0, ps<=100))[0]]
    xs = np.array([np.percentile(a, p) for p in ps])
    perc = ps[findNearestInd(x, xs)]
    return perc



round_to_nsf = lambda x, n: x if x == 0 else round(x, -int(floor(log10(abs(x)))) + (n - 1))

def load_boutOvDat(bsnms2load, imfTit, burstFolder, minchain=5, nShuff=100):
    #
    bsnmOv = {}
    bsnmRegDurs_ms = {}
    bsnmInfo = {}
    for bsnm in bsnms2load:
        vbf.printTime(bsnm+'     ')

        ### LOAD ###
        regBouts, regBoutTimes = load_regBouts(bsnm, imfTit, burstFolder)
        if regBouts is None:
            continue

        regs2run = np.array(list(regBouts.keys()))
        regs2runov = np.concatenate([regs2run, ['ov']])

        bsnmInfo[bsnm] = regBouts[regs2run[0]]['info']

        ### WORK ###
        bsnmRegDurs_ms[bsnm] = {}
        for region in regs2runov:
            bsnmRegDurs_ms[bsnm][region] = 1000*(np.subtract(regBoutTimes[region][:,1], regBoutTimes[region][:,0])/sr)

        # OVERLAP
        shuffRegBoutTimes = []
        for shi in range(nShuff):
            regBoutTimes_shuff = {}
            for region in regBoutTimes:
                regBoutTimes_shuff[region] = get_boutTimes_shuff(regBoutTimes[region])

            shuffRegBoutTimes.append(regBoutTimes_shuff)

        bsnmOv[bsnm] = {}
        for ai, rega in enumerate(list(regBouts.keys())):
            for bi, regb in enumerate(list(regBouts.keys())):
                if ai == bi or ai > bi:
                    continue

                regs4ov = [rega, regb]
                k = rega+'~'+regb
                ov  = get_boutOv4regs(regBoutTimes, regs4ov)

                obs_shuff = np.array([get_boutOv4regs(shuffRegBoutTimesi, regs4ov) for shuffRegBoutTimesi in shuffRegBoutTimes])

                bsnmOv[bsnm][k] = ov/obs_shuff
        #
        obs = np.concatenate([1000*(np.subtract(shuffRegBoutTimes[shi]['ov'][:,1], shuffRegBoutTimes[shi]['ov'][:,0])/sr) for shi in range(nShuff)])
        bsnmRegDurs_ms[bsnm]['ovShuff'] = np.random.choice(obs, shuffRegBoutTimes[shi]['ov'].shape[0], replace=True)

    return bsnmOv, bsnmRegDurs_ms, bsnmInfo




def load_regBouts(bsnm, imfTit, burstFolder, minchain=5, emdStr=emdStr, seshLab=None):

    if seshLab is None:
        path = burstFolder+bsnm+emdStr+'.regBoutData'+'.minchain.'+str(minchain)+'.'+imfTit+'.pkl'
    else:
        path = burstFolder+bsnm+'_'+seshLab+emdStr+'.regBoutData'+'.minchain.'+str(minchain)+'.'+imfTit+'.pkl'

    if not os.path.exists(path):
        regBouts, regBoutTimes = None, None
        return regBouts, regBoutTimes

    with open(path, 'rb') as h:
        regBouts = pk.load(h)

    regs2run = np.array(list(regBouts.keys()))
    regs2runov = np.concatenate([regs2run, ['ov']])

    regBoutTimes = {}
    for region in regs2run:
        regBoutTimes[region] = []
        for chainBout in regBouts[region]['chainBouts']:
            st, en = get_win4chainBout(chainBout, regBouts[region]['cycles'])
            regBoutTimes[region].append([st, en])
        try:
            regBoutTimes[region] = np.row_stack(regBoutTimes[region])
        except:
            regBoutTimes[region] = np.array([])
            

    # overlap
    regTimes = get_regTimes(regBouts, regs2run)

    regOv = regTimes[regs2run[0]]
    for region in regs2run[1:]:
        regOv = np.intersect1d(regOv, regTimes[region])
    try:
        regBoutTimes['ov'] = get_boutWindows(np.diff(regOv)==1)
    except:
        regBoutTimes['ov'] = None



    return regBouts, regBoutTimes




def get_regTimes(regBouts, regs2run):
    regTimes = {}
    for region in regs2run:
        regTimes[region] = []
        for chainBout in regBouts[region]['chainBouts']:
            st, en = get_win4chainBout(chainBout, regBouts[region]['cycles'])
            regTimes[region].append(np.arange(st, en))
        regTimes[region] = np.concatenate(regTimes[region])
    return regTimes


def plot_regBoutDurHists(bsnmRegDurs_ms, bsnmInfo, rateNorm, plotOv=False, nHistBins=30, percLim=(0, 99), min2zero=False, lw=2, alpha=0.5, nospace=False, smoothSDs=None, cum=False):


    if plotOv:
        regs2run = ['ov', 'ovShuff']
    else:
        regs2run = [r for r in np.array(list(bsnmRegDurs_ms[list(bsnmRegDurs_ms.keys())[0]].keys())) if r in regions]
    obs_ = np.concatenate([np.concatenate([bsnmRegDurs_ms[bsnm][region] for bsnm in bsnmRegDurs_ms]) for region in regs2run])

    binMin, binMax = [np.percentile(obs_, p) for p in percLim]

    if min2zero:
        binMin = 0

    histBins_ = np.linspace(binMin, binMax, nHistBins+1)

    reg_plotInfo = {}
    for ri, region in enumerate(regs2run):

        aniCounts = []
        for bsnm in bsnmRegDurs_ms:
            cens, counts, _ = hist(bsnmRegDurs_ms[bsnm][region], histBins_)

            if rateNorm:
                counts = np.divide(counts, float(bsnmInfo[bsnm]['seshLen']))
                ylab = '% Occurance'
            else:
                counts = np.divide(counts, counts.sum(dtype=float))
                ylab = 'Proportion'

            if cum:
                counts = np.cumsum(counts)
                ylab = 'Cum. '+ylab

            counts = smooth(counts, smoothSDs)

            aniCounts.append(counts)
        aniCounts = np.row_stack(aniCounts)
        m, se = [f(aniCounts, axis=0) for f in mse]

        x, y = get_xy_4_histFill(cens, m)
        x, se1 = get_xy_4_histFill(cens, m-se)
        x, se2 = get_xy_4_histFill(cens, m+se)

        reg_plotInfo[region] = x, y, se1, se2

    if nospace:
        space=0
    else:
        space = np.max([reg_plotInfo[region][1].max()*0.8 for region in reg_plotInfo])

    ytickYs=[]
    ytickLabs=[]
    for ri, region in enumerate(regs2run):
        if region == 'ovShuff':
            regTit, regCol = 'Overlap (shuff.)', gray3
        elif region == 'ov':
            regTit, regCol = 'Overlap', pastRed
        else:
            regTit, regCol = getTitColForReg(region)

        yShift = -ri*space

        x, y, se1, se2 = reg_plotInfo[region]


        y_ = round_to_nsf(y.max(), 2)
        ytickYs.append(np.array([0, y_])-ri*space)
        ytickLabs.append([0, y_])



        plt.plot(x, y+yShift, color=regCol, lw=lw, label=regTit)
        plt.fill_between(x, se1+yShift, se2+yShift, color=regCol, alpha=alpha)

        #plt.plot(cens, m, color=regCol, lw=lw, label=regTit)
        #plt.fill_between(cens, m-se, m+se, color=regCol, alpha=alpha)
    #plt.yscale('log')
    if not nospace:
        ytickYs = np.concatenate(ytickYs)
        ytickLabs = np.concatenate(ytickLabs)
        plt.yticks(ytickYs, ytickLabs)

    plt.legend ()
    plt.xlabel('Bout duration (ms)')
    plt.ylabel(ylab)



def get_boutOv4regs(regBoutTimes, regs4ov):
    regTimes = {}
    for region in regs4ov:
        regTimes[region] = np.concatenate([np.arange(st, en) for st, en in regBoutTimes[region]])
    ovTimes = regTimes[regs4ov[0]]
    for region in regs4ov[1:]:
        ovTimes = np.intersect1d(ovTimes, regTimes[region])

    ov = float(len(ovTimes)) / np.max([regBoutTimes[region].max() for region in regBoutTimes])

    return ov



def get_boutTimes_shuff(boutTimes):
    """ currently small overlap issue ... """

    boutDurs = np.subtract(boutTimes[:,1], boutTimes[:,0])
    boutInts = np.concatenate([[boutTimes[0,0]], np.array([boutTimes[i+1,0] - boutTimes[i,1] \
                                                           for i in range(boutTimes.shape[0]-1)])])
    nBouts = len(boutDurs)

    rdmInds1, rdmInds2 = [np.random.choice(np.arange(nBouts), nBouts, replace=False) for i in range(2)]

    rdmSts = np.cumsum([boutInts[i] for i in rdmInds1])

    boutTimes_shuff = np.row_stack([[st, st+boutDurs[j]] for st, j in zip(rdmSts, rdmInds2)])

    return boutTimes_shuff



def plot_stateStrengths(timeAx, stateStrengths, stateCols, window, focusStis=[], focus_zlwalpha=(3, 3, 1), unfocus_zlwalpha=(2, 1, 0.3)):
    st, en = window
    for sti, stateCol in enumerate(stateCols):
        if sti in focusStis:
            zorder, lw, alpha = focus_zlwalpha
        else:
            zorder, lw, alpha = unfocus_zlwalpha
        plt.plot(timeAx[st:en], stateStrengths[sti, st:en], color=stateCol, lw=lw, alpha=alpha, zorder=zorder)
    plt.xlim(timeAx[st], timeAx[en])


def plot_regSpkRaster(cellSpikes, cellInfo, timeAx, window, lw_spk=3, cellH=0.5, regions=regions):
    st, en = window

    plt.grid(False)
    plt.xticks([])
    plt.yticks([])
    plt.xlim(timeAx[st], timeAx[en])

    cellY=0
    for spkReg in regions:
        spkTit, spkCol = getTitColForReg(spkReg)
        spkInds = get_spkInds(spkReg, cellInfo)
        for traini in cellSpikes[spkInds, :]: #regSpikes[spkReg][:, st:en]:
            spkis = np.where(traini)[0]
            try:
                plt.vlines(x=timeAx[spkis], ymax=cellY, ymin=(cellY-cellH), color=spkCol, lw=lw_spk)
            except:
                continue
            cellY -= cellH





def get_contBoutObs(bsnmBoutData, obType):
    """ 
    obType : 'dur_ms', 'strength', 'rate'  
    """

    bsnm = list(bsnmBoutData.keys())[0]
    nConts = len(bsnmBoutData[bsnm]['contInfo'])

    contBoutObs = []
    for conti in range(nConts):
        region, imfTit = bsnmBoutData[bsnm]['contInfo'][conti]
        regTit, regCol = getTitColForReg(region)
        color = regCol
        label = regTit+'_'+imfTit

        if obType == 'rate':
            obs = np.array([bsnmBoutData[bsnm]['contBoutData'][conti]['chainBouts'].shape[0] / float(bsnmBoutData[bsnm]['info']['seshLen']/eegHz) for bsnm in bsnmBoutData])
            ais = np.arange(len(obs))
        else:
            obs = []
            ais = []
            for ai, bsnm in enumerate(bsnmBoutData):
                for chainBout, chainAmp in zip(bsnmBoutData[bsnm]['contBoutData'][conti]['chainBouts'], 
                                               bsnmBoutData[bsnm]['contBoutData'][conti]['chainAmps']):
                    if obType == 'dur_ms':
                        st, en = get_win4chainBout(chainBout, bsnmBoutData[bsnm]['contBoutData'][conti]['cycles'])
                        obs.append(1000.*((en-st)/eegHz))
                    elif obType == 'strength':
                        obs.append(chainAmp)
                    ais.append(ai)
        obs = np.array(obs)
        ais = np.array(ais)
        contBoutObs.append({'obs' : obs, 'ais' : ais, 'label' : label, 'color' : color})
    return contBoutObs





def get_rdmWindows(l0, sigLen_s, sr, n=10):
    st, en = [0, l0-1]
    ll = int(sigLen_s*sr)
    timeAx = np.linspace(0, sigLen_s, ll)
    rdmSts = np.random.choice(np.arange(st, en-ll), n, replace=False)
    rdmWindows = np.column_stack([rdmSts, rdmSts+ll])
    return timeAx, rdmWindows



def load_regImfTitData4bsnm(bsnm, imfTit, rootFolder, seshLabs2use='all', 
                            emdLoadStrs=['imf', 'ia', 'if'], cycleLoadStrs=['cycles'], regions2load=regions):

    tetRegs, tetLabs, seshLabs, tetUseInds = get_bsnmInfo(bsnm, rootFolder)

    if seshLabs2use == 'all':
        if bsnm in bsnms:
            seshLabs2use = [getSeshLab4testStage(bsnm, testStage, rootFolder, mnfs=mnfs) \
                            for testStage in ['rec', 'ext', 'ren']]
        else:
            seshLabs2use = [seshLabs[i] for i in get_seshIndsSleep(bsnm, rootFolder)]

    regCoupData = {}
    for region in regions2load:
        imfi = get_imfi4imfTit(imfTit, region)
        if imfi is not None:
            regCoupData[region] = {}
            for k in np.concatenate([emdLoadStrs, cycleLoadStrs]):
                regCoupData[region][k] = []

    seshTot = 0
    for seshLab in seshLabs2use:
        for region in regCoupData:
            tetLab = getFocusRegTetLab(region, tetUseInds, tetLabs)
            imfi = get_imfi4imfTit(imfTit, region)

            for i, loadID in enumerate(emdLoadStrs):
                x = np.load(rootFolder+bsnm+'/'+bsnm+'_'+seshLab+'.opt_mEMD.'+\
                            loadID+'.'+tetLab+'.npy')[:, imfi]
                if not i:
                    seshLen = len(x)
                regCoupData[region][loadID].append(x)
            with open(rootFolder+bsnm+'/'+bsnm+'_'+seshLab+'.'+region+'.opt_mEMD.'+imfTit+'.cycleStats'+'.pkl', 'rb') as h:
                cycleStats = pk.load(h)
                for cycleLoadStr in cycleLoadStrs:
                    if cycleLoadStr in ['cycles', 'troughs']:
                        regCoupData[region][cycleLoadStr].append(cycleStats[cycleLoadStr]+seshTot)
                    else:
                        regCoupData[region][cycleLoadStr].append(cycleStats[cycleLoadStr])
            #
        seshTot += seshLen

    #
    for region in regCoupData:
        if regCoupData[region] is not None:
            for k in emdLoadStrs:
                regCoupData[region][k] = np.concatenate(regCoupData[region][k])
            for cycleLoadStr in cycleLoadStrs:
                if cycleLoadStr in ['cycles']: 
                    regCoupData[region][cycleLoadStr] = np.row_stack(regCoupData[region][cycleLoadStr])
                else:
                    regCoupData[region][cycleLoadStr] = np.concatenate(regCoupData[region][cycleLoadStr])
    #
    return regCoupData, seshLabs2use



def get_fInds(freqAx, fmin=0, fmax=150):
    fInds = np.where(np.logical_and(freqAx>fmin, freqAx<fmax))[0]
    return fInds


def get_vta_optoCellTypeClusteringData(bsnms2load, cellTypes, rootFolder, mntStr_=None, nPCs2use=3, use_ccw_da_optotagging=True):
    ''' mntStr depreciated '''

    from sklearn.decomposition import PCA
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
    cellTypesFolder =  '/Dupret_Lab/analysis/ccwilliams_analysis/figure_4/vta_cellTypes/'
    #mntStr_+'/ripple2/data/ccwilliams_analysis/figure_4/vta_cellTypes/'

    optoCellInfo = {}
    for cellType in cellTypes:
        if cellType == 'glu':
            # LOAD opto-tagged Glu cells
            with open(cellTypesFolder+'optoTagged_VTA_glu.cellInfo'+'.pkl', 'rb') as h:
                optoCellInfo['glu'] = pk.load(h)
        elif cellType == 'da':
            # LOAD opto-tagged DA cells
            if use_ccw_da_optotagging:
                with open(cellTypesFolder+'optoTagged_VTA_da.cellInfo'+'.pkl', 'rb') as h:
                    optoCellInfo['da'] = pk.load(h)
            else:
                with open(cellTypesFolder+'rr_optoTagged_VTA_da.cellInfo'+'.pkl', 'rb') as h:
                    optoCellInfo['da'] = pk.load(h)
        elif cellType == 'gaba':
            print('ccw to code')


    ## LOAD ISIs / ACs ##
    isiInt_ms='log'
    spkReg='vta'
    binSize_ms = 2
    cellInfo = []
    cellIsis = []

    for bsnm in bsnms2load:
        try:
            cellInfoi = vbf.loadActivityMatrix(rootFolder+bsnm+'/'+bsnm+'_'+'1', NeuronIds='all')[1]
            cellIsisi = np.load(cellTypesFolder+bsnm+'.cell_ISIs_'+str(isiInt_ms)+'ms'+'.npy')
            isiCens = np.load(cellTypesFolder+bsnm+'.cell_ISIs_'+str(isiInt_ms)+'ms'+'.intCens'+'.npy')
        except:
            print(bsnm)
            continue
        cellInfoi['bsnm'] = np.repeat(bsnm, cellInfoi.shape[0])
        cellInfoi['ci'] = np.arange(cellInfoi.shape[0])
        if bsnm in bsnms_da:
            spkInds = []
            for i in range(cellInfoi.shape[0]):
                ci = cellInfoi.iloc[i]['ci']
                match = np.intersect1d(np.where(bsnm==optoCellInfo['da']['bsnm'].ravel())[0], 
                                       np.where(ci==optoCellInfo['da']['ci'].ravel())[0])
                if len(match):
                    spkInds.append(i)
            spkInds = np.array(spkInds)
        else:
            spkInds = get_spkInds(spkReg, cellInfoi)
        if len(spkInds):
            cellInfo.append(cellInfoi.iloc[spkInds])
            cellIsis.append(cellIsisi[spkInds, :])

    cellInfo = pd.concat(cellInfo) 
    cellIsis = np.row_stack(cellIsis)
    #isiEdges = ccw.cens2edges(isiCens)

    ## remove nan ISIs ##
    keepInds = np.where(np.array([not np.isnan(s) for s in cellIsis.sum(axis=1)]))[0]
    cellInfo = cellInfo.iloc[keepInds]
    cellIsis = cellIsis[keepInds]


    cellTypeUseInds = {}
    for cti, cellType in enumerate(cellTypes):
        cellTypeUseInds[cellType] = []
        for i in range(cellInfo.shape[0]):
            bsnm, ci = [cellInfo.iloc[i][k] for k in ['bsnm', 'ci']]

            match = np.intersect1d(np.where(optoCellInfo[cellType]['bsnm'].ravel()==bsnm)[0], 
                                   np.where(optoCellInfo[cellType]['ci'].ravel()==ci)[0])
            if len(match):
                cellTypeUseInds[cellType].append(i)
        cellTypeUseInds[cellType] = np.array(cellTypeUseInds[cellType])

    optoLabs = []
    optoUseInds = []
    for cellType in cellTypes:
        optoLabs.append(np.repeat(cellType, len(cellTypeUseInds[cellType])))
        optoUseInds.append(cellTypeUseInds[cellType])
    optoLabs = np.concatenate(optoLabs)
    optoUseInds = np.concatenate(optoUseInds)

    pca = PCA(n_components=nPCs2use)
    pca.fit(cellIsis.T)

    optoEmb = LDA()
    optoEmb = optoEmb.fit(pca.components_[:, optoUseInds].T, optoLabs)
    optoEmb2plot = optoEmb.fit_transform(pca.components_[:, optoUseInds].T, optoLabs)

    optoPCs = pca.components_[:, optoUseInds].T
    cellPCs = pca.components_.T

    return optoEmb, optoEmb2plot, optoPCs, optoLabs, cellPCs, cellInfo

def plot_beeswarm(obs, color='gray', nHistBins=100, sh=1, ms=2, histBins_=None, symetrical=True,
                  alpha=0.5, marker='o', xShift=0, smoothSDs=2, expFactor=1, zorder=2, plot=True):
    #
    if histBins_ is None:
        d = np.nanstd(obs)*0.1
        histBins_ = np.linspace(np.nanmin(obs)-d, np.nanmax(obs)+d, nHistBins+1)

    cens, counts, _ = hist(obs, histBins_)
    counts = smooth(counts, smoothSDs)

    counts = np.divide(counts, counts.max())

    binInds = np.digitize(obs, histBins_)-1

    x = np.zeros_like(obs)
    for bi in np.unique(binInds):
        if bi >= len(counts):
            continue
        for i in np.where(binInds==bi)[0]:
            if symetrical:
                x[i] = counts[bi]*sh*np.random.random()*randSign()
            else:
                x[i] = counts[bi]*sh*np.random.random()

    x += xShift
    if plot:
        plt.plot(x, obs, marker, lw=0, ms=ms, color=color, alpha=alpha, zorder=zorder)
    else:
        return x, obs

def plot_beeswarm2(obs, color='gray', nHistBins=100, sh=1, ms=2, histBins_=None, rdm_x=True,
                  alpha=0.5, marker='o', xShift=0, smoothSDs=2, expFactor=1, zorder=2):
    #
    print('add ccw.___  to fix')
    if histBins_ is None:
        d = obs.std()*0.1
        histBins_ = np.linspace(obs.min()-d, obs.max()+d, nHistBins+1)

    cens, counts, _ = ccw.hist(obs, histBins_)
    counts = ccw.smooth(counts, smoothSDs)

    counts = np.divide(counts, counts.max())

    binInds = np.digitize(obs, histBins_)-1

    x = np.zeros_like(obs)
    for bi in np.unique(binInds):
        if bi >= len(counts):
            continue
        if rdm_x:
            for i in np.where(binInds==bi)[0]:
                x[i] = counts[bi]*sh*np.random.random()*randSign()
        else:
            xs = np.linspace(-counts[bi]*sh, counts[bi]*sh, np.sum(binInds==bi))
            for i, xi in zip(np.where(binInds==bi)[0], xs):
                x[i] = xi

    x += xShift
    plt.plot(x, obs, marker, lw=0, ms=ms, color=color, alpha=alpha, zorder=zorder)


def load_everything(bsnm, barcodeStr, barcodeFolder, emdStr='.opt_mEMD', seshLabs2use=None):
    rootFolder = get_rootFolder4bsnm(bsnm)
    tetRegs, tetLabs, seshLabs, tetUseInds = get_bsnmInfo(bsnm, rootFolder)
    
    if seshLabs2use is None:
        seshLabs2use = get_seshLabs2use(bsnm, rootFolder)


    info = {}
    stateStrengths = []
    cellSpikes = []
    regImfs = {}
    regAmps = {}
    for region in regions:
        regImfs[region] = []
        regAmps[region] = []

    for seshLab in seshLabs2use:
        try:
            stateStrengths.append(np.load(barcodeFolder+bsnm+'_'+seshLab+barcodeStr+'.stateStrengths'+'.npy'))
        except:
            stateStrengths = None
        try:
            cellSpikes_, cellInfo = vbf.loadActivityMatrix(rootFolder+bsnm+'/'+bsnm+'_'+seshLab)
            cellSpikes.append(cellSpikes_)
        except:
            cellSpikes, cellInfo = None, None

        for region in regions:
            tetLab = getFocusRegTetLab(region, tetUseInds, tetLabs)
            regImfs[region].append(np.load(rootFolder+bsnm+'/'+bsnm+'_'+seshLab+emdStr+'.imf.'+tetLab+'.npy').T)
            regAmps[region].append(np.load(rootFolder+bsnm+'/'+bsnm+'_'+seshLab+emdStr+'.ia.'+tetLab+'.npy').T)

    if stateStrengths is not None:
        stateStrengths = np.column_stack(stateStrengths) #stats.zscore(np.column_stack(stateStrengths), axis=1)
    if cellSpikes is not None:
        cellSpikes = np.column_stack(cellSpikes)
    for region in regions:
        regImfs[region] = np.column_stack(regImfs[region])
        regAmps[region] = np.column_stack(regAmps[region])

    info['seshLabs2use'] = seshLabs2use
    info['cellInfo'] = cellInfo
    info['emdStr'] = emdStr
    info['barcodeStr'] = barcodeStr
    info['seshLen'] = regImfs[regions[0]].shape[1]
    del cellInfo

    seshLen = regImfs[regions[0]].shape[1]
    timeAx_secs = np.linspace(0, seshLen/eegHz, seshLen)

    return regImfs, regAmps, stateStrengths, cellSpikes, info, timeAx_secs








def load_spkTrigMets(bsnms2load, mets2load, spkTrigSpecsFolder, rootFolder, speedControl=False, ext='_mpdSpikes', spkRegs=['vta']):
    loadStr = '_cellTrigX'+ext
    regCellInfo = {}
    regSpkTrig = {}
    for region in spkRegs:
        regCellInfo[region] = []
        regSpkTrig[region] = {}
        for met in mets2load:
            regSpkTrig[region][met] = []

    speedStr = ['', '.speedControl'][speedControl]
    loadStrs4mets = {'lfp' : '.trigRegLFP', 
                     'amps' : '.trigRegAmps',
                     'spec' : '.trigRegSpecs',
                     'spikes' : '.spikeTimesUsed'}

    for bsnm in bsnms2load:
        try:
            cellInfoi = vbf.loadActivityMatrix(rootFolder+bsnm+'/'+bsnm+'_'+'1', NeuronIds='all')[1]
            nCells = cellInfoi.shape[0]
            cellInfoi['bsnm'] = np.repeat(bsnm, nCells)
            cellInfoi['ci'] = np.arange(nCells)
            for ci, des in enumerate(cellInfoi['des'].ravel()):
                try:
                    spkReg = des[1:]
                    if spkReg in spkRegs:
                        for met in mets2load:

                            orig = np.load(spkTrigSpecsFolder+bsnm+loadStr+'.ci.'+str(ci)+\
                                           loadStrs4mets[met]+''+'.npy')

                            if speedControl:
                                control = np.load(spkTrigSpecsFolder+bsnm+loadStr+'.ci.'+str(ci)+\
                                           loadStrs4mets[met]+'.speedControl'+'.npy')
                            else:
                                control = np.zeros_like(orig)
                            spkTrigMet = np.subtract(orig, control)

                            if met in ['spec', 'amps']:
                                spkTrigMet = np.array([stats.zscore(spkTrigMet_ri, axis=1) \
                                                       for spkTrigMet_ri in spkTrigMet])
                            elif met == 'lfp':
                                spkTrigMet = stats.zscore(spkTrigMet, axis=1)

                            #                      
                            regSpkTrig[spkReg][met].append(spkTrigMet)
                        regCellInfo[region].append(cellInfoi.iloc[ci])
                except:
                    #if ci != 0:
                        #print(bsnm, ci, 'error half way thru cells...')
                    continue
        except:
            continue

    for spkReg in spkRegs:
        regCellInfo[region] = pd.DataFrame(regCellInfo[region])
        for met in mets2load:
            if met != 'spikes':
                regSpkTrig[spkReg][met] = np.array(regSpkTrig[spkReg][met])

    with open(spkTrigSpecsFolder+bsnms2load[0]+loadStr+'.info'+'.pkl', 'rb') as h:
        trigInfo = pk.load(h)
    trigAx = trigInfo['trigAx']
    specFreqs = trigInfo['specFreqs']
    return regSpkTrig, regCellInfo, trigAx, specFreqs



def get_bsnmInfo(bsnm, rootFolder, multi_site=True):
    ''' 
    returns:
    tetRegs, tetLabs, seshLabs, tetUseInds
    '''
    seshLabs = []
    base = rootFolder+bsnm+'/'+bsnm
    tetRegs = np.array(vbf.LoadTrodes(base)['desel'])
    tetLabs = []
    for i in range(len(tetRegs)):
        tetLabs.append(str(i+1))
    tetLabs = np.array(tetLabs)
        
    #
    os.chdir(rootFolder+bsnm)
    sessions = vbf.LoadPar(base)['sessions']
    for sesh in sessions:
        st_ind = sesh.find('_')+1
        seshLabs.append(sesh[st_ind:len(sesh)])
    seshLabs = np.array(seshLabs)
    #
    if multi_site:
        tetUseInds = getAniTetUseInds(bsnm)
    else:
        tetUseInds = None
    #
    return tetRegs, tetLabs, seshLabs, tetUseInds


def getTetAndSeshInfo(bsnm, rootFolder):
    print('depreciated')
    seshLabsOut = []
    base = rootFolder+bsnm+'/'+bsnm
    tetIDsOut = np.asarray(vbf.LoadTrodes(base)['desel'])
    tetLabsOut = []
    for i in range(len(tetIDsOut)):
        tetLabsOut.append(str(i+1))
    #
    os.chdir(rootFolder+bsnm)
    sessions = vbf.LoadPar(base)['sessions']
    for sesh in sessions:
        st_ind = sesh.find('_')+1
        seshLabsOut.append(sesh[st_ind:len(sesh)])
    #
    nTetsOut = len(tetIDsOut)
    nSessionsOut = len(seshLabsOut)
    #
    #
    return tetIDsOut, tetLabsOut, nTetsOut, seshLabsOut, nSessionsOut





def get_seshLabs2use(bsnm, rootFolder=None, mnfs=True):
    from ccw_coc_closedLoop_functions import aniIDs2use, aniIDDayBsnms
    from rrMetaData import rr_bsnmL
    if rootFolder is None:
        rootFolder = get_rootFolder4bsnm(bsnm)

    _, _, seshLabs, _ = get_bsnmInfo(bsnm, rootFolder)
    

    if bsnm in bsnms:
        seshLabs2use = [getSeshLab4testStage(bsnm, testStage, rootFolder, mnfs=mnfs) \
                        for testStage in ['rec', 'ext', 'ren']]
    elif bsnm in [aniIDDayBsnms[aniID]['post'] for aniID in aniIDs2use]:
        seshLabs2use = [getSeshLab4testStage(bsnm, testStage, rootFolder, mnfs=mnfs) \
                        for testStage in ['rec', 'ext']]
    elif bsnm in bsnms_pre:
        seshLabs2use = [getSeshLab4testStage(bsnm, 'pre', rootFolder, mnfs=mnfs)]
    elif bsnm[1:] in rr_bsnmL()[1]:
        desen = vbf.LoadStages(rootFolder+bsnm+'/'+bsnm)['desen'].ravel()
        seshLabs2use = [seshLab for seshLab, d in zip(seshLabs, desen) if any(['2P' in d, '2L' in d]) and 's b' not in d]
    else:
        seshLabs2use = [seshLabs[i] for i in get_seshIndsSleep(bsnm, rootFolder)]

    return seshLabs2use

def get_condSeshLabs(bsnm, return_seshIDs=False, drug='coc'):
    rootFolder = get_rootFolder4bsnm(bsnm)
    tetRegs, tetLabs, seshLabs, tetUseInds = get_bsnmInfo(bsnm, rootFolder)
    desen = vbf.LoadStages(rootFolder+bsnm+'/'+bsnm)['desen'].ravel()
    seshLabs2use = np.array([seshLab for seshLab, d in zip(seshLabs, desen) if 'train' in d.lower()])
    ds2use = np.array([d for seshLab, d in zip(seshLabs, desen) if 'train' in d.lower()])
    if return_seshIDs:
        seshIDs = np.empty(len(seshLabs2use), dtype=object)
        for zone in ['sal', drug]:
            condi = [i for i, d in enumerate(ds2use) if zone in d.lower()]
            if len(condi):
                i = condi[0]
                seshIDs[i] = zone
                ledPaired = ds2use[i].lower().split('led')[1][0]
                prei = [i for i, d in enumerate(ds2use) if 'led'+ledPaired+'off' in d.lower()]
                if len(prei):
                    i = prei[0]
                    seshIDs[i] = zone+'0'
        #seshIDs = np.array(['sal0', 'sal', 'coc0', 'coc'])
        return seshLabs2use, seshIDs 
    return seshLabs2use

def get_condSeshLabs_cpp(bsnm, drug='coc'):
    rootFolder = get_rootFolder4bsnm(bsnm)
    tetRegs, tetLabs, seshLabs, tetUseInds = get_bsnmInfo(bsnm, rootFolder)
    desen = vbf.LoadStages(rootFolder+bsnm+'/'+bsnm)['desen'].ravel()
    condSeshLabs = {}
    for cond in ['sal', drug]:
        condSeshLabs[cond] = seshLabs[[i for i, d in enumerate(desen) if cond in d.lower()][0]]
    return condSeshLabs

def get_cond_bsnms2run(aniID, return_dayIDs=False):
    from ccw_coc_closedLoop_functions import aniIDDayBsnms
    bsnms2run = []
    dayIDs = []
    for k in ['pre', 'cond2', 'cond3']:
        match = np.flatnonzero([k in d for d in list(aniIDDayBsnms[aniID].keys())])
        if not len(match):
            continue
        bsnms2run.append(aniIDDayBsnms[aniID][list(aniIDDayBsnms[aniID].keys())[match[0]]])
        dayIDs.append(list(aniIDDayBsnms[aniID].keys())[match[0]])
    bsnms2run = np.array(bsnms2run)
    dayIDs = np.array(dayIDs)
    if return_dayIDs:
        return bsnms2run, dayIDs 
    return bsnms2run



def get_seshLabs2use_colin_noLaser(bsnm, rootFolder):
    if bsnm not in bsnms_da:
        print(bsnm, 'not found in colin bsnms', 'check m--???')
        return
    stages = vbf.LoadStages(rootFolder+bsnm+'/'+bsnm)
    #tetIDs, tetLabs, nTets, seshLabs, nSessions = getTetAndSeshInfo(bsnm, rootFolder)
    _, _, seshLabs, _ = get_bsnmInfo(bsnm, rootFolder)
    seshLabs2use = np.array(seshLabs)[np.array([not any([d in ['stm_vtaph', 'lasert'] \
                                                         for d in deseni.split()]) \
                                                for i, deseni in enumerate(stages['desen'].ravel())])]
    return seshLabs2use

def get_trk_withinDist(trk, x, y, p=0.1):
    maxDistance = np.mean([trk['x'].max()-trk['x'].min(), trk['y'].max()-trk['y'].min()])*p
    within_distance = np.array([eucDist([x, y], [x_, y_]) for x_, y_ in zip(trk['x'], trk['y'])]) < maxDistance
    return within_distance

def split_list(list2split, nPerSplit=3):
    splitInds = np.arange(0, len(list2split), nPerSplit)
    listSplit = []
    for i in range(len(splitInds)):
        if i != len(splitInds) -1:
            listSplit.append(list2split[splitInds[i]:splitInds[i+1]])
        else:
            listSplit.append(list2split[splitInds[i]:])
    return listSplit


def get_ylim4bar(y, prop=0.2, ignore_zeros=True):
    if ignore_zeros:
        y[y==0.] = np.nan
    y_min, y_max = [f(y) for f in [np.nanmin, np.nanmax]]
    yRange = y_max-y_min
    return y_min - yRange*prop, y_max + yRange*prop



def get_winGradient(signal, winLen):
    signal_pad = np.concatenate([np.zeros(int(winLen/2)), signal, np.zeros(int(winLen/2))])
    return np.array([signal_pad[i+winLen-1]-signal_pad[i] for i in range(len(signal))])

def measure_switching(signals, winLen_ms=20, sr=1250.):
    winLen = int(sr*(winLen_ms/1000.))

    signals = stats.zscore(signals, axis=1)
    x2y = np.subtract(get_winGradient(signals[0], winLen), get_winGradient(signals[1], winLen))
    return stats.zscore(x2y)


def get_gaussian_filter(length=1250):
    counts, cens, _ = vbf.hist([np.random.randn(100) for i in range(100000)], 101)
    counts = np.concatenate([counts[:counts.argmax()+1], counts[:counts.argmax()][::-1]]) 
    counts = np.divide(counts, counts.sum(dtype=float))
    gaussian_filter = sig.resample(counts, length)
    return gaussian_filter/gaussian_filter.max()

def gaussian_spikeTrain(spikeTrain, filterLength_ms=50, sr=1250.):
    filter_length=int(sr*(filterLength_ms/1000.))
    if filter_length%2:
        filter_length += 1
    kernal = get_gaussian_filter(filter_length)
    train4calc = np.concatenate([np.zeros(int(filter_length/2)), spikeTrain, np.zeros(int(filter_length/2))])
    smooTrain = np.array([np.dot(kernal, train4calc[i:(i+filter_length)]) for i in range(len(spikeTrain))])
    return smooTrain

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y

def get_spkInds(spkReg, cellInfo):

    spkInds = np.array([i for i in range(cellInfo.shape[0]) if cellInfo.iloc[i]['des'][1:]==spkReg])

    if 'bsnm' in list(cellInfo.keys()):
        if len(np.intersect1d(bsnms_da, np.unique(cellInfo['bsnm']))): # if colins mice are present
            spkInds2add = np.where([bool(len(locateMotif(np.array([i for i in spkReg]), 
                                                         np.array([i for i in cellInfo.iloc[i]['des']])))) \
                                    for i in range(cellInfo.shape[0])])[0]
            spkInds = np.unique(np.concatenate([spkInds, spkInds2add]))

    spkInds = np.array(spkInds, dtype=int)
    return spkInds


def plot_isi(isi2plot, isiCens, color='gray', alpha=1):
    x, y = get_xy_4_histFill(isiCens, isi2plot)
    plt.fill_between(x, np.zeros_like(x), y, color=color, alpha=alpha)
    plt.xscale('log')


def plot_ac(acAx0, ac0, color, already_binned=False, alpha=1, binSize_ms=5, multiple_fRanges=False, multiple_fCols=False, multiple_fAlphas=False,
            fRange=None, col_fRange='gray', alpha_fRange=0.3, max2zero=False, norm2max=False,
            xlim=(-500, 500), fontsize=12, bar=True, lw=1, zorder=2):
    '''acAx0 in ms

    multiple_fRanges, multiple_fCols, multiple_fAlphas = list to zip
    if False;
    fRange: if you want to highlight a freq interval as an ac interval
    col_fRange: color for this single range

    '''
    if already_binned:
        ac = ac0
        acAx = acAx0
    else:
        ac0[findNearestInd(0, acAx0)] = 0.
        ac = binSpikes(ac0, binSize_ms, np.mean)
        acAx = binSpikes(acAx0, binSize_ms, np.median)

    if multiple_fRanges and multiple_fCols and multiple_fAlphas:
        for fRange, col_fRange, alpha_fRange in zip(multiple_fRanges, multiple_fCols, multiple_fAlphas):
            fRange_ts = np.array([1000./f for f in (fRange)])
            for s in [-1, 1]:
                  plt.fill_between(fRange_ts*s, [0]*2, [np.max(ac)*1.2]*2, color=col_fRange, alpha=alpha_fRange)
    else:
        if fRange is not None:
            fRange_ts = np.array([1000./f for f in (fRange)])
            for s in [-1, 1]:
                plt.fill_between(fRange_ts*s, [0]*2, [np.max(ac)*1.2]*2, color=col_fRange, alpha=alpha_fRange)

    if max2zero:
        ac[ac.argmax()] = 0.

    if norm2max:
        ac /= ac.max()
    if bar:
        plt.bar(acAx, ac, width=acAx[1]-acAx[0], color=color, alpha=alpha, zorder=zorder)
    else:
        plt.plot(acAx, ac, color=color, alpha=alpha, lw=lw, zorder=zorder)
    plt.xlim(xlim)





def get_nonLinOsc(freq=5, nonlin_deg=.25, nonlin_phi=-np.pi/4, seconds=5, sample_rate=eegHz):
    '''
    generate a non-sinosoidal oscillation
    nonlin_deg : change extent of deformation from sinusoidal shape [-1 to 1]
    nonlin_phi : change left-right skew of deformation [-pi to pi]
    '''
    import emd
    num_samples = int(sample_rate*seconds)
    timeAx = np.linspace(0, seconds, num_samples)
    osc = emd.utils.abreu2010( freq, nonlin_deg, nonlin_phi, sample_rate, seconds )
    return timeAx, osc



def LoadPar(b):
    '''LoadPar
    Parses block or stage level par file into <dict>-object.

    INPUT:
    - [b]:   <str> containing "block base" (= path + file-base; e.g., '/mnfs/swrd3/data/ddLab_merged/mdm96-2006-0124/mdm96-2006-0124')

    OUTPUT:
    - [par]: <dict> containing important info from par-file
               - 'nch' = <int> total number of recorded channels
               - 'ref_ch' = <int> channelID of 'reference'
               - 'ref_trode' = <int> tetrode-number of 'reference'
               - 'trode_ch' = <list> with for each tetrode a <list> with its channelIDs
               - 'sessions' = <list> with all session-names'''

    ## Read par-file, returns a <list> with each row converted into a <str>
    lines = open(b+'.par').readlines()

    ## Create an "anonymous" function to split a <str> into its constituent integers
    #to_ints = lambda x:map(int, x.split())
    #
    def to_ints(x):
            return x.split()
    #
    ## Extract total number of channels, number of tetrodes and channelID of "reference"
    nch, bits = to_ints(lines[0])
    num_trodes, ref_ch = np.array(lines[2].split()) #.astype(int)
    num_trodes = int(num_trodes)

    ## Create <list> with for each tetrode a <list> with its  channelIDs
    trode_ch = []
    # -loop over all tetrodes
    for l in lines[3:3+num_trodes]:
        l = to_ints(l)
        t = l[1:]
        # -check whether for this tetrode correct number of channels is listed in par-file
        assert int(l[0])==len(t), 'par error: n ch in trode'
        trode_ch.append(t)

    ## Find tetrode-number of "reference"
    for tetrodeIndex in range(1, num_trodes+1):
        if ref_ch in trode_ch[tetrodeIndex-1]:
            ref_trode = tetrodeIndex

    ## Create <list> with all session-names
    #sessions = map(str.strip, lines[4+num_trodes:])
    sessions = lines[4+int(num_trodes):]
    sessions = [s[:-1] for s in sessions]
    ##
    #
    # -check whehter correct number of sessions is listed in par-file
    #assert int(lines[3+num_trodes])==len(sessions), 'par error: n sessions'
    assert int(lines[3+int(num_trodes)])==len(sessions), 'par error: n sessions'
    ##

    ## Create <dict>-object and return it
    if 'ref_trode' in globals():
        par = {'nch':nch, 'ref_ch':ref_ch, 'ref_trode':ref_trode, 'trode_ch':trode_ch, 'sessions':sessions}
    else:
        par = {'nch':nch, 'ref_ch':ref_ch, 'trode_ch':trode_ch, 'sessions':sessions}
    return par



def MapLFPs(path, nch, dtype=np.int16,order='F'):
    '''Returns a 2D numpy <memmap>-object to a binary file, which is indexable as [channel, sample].

    INPUT:
    - [path]:              <str> containing full path to binary-file
    - [nch]:               <int> number of channels in binary file
    - [dtype]=np.int16:    <numpy-type> of binary data points'''

    ## Calculate the total number of data points in the provided binary-file
    size = os.path.getsize(path)
    size = size/np.dtype(dtype).itemsize

    ## Create and return the 2D memory-map object
    memMap = np.memmap(path, mode='r', dtype=dtype, order=order, shape=(nch, int(size/nch)))
    return memMap


def GetLFPis(base,ChAreas='1'):

    Info = LoadTrodes(base)
    Info['desel'] = Info['desel'].astype('str')

    if type(ChAreas) is str:
        if ChAreas is 'all':
            ChLabels = np.array(Info['desel'].index).astype(int)
            LFPis = np.array(Info['lfp_ch'][ChLabels]).astype(int)
        else:
            aux = np.where(Info['desel']==ChAreas)[0]
            ChLabels = np.array(Info['desel'].index[aux]).astype(int) # channel labels (index in panda sheet)
            LFPis = np.array(Info['lfp_ch'][ChLabels]).astype(int) # channel rows in lfp matrix
    elif type(ChAreas) is list:
        ChLabels = list()
        LFPis = list()
        for ChArea in ChAreas: 
            aux = np.where(Info['desel']==ChArea)[0]
            ChLabels = np.concatenate((ChLabels,np.array(Info['desel'].index[aux]))) # channel labels (index in panda sheet)
            LFPis = np.concatenate((LFPis,np.array(Info['lfp_ch'][ChLabels]).astype(int))) # channel rows in lfp matrix
        ChLabels = np.array(ChLabels).astype(int)
        LFPis = np.array(LFPis).astype(int)

    return LFPis,ChLabels


def LoadTrodes(b, par=None):
    '''Load "trodes" information (mostly from desel- and par-file).

    INPUT:
    - [b]:       <str> containing "block base"

    OUTPUT:
    - [trodes]:  <DataFrame>'''

    ## If not provided, load the par-file information
    if par is None:
        par = LoadPar(b)

    ## Read desel-file and store as <DataFrame> using "pandas"-package
    trodes = pd.read_csv(b+'.desel', header=None, names=['desel'])

    ## Add for each tetrode its # of channels, the channelID of its "main" channel and a list of all its channelIDs to the "trodes"-<DataFrame>
    trodes['n_tr_ch'] = [len(chs) for chs in par['trode_ch']]
    trodes['lfp_ch'] = [chs[0] for chs in par['trode_ch']]
    trodes['ch'] = par['trode_ch']

    ## Let the index of this <DatLoadStagesaFrame> start from 1 instead of 0
    trodes.index += 1

    ## Add "trode" as column, and set name of column-index to "trode" (NOTE: not sure why this is needed!?)
    #trodes['trode'] = trodes.index
    #trodes.index.set_names('trode', inplace=True)

    ## Return the "trodes"-information as <DataFrame>
    return trodes


def loadLFPs(bsnm, seshLab, rootFolder, region='all', conv2microv=1.95):

    par = LoadPar(rootFolder+bsnm+'/'+bsnm)
    lfpInds, _ = GetLFPis(rootFolder+bsnm+'/'+bsnm, ChAreas='all')
    lfps = MapLFPs(rootFolder+bsnm+'/'+bsnm+'_'+seshLab+'.eeg', int(par['nch']))
    lfps = lfps[np.array(lfpInds), :]
    lfps = np.multiply(lfps, conv2microv)
    if region=='all':
        return lfps
    elif region=='tetUseInds':
        #try:
        tetUseInds = getAniTetUseInds(bsnm)
        #except:
            #print('no tetUseInds found for', bsnm)
            #return lfps
        return lfps[tetUseInds, :]
    else:
        regInd = getAniTetUseInds(bsnm)[np.where(regions==region)[0][0]]
        return lfps[regInd, :]



def get_aligned_wvform(wvform, troughi=16, pad_val=0):

    """
    troughi : the trough will be this index
    pad_val: if original_troughi != troughi, the waveform will be shifted. the "gap" will be padded with this value 
    """    

    troughi0 = wvform.argmin()

    n2shift = troughi0 - troughi

    if np.sign(n2shift) == -1:
        wvform = np.concatenate([np.repeat(pad_val, abs(n2shift)), wvform[:n2shift]])    
    elif np.sign(n2shift) == 1:
        wvform = np.concatenate([wvform[abs(n2shift):], np.repeat(pad_val, abs(n2shift))])    
    return wvform

def get_wvform4cellInfoi(cellInfoi, norm_indiv=False, norm_mean=True, return_indiv=False, return_chWvs=False, rootFolder=None, seshLab=None, 
                         align_trough=False, filter_pulseSpks=False, inside_pulse=False, pulseStr='.laser_pulse', first_pulse_in_burst=[False, 500]):


    '''
    filter_pulseSpks - will read pulse file and filter spikes
    
    if first_pulse_in_burst[0], then pulses which commence {first_pulse_in_burst[1]-milliseconds} after a pulse will be ignored.   

    RETURNS:
    
        if return_chWvs:    wvAx_ms, chWvs
        
        if return_indiv:    wvAx_ms, indivWvs, indivTimes

        else:               wvAx_ms, wvform

    '''


    dtype='int16'
    #nCh=4
    nSamples=32

    bsnm, ci, tetLab, tetClu = [cellInfoi[k] for k in ['bsnm', 'ci', 'trode', 'trode_unit']]
    tetLab = str(tetLab)
    if rootFolder is None:
        rootFolder = get_rootFolder4bsnm(bsnm)

    b = rootFolder+bsnm+'/'+bsnm
    if seshLab is not None:
        b += '_'+seshLab

    clu = np.array(pd.read_csv(b+'.clu.'+tetLab, header=None)).ravel()[1:]
    res = np.array(pd.read_csv(b+'.res.'+tetLab, header=None)).ravel()
    cluInds = np.where(clu==tetClu)[0]
    #spkTimes_res = res[cluInds]

    datfname = b+'.spk.'+tetLab
    size = os.path.getsize(datfname)
    size = int(size/np.dtype(dtype).itemsize)

    #spks = np.memmap(datfname, mode='r', dtype=dtype, order='F', shape=(nCh, nSamples, len(res)))
    #else:
    array = np.fromfile(datfname, dtype=np.int16)
    spks = array.reshape(len(res), nSamples, -1).T[1:]
    spks = spks[:, :, cluInds]

    if filter_pulseSpks:
        spkTimes = res[cluInds]

        if seshLab is None:
            seshInfo = vbf.LoadStages(b)
            tetRegs, tetLabs, seshLabs, tetUseInds = get_bsnmInfo(bsnm, rootFolder, False)
            
            pulses = []
            for seshLab, t in zip(seshLabs, seshInfo['start_t'].ravel()):
                try:
                    pulses_ = pd.read_csv(b+'_'+seshLab+pulseStr, sep='\s+', header=None)
                except:
                    continue
                if pulses_.size:
                    pulses.append(pulses_+t)
            pulses = np.row_stack(pulses)
        
        else:
            pulses = np.array(pd.read_csv(b+pulseStr, sep='\s+', header=None))
        
        pulseTimes = np.concatenate([np.arange(st, en) for st, en in pulses])
        if inside_pulse:
            if first_pulse_in_burst[0]:
                ipis_ms = 1000*(np.diff(pulses[:,0])/20000.)
                pulses = pulses[1:, :][np.flatnonzero(ipis_ms > first_pulse_in_burst[1])]
                pulseTimes = np.concatenate([np.arange(st, en) for st, en in pulses])
                
            inds = np.where(np.in1d(spkTimes, pulseTimes))[0]
        else:
            inds = np.where(np.in1d(spkTimes, pulseTimes)==False)[0]
            
            

        spks = spks[:, :, inds]



    chIndivWvs = np.moveaxis(spks, 1, 2) # [tet X num. X time]

    chWvs = chIndivWvs.mean(axis=1)
    
    t = chWvs.shape[1]
    tWv = 1000 * t * 1/20000.
    wvAx_ms = np.linspace(-tWv, tWv, t)
    
    if return_chWvs:
        return wvAx_ms, chWvs
    
    chi = np.abs(chWvs).max(axis=1).argmax() # use tet with largest waveform amp

    indivWvs = chIndivWvs[chi]
    if norm_indiv:
        indivWvs = np.row_stack([np.divide(wvform, np.abs(wvform).max()) for wvform in indivWvs])



    if return_indiv:
        if align_trough:
            print('WARNING: gonna re-align all troughs for individual spikes (...ccw to think)')
            indivWvs = np.row_stack([get_aligned_wvform(wv) for wv in indivWvs])
        indivTimes =  res[cluInds]
        indivTimes = np.multiply(sr, np.divide(indivTimes, 20000.)).astype(int)
        return wvAx_ms, indivWvs, indivTimes
    else:
        wvform = indivWvs.mean(axis=0)
        if norm_mean:
            wvform /= np.abs(wvform).max()
        if align_trough:
            wvform = get_aligned_wvform(wvform)


        return wvAx_ms, wvform
    
'''
def get_wvform4cellInfoi(cellInfoi, norm_indiv=False, norm_mean=True, return_indiv=False, rootFolder=None, seshLab=None, 
                         align_trough=False):


    dtype='int16'
    #nCh=4
    nSamples=32

    bsnm, ci, tetLab, tetClu = [cellInfoi[k] for k in ['bsnm', 'ci', 'trode', 'trode_unit']]
    tetLab = str(tetLab)
    if rootFolder is None:
        rootFolder = get_rootFolder4bsnm(bsnm)

    b = rootFolder+bsnm+'/'+bsnm
    if seshLab is not None:
        b += '_'+seshLab

    clu = np.array(pd.read_csv(b+'.clu.'+tetLab, header=None)).ravel()[1:]
    res = np.array(pd.read_csv(b+'.res.'+tetLab, header=None)).ravel()
    cluInds = np.where(clu==tetClu)[0]
    #spkTimes_res = res[cluInds]

    datfname = b+'.spk.'+tetLab
    size = os.path.getsize(datfname)
    size = int(size/np.dtype(dtype).itemsize)

    #spks = np.memmap(datfname, mode='r', dtype=dtype, order='F', shape=(nCh, nSamples, len(res)))
    #else:
    array = np.fromfile(datfname, dtype=np.int16)
    spks = array.reshape(len(res), nSamples, -1).T[1:]
    
    #wvfms = wvfms[1:]

    #except:
        #print('too many spikes to load!')
        #return None, None
    spks = spks[:, :, cluInds]
    tetIndivWvs = np.moveaxis(spks, 1, 2) # [tet X num. X time]

    tetWvs = tetIndivWvs.mean(axis=1)
    teti = np.abs(tetWvs).max(axis=1).argmax() # use tet with largest waveform amp

    indivWvs = tetIndivWvs[teti]
    if norm_indiv:
        indivWvs = np.row_stack([np.divide(wvform, np.abs(wvform).max()) for wvform in indivWvs])


    t = indivWvs.shape[1]
    tWv = 1000 * t * 1/20000.
    wvAx_ms = np.linspace(-tWv, tWv, t)

    if return_indiv:
        if align_trough:
            print('WARNING: gonna re-align all troughs for individual spikes (...ccw to think)')
            indivWvs = np.row_stack([get_aligned_wvform(wv) for wv in indivWvs])
        return wvAx_ms, indivWvs
    else:
        wvform = indivWvs.mean(axis=0)
        if norm_mean:
            wvform /= np.abs(wvform).max()
        if align_trough:
            wvform = get_aligned_wvform(wvform)
            
        
        return wvAx_ms, wvform
'''


'''    
def get_wvform(bsnm, ci, rootFolder, norm=False):
    from tetProf_figs import getwv

    b = rootFolder+bsnm+'/'+bsnm

    units = vbf.LoadUnits(b)

    meanwv, stdwv, duration = getwv(b, ci, units)
    chi = np.abs(meanwv).max(axis=0).argmax()
    wvform = meanwv[:, chi]

    if norm:
        wvform /= np.abs(wvform).max()

    tWv = 1000 * len(wvform) * 1/20000.
    wvAx_ms = np.linspace(-tWv, tWv, len(wvform))

    return wvAx_ms, wvform


def get_wvform4cellInfoiOLD(cellInfoi, rootFolder, normWv=True):


    #returns 
    #wvAx_ms, wvform


    bsnm, tetLab, cluID = [cellInfoi[k] for k in ['bsnm', 'trode', 'trode_unit']]

    # clu, res
    clu = pd.read_csv(rootFolder+bsnm+'/'+bsnm+'.clu.'+str(tetLab), sep = ' ', header=None)[1:]
    clu.index = np.arange(0,len(clu))
    res = pd.read_csv(rootFolder+bsnm+'/'+bsnm+'.res.'+str(tetLab), sep = ' ', header=None)

    # get nWires
    #nWires_par = int(pd.read_table(rootFolder+bsnm+'/'+bsnm+'.par.'+str(tetLab), 
                                   #sep=' ', nrows=1, header = None)[1][0])
    nWires_par = int(pd.read_table(rootFolder+bsnm+'/'+bsnm+'.par', 
                                   sep=' ', nrows=1, header = None)[1][0])

    # read in binary
    array = np.fromfile(rootFolder+bsnm+'/'+bsnm+'.spk.'+str(tetLab), dtype=np.int16)


    #assert(len(array)/32/nWires_par == len(res))  ## PROBLEM WITH THIS 05.11.21 - ask vitor

    # number spikes X number samples per waveform X number channels
    n = int(len(array) / (len(res)*nWires_par)) #################### << fix - now gives downsampled wv vector.....
    waveforms = array.reshape(len(res),n,nWires_par)
    #waveforms = array.reshape(len(res),32,nWires_par)
    clu = clu[0].ravel()
    chWvForms = waveforms[np.where(cluID==clu)[0], :, :].mean(axis=0).T
    wvform = chWvForms[np.argmax(np.abs(chWvForms).max(axis=1))]
    if normWv:
        wvform /= np.abs(wvform).max()

    tWv = 1000 * len(wvform) * 1/20000.
    wvAx_ms = np.linspace(-tWv, tWv, len(wvform))


    return wvAx_ms, wvform
''' 


emd_variants = ['EMD', 'eEMD', 'ceEMD', 'itEMD', 'mEMD_zc', 'mEMD']
def run_emd(X, sample_rate, variant, args=None, max_imfs=9):
    # fix incorperate set maskfreqs - new arg(s)... **kwargs
    
    import emd
    
    if X is None:
        print(emd_variants)
        return
    if variant == 'EMD':
        imfs = emd.sift.sift(X, max_imfs=max_imfs)
    elif variant == 'eEMD':
        imfs = emd.sift.ensemble_sift(X, max_imfs=max_imfs)
    elif variant == 'ceEMD':
        print('ccw to sort!')
        # imfs, noise = emd.sift.complete_ensemble_sift(X, max_imfs=max_imfs)
        imfs = None
    elif variant == 'itEMD':
        from ccw_it_emd import it_emd
        imfs = it_emd(X, sample_rate, N_imf=max_imfs)[0]
    elif variant == 'mEMD_zc':
        imfs = emd.sift.mask_sift(X, max_imfs=max_imfs)
    elif variant == 'mEMD':
        mask_freqs = args[0]
        imfs = emd.sift.mask_sift(X, mask_freqs=mask_freqs/sample_rate, max_imfs=max_imfs)
    else:
        print('method not recognised')
        imfs = None

    return imfs

def get_scores4emd(Xs, sample_rate, variant, psd_func, args=None):
    # scores2compute=['m_corr', 'm_psd', 'c_psd']
    # fix - move to tmEMD.py?
    import ccw_tmEMD as temd
    
    kScores = {'m_corr' : [], 
               'm_psd' : []
              }
    
    X_imfPSDs = []
    for X in Xs:
        imfs = run_emd(X, sample_rate, variant, args)
        imfis_4_scoring = np.arange(imfs.shape[1])
        
        kScores['m_corr'].append(temd.get_modeMixScore(imfs, imfis_4_scoring))
        
        freqAx_psd, imfPSDs = temd.get_imfPSDs(imfs, sample_rate, psd_func)
        
        kScores['m_psd'].append(temd.get_modeMixScore_4_imfPSDs(imfPSDs, imfis_4_scoring, sample_rate))
        
        X_imfPSDs.append(imfPSDs)
    X_imfPSDs = np.array(X_imfPSDs)
    
    for k in kScores:
        kScores[k] = np.array(kScores[k])
    
    kScores['c_psd'] = temd.get_consistencyScores(freqAx_psd, X_imfPSDs, imfis_4_scoring)
    
    return kScores
        
        

    
    
    
def run_mEMD_JOSS(X, mask_freqs, sample_rate):
    import emd
    mask_amp=3
    mask_amp_mode='ratio_imf'
    sift_thresh=1e-08
    nphases=4
    imf_opts={}
    envelope_opts={}
    extrema_opts={}
    # mask_freq optimisation
    mask_args = {'mask_amp' : mask_amp,
                 'mask_amp_mode' : mask_amp_mode,
                 'sift_thresh' : sift_thresh, 
                 'nphases' : nphases,
                 'imf_opts' : imf_opts,
                 'envelope_opts' : envelope_opts,
                 'extrema_opts' : extrema_opts
                }

    sift_config = emd.sift.get_config('mask_sift')
    sift_config['mask_freqs'] = mask_freqs/sample_rate
    sift_config['max_imfs'] = len(mask_freqs)
    for k in mask_args:
        if mask_args[k] is not None:
            sift_config[k] = mask_args[k]

    imfs = emd.sift.mask_sift(X, **sift_config)
    
    return imfs


def get_bsnm_maxActivWins(bsnms, window, testStage, rootFolder, activSpeedMin=2):
    bsnmWins = {}
    for bsnm in bsnms:
        seshLab=getSeshLab4testStage(bsnm, testStage, rootFolder) #, mntStr=mntStr, mnfs=mnfs)
        lfp_ = loadLFPs(bsnm, seshLab, rootFolder, '1')
        seshLen = len(lfp_)
        speed = sig.resample(np.nan_to_num(getEEGtracking(
            rootFolder+bsnm+'/'+bsnm+'_'+seshLab)['speed'].ravel()), seshLen)
        #
        stInds = np.arange(0, seshLen-window, int(eegHz))
        winPropActiv = np.array([len(np.where(speed[st:(st+window)]>activSpeedMin)[0])/window \
                                 for st in stInds])
        st = stInds[winPropActiv.argmax()]
        en = st+window
        bsnmWins[bsnm] = np.array([st, en])
    return bsnmWins


def minuspi_2_2pi(phase_):
    phaseAdd = np.zeros_like(phase_)
    negInds = np.where(phase_<0.)[0]
    phaseAdd[negInds] = np.repeat(2*np.pi, len(negInds))
    phase_ += phaseAdd
    return phase_

def get_instPhase(coupTrace):
    #phase_ = np.angle(vbf.runHilbert(coupTrace))
    #phaseAdd=np.zeros_like(phase_)
    #negInds = np.where(phase_<0.)[0]
    #phaseAdd[negInds] = np.repeat(2*np.pi, len(negInds))
    #return phase_+phaseAdd
    phase = minuspi_2_2pi(np.angle(vbf.runHilbert(coupTrace)))
    return phase

    


#oldMid=(.99992311, .9976163, .26502115)

def get_spectral_freqCols(nMainFreqs, change_mid=[True, (0.94509804, 0.91372549, 0.09411765)]):
    if nMainFreqs==2:
        freqCols = [sb.color_palette('Spectral_r', 9)[i] for i in [0, -1]]
    elif nMainFreqs==3:
        freqCols = [sb.color_palette('Spectral_r', 9)[i] for i in [0, 4, -1]]
    elif nMainFreqs==4:
        freqCols = [sb.color_palette('Spectral_r', 10)[i] for i in [0, 3, 6, -1]]
    elif nMainFreqs==5:
        freqCols = [sb.color_palette('Spectral_r', 13)[i] for i in [0, 3, 6, 9, -1]]
    else:
        freqCols = sb.color_palette('Spectral_r', nMainFreqs)
    if change_mid[0] and nMainFreqs%2:
        freqCols[int(nMainFreqs/2)] = change_mid[1]
    return freqCols


def switch_colInPalette(palette, oldCol=(1.0, 1.0, 0.2), replaceCol=(0.94509804, 0.91372549, 0.09411765)):
    inds2replace = [i for i, rgb in enumerate(palette) if np.array_equal(rgb, oldCol)]
    for i in inds2replace:
        palette[i] = replaceCol
    return palette




def getColForStage(stage,
                   stageCols = {
                       'pre' : '#3dde8bff',
                       'post' : '#db3b6cff',
                       'rec' : '#db3b6cff',
                       'iext' : '#f5a2fff0',
                       'ext' : '#00aefff0',
                       'ren' : '#db9b3bff'}):
    return stageCols[stage]


def getMarkerForCellType(cellType):
    if cellType == 'p': markerOut = '^'
    elif cellType == 'b': markerOut = 'o'
    else: markerOut = '+'
    return markerOut


def get_2D_binned_heatmap(x, y, z, nBins_x=5, nBins_y=5, zfunc=np.median, emptyVal=np.nan):

    ''' each element in x and y is assigned a bin index. 
    a [nBins_x, nBins_y] matrix is returned, where each matrix element is the zfunc() of all z values which match these respective bins.
    if there are no entries to fill an x,y element, this value := emptyVal
    at the moment it uses percentile bins, but other binning methods can be implemented if necessary
    '''

    dimsxyz = np.row_stack([x, y, z])

    percs_ = {}
    for dim, nBins in zip(['x', 'y'], [nBins_x, nBins_y]):
        percs_[dim] = np.linspace(0, 100, nBins+1)


    binEdges = {}
    binInds = {}
    for di, dim in enumerate(['x', 'y']):
        binEdges[dim] = np.array([np.percentile(dimsxyz[di], p) for p in percs_[dim]])
        binInds[dim] = np.digitize(dimsxyz[di], binEdges[dim])-1

    hmap = np.zeros([nBins_x, nBins_y])
    for bx in range(nBins_x):
        indsx = np.where(binInds['x']==bx)[0]
        for by in range(nBins_y):
            indsy = np.where(binInds['y']==by)[0]
            if all([len(indsx), len(indsy)]):
                useInds = np.intersect1d(indsx, indsy)
                hmap[bx, by] = zfunc(dimsxyz[2][useInds])
            else:
                hmap[bx, by] = emptyVal
    xCens = edges2cens(binEdges['x'])
    yCens = edges2cens(binEdges['y'])

    return hmap, xCens, yCens



def get_2D_binned_probmap(x, y, nBins_x=5, nBins_y=5, zfunc=np.median, emptyVal=np.nan):

    ''' each element in x and y is assigned a bin index. 
    a [nBins_x, nBins_y] matrix is returned, where each matrix element is the zfunc() of all z values which match these respective bins.
    if there are no entries to fill an x,y element, this value := emptyVal
    at the moment it uses percentile bins, but other binning methods can be implemented if necessary
    '''

    dimsxy = np.row_stack([x, y])
    percs_ = {}
    for dim, nBins in zip(['x', 'y'], [nBins_x, nBins_y]):
        percs_[dim] = np.linspace(0, 100, nBins+1)

    binEdges = {}
    binInds = {}
    for di, dim in enumerate(['x', 'y']):
        binEdges[dim] = np.array([np.percentile(dimsxy[di], p) for p in percs_[dim]])
        binInds[dim] = np.digitize(dimsxy[di], binEdges[dim])-1

    hmap = np.zeros([nBins_x, nBins_y])
    for bx in range(nBins_x):
        indsx = np.where(binInds['x']==bx)[0]
        for by in range(nBins_y):
            indsy = np.where(binInds['y']==by)[0]
            if all([len(indsx), len(indsy)]):
                useInds = np.intersect1d(indsx, indsy)
                hmap[bx, by] = len(useInds)
    xCens = edges2cens(binEdges['x'])
    yCens = edges2cens(binEdges['y'])

    hmap /= hmap.sum()

    return hmap, xCens, yCens



def plotTitle(tit, col='black', fontsize=50, rotation=0, fontweight='bold', ha='center'):
    plotNothing()
    plt.text(0, 0, tit, color=col, fontsize=fontsize, ha=ha, 
             verticalalignment='center', rotation=rotation, fontweight=fontweight)


def plot_textInMid(text, color, fontweight='bold', fontsize=60):
    plt.text(0, 0, text, color=color, ha='center', va='center', fontweight=fontweight, fontsize=fontsize)
    plt.xlim(-0.1, 0.1)
    plt.ylim(-0.1, 0.1)
    plt.xticks([])
    plt.yticks([])

def getColMarkerForZone(zone):
    if zone == 'sal':
        colOut = salCol
        markerOut = 's'
    if zone == 'coc':
        colOut = cocCol
        markerOut = 'v'
    #
    return colOut, markerOut

def getTitleForZone(zone):
    if zone == 'sal':
        titOut = 'SAL'
    if zone == 'coc':
        titOut = 'COC'
    #
    return titOut

def get_simFunc(simMethod):
    ''' cos | euc | dot'''
    if simMethod=='cos':
        simFunc = cosine_similarity
    elif simMethod=='euc':
        simFunc = eucCloseness
    elif simMethod=='dot':
        simFunc = np.dot
    else:
        print('sim func not considered')
        simFunc = None        
    return simFunc

def cosine_similarity(x,y):
    #
    def square_rooted(x):
        #return round(sqrt(sum([a*a for a in x])),3)
        return sqrt(sum([a*a for a in x]))
    #
    numerator = sum(a*b for a,b in zip(x,y))
    denominator = square_rooted(x)*square_rooted(y)
    #return round(numerator/float(denominator),3)
    return numerator/float(denominator)

def plotNothing():
    plt.plot()
    plt.xticks([])
    plt.yticks([])

def splitX(X):
    """randomly splits elements in X to two sub-arrays"""
    X = np.array(X)
    split1Inds = np.random.choice(range(len(X)), int(len(X)/2.), replace=False)
    split2Inds =np.array([i for i in range(len(X)) if i not in split1Inds])
    split1Inds = np.array(sorted(split1Inds))
    #
    x1 = X[split1Inds]
    x2 = X[split2Inds]
    return x1, x2

'''
def splitX(X):
    split1Inds = np.random.choice(range(len(X)), len(X)/2, replace=False)
    split2Inds = [i for i in range(len(X)) if i not in split1Inds]
    #
    x1 = np.array([X[i] for i in split1Inds])
    x2 = np.array([X[i] for i in split2Inds])
    return x1, x2
'''


def randSign():
    return np.random.choice([-1, 1], 1)[0]





def plot_regTit(region, ext_pre='', ext_post='', fontsize=None, loc='left', 
                fontweight='bold', pad=10):
    regTit, regCol = getTitColForReg(region)
    plt.title(ext_pre+regTit+ext_post, color=regCol, loc=loc, fontweight=fontweight, fontsize=fontsize, pad=pad)



def get_tit4cellType(cellType):
    if cellType == 'glu':
        return 'Glu'
    elif cellType == 'da':
        return 'DA'
    elif cellType == 'gaba':
        return 'GABA'
    else:
        print('cellType not recognised')
        return cellType



def run_specICA(spec, nICs, dsHz=4):
    dsInt = int(eegHz/dsHz)
    dsInds = np.arange(np.random.choice(range(dsInt)), spec.shape[1], dsInt)
    from ccw_imfBarcode_functions import runICA
    ics_ = runICA(spec[:,dsInds].T, nICs)
    icCols = sb.color_palette('husl', nICs)[::-1]
    sortInds = np.argsort(np.argsort(np.abs(ics_).argmax(axis=1)))
    ics = np.zeros_like(ics_)
    for i, sortInd in enumerate(sortInds):
        if ics_[i].max() != np.abs(ics_[i]).max():
            ics[sortInd, :] = -ics_[i]
        else:
            ics[sortInd, :] = ics_[i]
    return ics


def flipSign(vec):
    """ flip based on mean of all neg vs mean of all pos """
    neg = np.mean(np.abs([i for i in vec.ravel() if i < 0]))
    pos = np.mean(np.abs([i for i in vec.ravel() if i > 0]))
    if neg > pos:
        return -vec
    else:
        return vec

def flipSign2(arr):
    ''' flip if sum is negative '''
    if np.sign(np.sum(arr)) > 0:
        return arr
    else:
        return -arr


def get_imfTit4imfi(imfi, region=None):
    if region is None:
        tit = tit = 'IMF'+str(imfi+1)
    elif imfi in [0, 7, 8]:
        if region in ['1']:
            if imfi == 7:
                tit = '2Hz'
            elif imfi == 8:
                tit = '1Hz'
            else:
                tit = 'IMF'+str(imfi+1)
        else:
            tit = 'IMF'+str(imfi+1)
    elif imfi==1:
        tit = 'fast-gamma'
    elif imfi==2:
        tit = 'mid-gamma'
    elif imfi==3:
        tit = 'slow-gamma'
    elif imfi==4:
        if region in ['vta']:
            tit = 'theta'
        else:
            tit = 'beta'
    elif imfi==5:
        if region in ['vta']:
            tit = '4Hz'
        else:
            tit = 'theta'
    elif imfi==6:
        if region in ['1']:
            tit = '4Hz'
        else:
            tit = '1Hz' #tit = 'IMF'+str(imfi+1)
    #
    return tit



def get_imfi4imfTit(imfTit, region, nImfs=9):
    try:
        imfiOut = np.where(np.array([get_imfTit4imfi(imfi, region) \
                                     for imfi in range(nImfs)])==imfTit)[0][0]
    except:
        imfiOut = None
    return imfiOut

def get_imfTits4reg(region, nImfs=9):
    imfTits = np.array([get_imfTit4imfi(imfi, region) for imfi in range(nImfs)])
    return imfTits

def get_imfi4freq(f, imfPSDs, freqAx_psd):
    imfi = findNearestInd(freqAx_psd[imfPSDs.argmax(axis=1)], f)
    return imfi


regImfUseInds_v1 = {'1' : np.array([4, 3, 2, 1]), 
                    'bla' : np.array([4, 3, 2, 1]), 
                    'nac' : np.array([4, 3, 2, 1]), 
                    'pfc' : np.array([4, 3, 2, 1]), 
                    'vta' : np.array([3, 2, 1])}

regImfUseInds_v2 = {'1' : np.array([4, 3, 2]),
                    'bla' : np.array([4, 3, 2]), 
                    'nac' : np.array([4, 3, 2]), 
                    'pfc' : np.array([4, 3, 2]), 
                    'vta' : np.array([3, 2])}

regImfUseInds_v3 = {'1' : np.array([5, 4, 3, 2, 1]), 
                    'bla' : np.array([5, 4, 3, 2, 1]), 
                    'nac' : np.array([5, 4, 3, 2, 1]), 
                    'pfc' : np.array([5, 4, 3, 2, 1]), 
                    'vta' : np.array([4, 3, 2, 1])}

regImfUseInds_v4 = {'1' : np.array([6, 5, 4, 3, 2, 1]), 
                    'bla' : np.array([5, 4, 3, 2, 1]), 
                    'nac' : np.array([5, 4, 3, 2, 1]), 
                    'pfc' : np.array([5, 4, 3, 2, 1]), 
                    'vta' : np.array([5, 4, 3, 2, 1])}


regImfUseInds_beta = {}
for region in regions:
    imfi = get_imfi4imfTit('beta', region)
    if imfi is not None:
        regImfUseInds_beta[region] = np.array([imfi])
    else:
        continue


def get_kSeshLabs_sucCPP(bsnm):
    rootFolder = get_rootFolder4bsnm(bsnm)
    desen = vbf.LoadStages(rootFolder+bsnm+'/'+bsnm)['desen'].ravel()
    tetRegs, tetLabs, seshLabs, tetUseInds = get_bsnmInfo(bsnm, rootFolder, multi_site=False)

    kSeshLabs = {}
    for k in ['pre', 'suc', 'wat', ' test']:
        kSeshLabs[k.split(' ')[-1]] = [seshLab for seshLab, d in zip(seshLabs, desen) if k in d][0]

    return kSeshLabs

bsnm_sucCPP_compartmentPaired = {}
bsnm_sucCPP_compartmentPaired['mdm220-170731'] = 'r'
bsnm_sucCPP_compartmentPaired['mdm227-171031'] = 'r'
bsnm_sucCPP_compartmentPaired['mdm168-161102'] = 'r'
bsnm_sucCPP_compartmentPaired['mdm227-171110'] = 'l'
bsnm_sucCPP_compartmentPaired['mdm169-161104'] = 'r'
bsnm_sucCPP_compartmentPaired['mst1095-170927'] = 'r'
bsnm_sucCPP_compartmentPaired['mdm195-170308'] = 'l'
bsnm_sucCPP_compartmentPaired['mdm227-171102'] = 'r'
bsnm_sucCPP_compartmentPaired['mst1095-170928'] = 'l'
bsnm_sucCPP_compartmentPaired['mst1095-170926'] = 'l'
bsnm_sucCPP_compartmentPaired['mdm169-161108'] = 'r'
bsnm_sucCPP_compartmentPaired['mdm227-171101'] = 'l'

def get_trk_sucCPP(bsnm, seshLab, rootFolder, n_xBins=201):
    trk = getEEGtracking(rootFolder+bsnm+'/'+bsnm+'_'+seshLab)

    xEdges = np.linspace(trk['x'].min(), trk['x'].max(), n_xBins+1)
    trk_x_binInds = np.digitize(trk['x'].ravel(), xEdges)-1

    x_yRange = np.row_stack([[f(trk['y'].iloc[np.flatnonzero(trk_x_binInds == bi)]) for f in [np.nanmin, np.nanmax]] for bi in range(n_xBins)])
    x_yRange = np.subtract(x_yRange[:, 1], x_yRange[:, 0])
    n_mask = int(n_xBins/3)
    x_yRange[:n_mask] = np.repeat(np.nan, n_mask)
    x_yRange[-n_mask:] = np.repeat(np.nan, n_mask)
    xCens = edges2cens(xEdges)

    x_ms = np.abs(np.diff(x_yRange))
    x_ms = np.nan_to_num(x_ms, -2)
    centre_i = int(n_xBins/2)
    shift = int(centre_i - x_ms.argmax())
    if np.sign(shift) == 1:
        st, en = centre_i-shift, centre_i+shift
    elif np.sign(shift) == -1:
        st, en = centre_i+shift, centre_i-shift

    if bsnm_sucCPP_compartmentPaired[bsnm] == 'l':
        def in_suc(x):
            return x < xCens[st]
        def in_wat(x):
            return x > xCens[en]
    elif bsnm_sucCPP_compartmentPaired[bsnm] == 'r':
        def in_wat(x):
            return x < xCens[st]
        def in_suc(x):
            return x > xCens[en]

    trk['suc'] = in_suc(trk['x'].ravel())
    trk['wat'] = in_wat(trk['x'].ravel())
    
    return trk

def get_kVisitTimes_sucCPP(trk):
    kVisitTimes = {}
    for k in ['suc', 'wat']:
        bouts = getBoutTimes(trk[k])
        kVisitTimes[k] = np.row_stack([[int(trk.iloc[i]['t']) for i in boutsi] for boutsi in bouts])
    return kVisitTimes


def get_trk_cpp(bsnm, seshLab, rootFolder, side_conditioned, n_xBins=201):
    '''
    side_conditioned : str  ['l' | 'r']
    '''
    trk = getEEGtracking(rootFolder+bsnm+'/'+bsnm+'_'+seshLab)

    xEdges = np.linspace(trk['x'].min(), trk['x'].max(), n_xBins+1)
    trk_x_binInds = np.digitize(trk['x'].ravel(), xEdges)-1

    x_yRange = np.row_stack([[f(trk['y'].iloc[np.flatnonzero(trk_x_binInds == bi)]) for f in [np.nanmin, np.nanmax]] for bi in range(n_xBins)])
    x_yRange = np.subtract(x_yRange[:, 1], x_yRange[:, 0])
    n_mask = int(n_xBins/3)
    x_yRange[:n_mask] = np.repeat(np.nan, n_mask)
    x_yRange[-n_mask:] = np.repeat(np.nan, n_mask)
    xCens = edges2cens(xEdges)

    x_ms = np.abs(np.diff(x_yRange))
    x_ms = np.nan_to_num(x_ms, -2)
    centre_i = int(n_xBins/2)
    shift = int(centre_i - x_ms.argmax())
    if np.sign(shift) == 1:
        st, en = centre_i-shift, centre_i+shift
    elif np.sign(shift) == -1:
        st, en = centre_i+shift, centre_i-shift

    if side_conditioned == 'l':
        def in_coc(x):
            return x < xCens[st]
        def in_sal(x):
            return x > xCens[en]
    elif side_conditioned == 'r':
        def in_sal(x):
            return x < xCens[st]
        def in_coc(x):
            return x > xCens[en]

    trk['coc'] = in_coc(trk['x'].ravel())
    trk['sal'] = in_sal(trk['x'].ravel())
    
    return trk







aniIDtetUseInds = {
    'mdm168' : {'1' : 1, 
                'nac' : 5,
                'nacCaudal' : 3},
    
    'mdm169' : {'1' : 2,
                'nac' : 10,
                'nacCaudal' : 3},
    
    'mdm195' : {'1' : 1, 
                'nac' : 2, 
                'nacCaudal' : 3},
    
    'mdm220' : {'1' : 7, 
                'nac' : 2, 
                'nacCaudal' : 3, 
                'nacLateral' : 0},
    
    'mdm227' : {'1' : 8, 
                'nac' : 4, 
                'u' : 3},
    
    'mst1095' : {'1' : 2, 
                 'nac' : 10}
}



    


def getAniTetUseInds(bsnm):
    ''' THIS FUNC TO BE DEPRECIATED ...'''
    aniID = getAniID(bsnm)
    #
    # 'pfc', 'nac', 'bla', '1', 'vta'
    #
    if aniID == 'hb17' or aniID == 'mhb17':
        tetUseInds_ = [8, 1, 3, 5, 6] #1, b, n, p, v
    elif aniID == 'rr02' or aniID == 'mrr02':
        tetUseInds_ = [0, 1, 3, 11, 7] #1, b, n, p, v
    elif aniID == 'rr04' or aniID == 'mrr04':
        tetUseInds_ = [8, 2, 4, 5, 6] #1, b, n, p, v
    elif aniID == 'ccw03' or aniID == 'mccw03':
        tetUseInds_ = [0, 2, 3, 5, 6] #1, b, n, p, v
    elif aniID == 'ccw05' or aniID == 'mccw05':
        tetUseInds_ = [0, 1, 4, 11, 6] #1, b, n, p, v
    elif aniID == 'am10' or aniID == 'mam10':
        tetUseInds_ = [0, 1, 4, 5, 6] #1, b, n, p, v
    elif aniID in ['am9', 'mam9', 'am5', 'mam5']:
        tetUseInds_ = [0, 1, 3, 5, 6] #1, b, n, p, v
    elif aniID == 'rr03' or aniID == 'mrr03':
        tetUseInds_ = [0, 3, 4, 12, 7] #1, b, n, p, v
    elif aniID == 'rr05' or aniID == 'mrr05':
        tetUseInds_ = [10, 0, 4, 12, 8] #1, b, n, p, v
    elif aniID == 'rr06' or aniID == 'mrr06':
        tetUseInds_ = [1, 3, 6, 11, 7] #1, b, n, p, v
    elif aniID in ['rr07', 'mrr07']:
        tetUseInds_ = [1, 3, 5, 13, 7] #1, b, n, p, v
    elif aniID in ['rr08', 'mrr08']:
        tetUseInds_ = [10, 2, 4, 12, 8] #1, b, n, p, v
    elif aniID in ['rr09', 'mrr09']:
        tetUseInds_ = [1, 2, 4, 11, 7] #1, b, n, p, v
    elif aniID in ['rr10', 'mrr10']:
        tetUseInds_ = [9, 0, 6, 13, 7] #1, b, n, p, v
    elif aniID in ['rr11', 'mrr11', 'rr11a', 'mrr11a', 'rr11b', 'mrr11b', 'rr11c', 'mrr11c']:
        tetUseInds_ = [1, 3, 4, 11, 7] #1, b, n, p, v
    elif aniID in ['rr12', 'mrr12', 'rr12a', 'mrr12a', 'rr12b', 'mrr12b', 'rr12c', 'mrr12c']:
        tetUseInds_ = [10, 0, 5, 13, 7] #1, b, n, p, v
    elif aniID in ['ccw06', 'mccw06']:
        tetUseInds_ = [9, 2, 4, 13, 8] #1, b, n, p, v
    elif aniID in ['ccw07', 'mccw07']:
        tetUseInds_ = [10, 2, 4, 11, 7] #1, b, n, p, v
    elif aniID in ['ccw09', 'mccw09']:
        tetUseInds_ = [9, 1, 3, 5, 8] #1, b, n, p, v
    elif aniID in ['ccw10', 'mccw10']:
        tetUseInds_ = [9, 3, 5, 11, 8] #1, b, n, p, v
    elif aniID in ['ccw11', 'mccw11']:
        tetUseInds_ = [9, 1, 5, 12, 7] #1, b, n, p, v
    elif aniID in ['ccw12', 'mccw12']:
        tetUseInds_ = [0, 1, 4, 10, 8] #1, b, n, p, v
    elif aniID in ['ccw13', 'mccw13']:
        tetUseInds_ = [0, 2, 6, 11, 7] #1, b, n, p, v
    elif aniID in ['ccw14', 'mccw14']:
        tetUseInds_ = [9, 3, 5, 12, 8] #1, b, n, p, v
    elif aniID in ['ccw16', 'mccw16', 'mccw16s']:
        tetUseInds_ = [0, 2, 5, 12, 8] #1, b, n, p, v
    elif aniID in ['ccw17', 'mccw17', 'mccw17s']:
        tetUseInds_ = [0, 2, 6, 12, 8] #1, b, n, p, v
    elif aniID in ['ccw18', 'mccw18', 'mccw18s']:
        tetUseInds_ = [9, 2, 5, 11, 8] #1, b, n, p, v
    elif aniID in ['ccw20', 'mccw20', 'mccw20s']:
        tetUseInds_ = [9, 2, 5, 11, 8] #1, b, n, p, v
    elif aniID in ['ccw23', 'mccw23']:
        tetUseInds_ = [0, 2, 6, 12, 8] #1, b, n, p, v
    elif aniID in ['mccw21']:
        tetUseInds_ = [0, 1, 5, 12, 8] #1, b, n, p, v
    elif aniID in ['mccw22', 'mccw22s']:
        tetUseInds_ = [10, 2, 5, 12, 8] #1, b, n, p, v
        
    
    else:
        print('no tetUseInds  found for', bsnm)
    #
    # !!!! swicth order to nre reg order as of 10.06.2021
    tetUseInds = []
    for ri, region in enumerate(regions):
        ind2match = np.where(region==regions_old)[0][0]
        tetUseInds.append(tetUseInds_[ind2match])
    tetUseInds = np.array(tetUseInds)
    return tetUseInds

def get_eegspeed(bsnm, seshLab, rootFolder):
    b = rootFolder+bsnm+'/'+bsnm+'_'+seshLab
    filename = b+'.eegspeed'
    if os.path.exists(filename+'.npy'):
        speed = np.load(filename+'.npy')
    else:
        trk = getEEGtracking(b)
        seshLen = get_seshLen4seshLab(bsnm, seshLab, rootFolder)
        speed = sig.resample(trk['speed'].ravel(), seshLen)
        os.chdir(rootFolder+bsnm)
        np.save(filename, speed)
    return speed
    
def getFocusRegTetLab(focusReg, tetUseInds, tetLabs, regions=regions):
    useInd = [i for i, reg in enumerate(regions) if reg == focusReg][0]
    return tetLabs[tetUseInds[useInd]]


def rotate(mat, nRots):
    '''rotates clock-wise'''
    def rotate90(m):
        return list(zip(*reversed(m)))
    if nRots != 0:
        for n in range(nRots):
            mat = rotate90(mat)
    return np.array(mat)

def get_pairwiseCompInds(N, maskWithin=True):
    compsOut = []
    for ai in range(N):
        try:
            compsOut.append(np.row_stack([[ai, aj] for aj in range(N)[::-1] if aj > ai]))
        except: continue
    #
    if not maskWithin:
        compsOut.append(np.column_stack([np.arange(N)]*2))

    return np.row_stack(compsOut)


def binInd2EEG(b, binLen_ms=25, eegHz=1250.):
    bHz = 1000./binLen_ms
    eegPerBin = eegHz/bHz
    return int(round(b*eegPerBin))


def eeg2binInd(e, binLen_ms=25, eegHz=1250.): # now in ccwBaseFunctions
    bHz = 1000./binLen_ms
    binPerEEG = bHz/eegHz
    return int(round(e*binPerEEG))



def eucDist(x, y):
    return np.sqrt(np.dot(x, x) - 2 * np.dot(x, y) + np.dot(y, y))

def eucCloseness(x, y):
    return eucDist(np.zeros_like(x), np.subtract(x, y))

def extractPatterns(data,dimrec,nPCs,nreplicates=50):

    if dimrec is 'ica':
        from sklearn.decomposition import FastICA
        model = FastICA(n_components=nPCs)
    if dimrec is 'nmf':
        from sklearn.decomposition import NMF
        model = NMF(n_components=nPCs, init='random', random_state=0)

    cs = np.zeros(nreplicates)
    candidatePatterns = [None]*nreplicates
    for runi in range(nreplicates):
        aux = model.fit(data).transform(data)
        c = np.corrcoef(pow(aux,2).T)
        cs[runi] = np.mean(np.abs(c-np.diag(np.diag(c))))
        candidatePatterns[runi] = np.copy(model)
    model = np.copy(candidatePatterns[np.argmin(cs)]).item()

    return model


def smooth(x, smoothSDs=1, axis=0):
    #
    if smoothSDs is None:
        return x
    return scipy.ndimage.filters.gaussian_filter1d(x, smoothSDs, axis=axis)






def get_cycleChunkInds(cycles, chSt, chLen):
    ''' cycles should be nCycles x 2 matrix '''
    cycles_ch = cycles - chSt
    cyStInds, cyEnInds = [np.where(np.logical_and(cycles_ch[:,i]>=0, 
                                                  cycles_ch[:,i]<chLen))[0] for i in [0,1]]
    cyUseInds = np.intersect1d(cyStInds, cyEnInds)
    return cycles_ch[cyUseInds, :]




def plot_pointLineSE(m, se, col, x=None, label=None, lw=2, marker='s', err=0.1, ls='-', ms=None, alpha=1):
    #
    if x is None:
        x = range(len(m))
    plt.plot(x, m, marker=marker, color=col, lw=lw, ls=ls, label=label, ms=ms, alpha=alpha)
    plt.vlines(x=x, ymin=m-se, ymax=m+se, color=col, lw=2, alpha=alpha, linestyles=ls)
    for i, xi in enumerate(x):
        for y in [m[i]-se[i], m[i]+se[i]]:
            plt.hlines(y=y, xmin=xi-err, xmax=xi+err, color=col, lw=lw, alpha=alpha, linestyles=ls)

def plot_pointLineSE_ax2(m, se, col, x=None, label=None, lw=2, marker='s', err=0.1, ms=None, alpha=1):
    #
    if x is None:
        x = range(len(m))
    ax2.plot(x, m, marker=marker, color=col, lw=lw, label=label, ms=ms, alpha=alpha)
    ax2.vlines(x=x, ymin=m-se, ymax=m+se, color=col, lw=2, alpha=alpha)
    for i, xi in enumerate(x):
        for y in [m[i]-se[i], m[i]+se[i]]:
            ax2.hlines(y=y, xmin=xi-err, xmax=xi+err, color=col, lw=lw, alpha=alpha)




def getNotNA(arr, returnInds=False):
    #
    indsOut = []
    arr = np.asarray(arr)
    #
    if returnInds is True:
        for ind, val in enumerate(arr):
            if ~np.isnan(val):
                indsOut.append(ind)
    #
    if returnInds is False:
        for val in arr:
            if ~np.isnan(val):
                indsOut.append(val)
    indsOut = np.array(indsOut)
    return indsOut


def makeUniform(arr):
    unifOut = []
    aMin = np.min(arr)
    aMax = np.max(arr)
    for val in enumerate(arr):
        unifOut.append(rdm.uniform(aMin, aMax))
    #
    return unifOut


def locateMotif(Motif, Arr):
    '''returns the start inds of a motif found in an array/list'''
    #
    #startInds = [ind for ind, val in enumerate(Arr) if val == Motif[0]]
    #startInds = [ind for ind in startInds if ind <= len(Arr)-len(Motif)]
    startInds = np.where(Arr==Motif[0])[0]
    startInds = startInds[np.where(startInds<=len(Arr)-len(Motif))]
    #
    motifStartIndsOut = []
    #
    for candidateInd in startInds:
        #
        if Arr[(candidateInd+(len(Motif)-1))] == Motif[-1]:
            candidateMotif = Arr[candidateInd:(candidateInd+(len(Motif)))]
            if np.array_equal(candidateMotif, Motif):
                motifStartIndsOut.append(candidateInd)
    #
    motifStartIndsOut = np.array(motifStartIndsOut)
    return motifStartIndsOut


def getCellTypeName(cellType):
    #
    if cellType == 'all':
        typeOut = 'All cell'
    else:
        reg = cellType[1:]
        if reg == '1':
            if cellType[0] == 'i':
                typeOut = 'Interneuron'
            else:
                typeOut = 'Pyramidal Cell'
        else:
            if cellType[0] == 'b':
                typeOut = 'Basket Cell'
            else:
                typeOut = 'Principal Cell'
    #
    return typeOut


def get_phase_est(b, seshLen):
    datfname, nCh, dtype = b+'.analogin', 3, np.int16
    size = os.path.getsize(datfname)
    size = int(size/np.dtype(dtype).itemsize)
    anCh = np.memmap(datfname, mode='r', dtype=dtype, order='F', shape=(nCh, int(size/nCh)))
    phase_est = sig.resample(anCh[1, :], seshLen)
    d = np.abs(phase_est).max() / np.pi
    phase_est /= d
    return phase_est

def convRadDeg(x, toDeg=True):
    #
    if toDeg:
        conv = (180./np.pi)
    else:
        conv = (np.pi/180.)
    return np.dot(x, conv) # i have no idea why this is dot??


def get_phaseCens4plot(phaseCens, nCy):
    return np.concatenate([np.array([convRadDeg(r) for r in phaseCens])+360.*cyi for cyi in range(nCy)])

def cyPhaseLags_02pi(pDiffs):
    ''' for phase lags between -2pi to 2pi, will shift them within 0 to pi '''
    edges = np.arange(-2*np.pi, 2*np.pi+1, np.pi)
    edges[-1] = edges[-1] + 0.0001
    binInds = np.digitize(pDiffs, edges)-1
    for bi, pShift in zip(range(4), [2*np.pi, np.pi, 0, -np.pi]):
        inds = np.flatnonzero(binInds == bi)
        pDiffs[inds] += np.repeat(pShift, len(inds))
    return pDiffs


def groupStrings(ids):
    '''returns ordered labs with corresponding inds from orig tetlist'''
    #
    orderedOut = sorted(ids)
    unique = set(orderedOut)
    indsOut = []
    #
    for ID in unique:
        matchInds = [ ind for ind, match, in enumerate(ids) if match == ID ]
        for matchInd in matchInds:
            indsOut.append(matchInd)
    #
    return orderedOut, indsOut


def LoadSpikeTimes(b, trode=None, MinCluId=2, res2eeg=(1250./20000)):
    # modified vitors so no error if no spikes
    t = '' if trode is None else '.'+str(int(trode))
    if os.stat(b+'.res'+t).st_size != 0:
        res = pd.read_csv(b+'.res'+t, header=None, squeeze=True).values
        clu = pd.read_csv(b+'.clu'+t, squeeze=True).values
        if MinCluId is not None:
            mask = clu >= MinCluId
            clu = clu[mask]
            res = res[mask]
        res = np.round(res*res2eeg).astype(int)
    #
    if os.stat(b+'.res'+t).st_size == 0:
        print('No spikes for tetrode '+t)
        res = []
        clu = []
    #
    return res,clu


def getCellFRs(spikes, samplingRate=1250.):
    #
    avFRs = []
    modeFRs = []
    for traini in spikes:
        if np.sum(traini) > 1:
            avFRs.append(float(np.sum(traini))/(np.size(spikes,1)/samplingRate))
            spInds = [i for i, s in enumerate(traini) if s == 1]
            #
            modeFRs.append(1./((stats.mode(np.diff(spInds))[0][0])/eegHz))
        else:
            avFRs.append(np.nan)
            modeFRs.append(np.nan)
    #
    return avFRs, modeFRs


def apply_minIndDist(inds, mpd):

    mpdInds = np.arange(0, inds.max()+mpd*2, mpd*2)
    return inds[np.unique([findNearestInd(mpdi, inds) for mpdi in mpdInds])]


'''
def apply_minIndDist(inds, mpd):

    ind_is2keep = []
    for i, ind in enumerate(inds):
        if all([i>0, i<(len(inds)-1)]):
            ind_is2keep.append(all([((inds[i+1] - ind) > mpd), ((ind - inds[i-1]) > mpd)]))
        elif i==0:
            ind_is2keep.append((inds[i+1] - ind) > mpd)
        else:
            ind_is2keep.append((ind - inds[i-1]) > mpd)
    ind_is2keep = np.array(ind_is2keep)

    stInds = list(np.array(locateMotif([True, False], ind_is2keep))+1)
    enInds = list(np.array(locateMotif([False, True], ind_is2keep)))
    #
    if all([len(stInds), len(enInds)]):
        if enInds[0] < stInds[0]:
            enInds = enInds[1:]
        if stInds[-1] > enInds[-1]:
            enInds.append(len(ind_is2keep)-1)
        #
        inds2sel1 = np.column_stack([stInds, enInds])

        mdpInds = inds[np.array([np.random.choice(np.arange(st, en)) for st, en in inds2sel1])]
    else:
        mdpInds = []

    inds = np.array(sorted(np.concatenate([mdpInds, inds[ind_is2keep]])))

    return inds
'''



def getLEDtimes(bsnm, seshLab, rootFolder, filtWins=False):
        #
        base = rootFolder+bsnm+'/'+bsnm+'_'+seshLab
        ledPaired = getLedPaired(bsnm)

        if ledPaired == 1:
            cocSuffix = led1Suffix
            salSuffix = led2Suffix
        else:
            cocSuffix = led2Suffix
            salSuffix = led1Suffix
        #
        cocLED = np.array(vbf.LoadIntervals(base, cocSuffix))
        salLED = np.array(vbf.LoadIntervals(base, salSuffix))
        #
        return {'sal' : salLED, 
                'coc' : cocLED }

def getLEDtimes2(bsnm, seshLab, rootFolder, filtWins=False):
        base = rootFolder+bsnm+'/'+bsnm+'_'+seshLab
        return {1 : np.array(vbf.LoadIntervals(base, led1Suffix)), 
                2 : np.array(vbf.LoadIntervals(base, led2Suffix)) }
    
    
def get_pseudo_ledTrigs_zoneStim(bsnm, seshLab, rootFolder, gapProp=0.005, nQuadrants=4):
    from ccw_coc_closedLoop_functions import bsnm_zoneStim_paired
    b = rootFolder+bsnm+'/'+bsnm+'_'+seshLab
    zonePaired = bsnm_zoneStim_paired[bsnm]
    zoneVisits, trk, zone_inFuncs, zoneBounds = get_zoneVisits_zoneStim(
        bsnm, seshLab, rootFolder, return_everything=True)
    zones = list(zoneBounds.keys())

    ledTrigs = {}
    for zone in zones:
        if zonePaired == 1:
            yWin = [[trk['y'].min()-1, zoneBounds[zone]], [zoneBounds[zone], trk['y'].max()+1]][zone == 'stim']
        elif zonePaired == 2:
            yWin = [[zoneBounds[zone], trk['y'].max()+1], [trk['y'].min()-1, zoneBounds[zone]]][zone == 'stim']
        xs = []
        for visit in zoneVisits[zone]:
            st, en = [findNearestInd(t, trk['t'].ravel()) for t in visit]
            xs.append(trk.iloc[st:en]['x'].ravel())
        xs = np.concatenate(xs)

        left, right = [np.percentile(xs, p) for p in [1, 99]]
        qEdges = np.concatenate([[xs.min()-1], [left + p*(right-left) for p in [0.25, 0.5, 0.75]], [xs.max()+1]])
        gap = gapProp*(right-left)

        qStatus = {}
        for qi in range(nQuadrants):
            if qi == 0:
                xWin = [qEdges[qi], qEdges[qi+1]-gap]
            elif qi != (nQuadrants-1):
                xWin = [qEdges[qi]+gap, qEdges[qi+1]-gap]
            else:
                xWin = [qEdges[qi], qEdges[qi+1]-gap]

            qStatus[qi] = np.add(np.array(np.digitize(trk['x'], xWin) == 1, dtype=int), 
                                 np.array(np.digitize(trk['y'], yWin) == 1, dtype=int)) == 2

        inds = np.array(sorted(np.concatenate([locateMotif([False, True], qStatus[qi]) for qi in qStatus])))
        ledTrigs[zone] = trk['t'].iloc[inds].ravel()
    return ledTrigs

    
def load_EEGtracking(b, smoothing=1,ext='.whl',trkHz=512.):
    '''Load position data (whl)'''
    trk = pd.read_csv(b+ext, sep='\s+', header=None).values
    trk[trk<=0] = np.nan
    if smoothing is not None:
        trk = scipy.ndimage.filters.gaussian_filter1d(trk, smoothing, axis=0)
    trk = pd.DataFrame(trk, columns=['x','y'])
    trk['x'] = vbf.linearInterpolate(trk['x'].ravel())
    trk['y'] = vbf.linearInterpolate(trk['y'].ravel())
    trk['t'] = convert20KHzToEEG(np.dot(np.arange(1, np.shape(trk)[0]+1), trkHz))
    trk['speed'] = np.concatenate([[np.nan], getInstSpeed(trk)])
    return trk






def getVisitEEGtimestamps(bsnm, seshLab, rootFolder, returnTrk=False, returnZone4y=False, mnfs=False,
                          mntStr=None):
    #
    #

    def getCoordsForOnTimes(onTimes, trkEEG):
        #
        trkEEG = trkEEG.dropna(axis=0)
        def getLocForTimeEEG(t):
            locInd = findNearestInd(t, trkEEG['t'].ravel())
            return np.array([trkEEG['x'].ravel()[locInd], trkEEG['y'].ravel()[locInd]])
        #
        coordsOut = []
        for ind, onTime in enumerate(onTimes):
            coords = getLocForTimeEEG(onTime)
            if not np.isnan(coords.sum()):
                coordsOut.append(coords)
        #
        return np.row_stack(coordsOut)
    ###################################################################################################
    def getLedZoneBoudary(coordsList, ledID, rankToUse=2): # currently takes the third most central LED stamp as boundary co-ordinate
        if bsnm in ['mccw05-180525'] and getSeshLab4testStage(bsnm, 'ext', rootFolder)==seshLab: # SPECIAL EXCEPTION AS BUGGY DATA....
            rankToUse=20
        if ledID == 1:
            boundOut = sorted(coordsList[:, 1], reverse=True)[rankToUse]
        if ledID == 2:
            boundOut = sorted(coordsList[:, 1], reverse=False)[rankToUse]
        return boundOut
    ###################################################################################################
    #
    b = rootFolder+bsnm+'/'+bsnm+'_'+seshLab
    ledPaired = getLedPaired(bsnm) #ledPaired = getLedPaired(bsnm, mntStr=mntStr)# this note should be here as of 07/04/2020
    ledTimes = getLEDtimes(bsnm, seshLab, rootFolder)
    eegTrack = load_EEGtracking(b)
    #
    salLoc = getCoordsForOnTimes(ledTimes['sal'][:, 0], eegTrack)
    cocLoc = getCoordsForOnTimes(ledTimes['coc'][:, 0], eegTrack)
    #
    if ledPaired == 1:
        salBound = getLedZoneBoudary(salLoc, 2)
        cocBound = getLedZoneBoudary(cocLoc, 1)
        def inCoc(y, cocBound=cocBound):
            return y < cocBound
        def inSal(y, salBound=salBound):
            return y > salBound
        def outOfZone(y, salBound=salBound, cocBound=cocBound):
            return cocBound <= y <= salBound 
    else:
        salBound = getLedZoneBoudary(salLoc, 1)
        cocBound = getLedZoneBoudary(cocLoc, 2)
        def inCoc(y, cocBound=cocBound):
            return y > cocBound
        def inSal(y, salBound=salBound):
            return y < salBound
        def outOfZone(y, salBound=salBound, cocBound=cocBound):
            return cocBound >= y >= salBound
    #
    #
    aniInSal = np.array([ inSal(locY) for locInd, locY in enumerate(eegTrack['y']) ])
    aniInCoc = np.array([ inCoc(locY) for locInd, locY in enumerate(eegTrack['y']) ])
    aniOut = np.array([ outOfZone(locY) for locInd, locY in enumerate(eegTrack['y']) ])
    #
    #
    def getVisitTimes(inStatus, eegTrack):
        #
        visitStartInds = list(np.array(locateMotif([False, True], inStatus))+1)
        visitEndInds = list(np.array(locateMotif([True, False], inStatus)))
        if not all([len(visitStartInds), len(visitEndInds)]):
            return np.array([])
        #
        if visitEndInds[0] < visitStartInds[0]:
            visitEndInds = visitEndInds[1:]
        if visitStartInds[-1] > visitEndInds[-1]:
            visitEndInds.append(np.size(eegTrack,0)-1)
        #
        visitStartTimes = [ int(round(eegTrack['t'][ind])) for ind in visitStartInds ]
        visitEndTimes = [ int(round(eegTrack['t'][ind])) for ind in visitEndInds ]
        #
        return np.column_stack([visitStartTimes, visitEndTimes])
    #
    #
    entriesOut = {'sal' : getVisitTimes(aniInSal, eegTrack),
                  'coc' : getVisitTimes(aniInCoc, eegTrack),
                  'out' : getVisitTimes(aniOut, eegTrack)}
    #
    if returnZone4y:
        return entriesOut, eegTrack, {'sal': inSal, 'coc': inCoc, 'out' : outOfZone}
    else:
        if returnTrk:
            return entriesOut, eegTrack
        else:
            return entriesOut

        
def get_zoneVisits_zoneStim(bsnm, seshLab, rootFolder, returnTrk=True, return_everything=False, perc=99):
    
    from ccw_coc_closedLoop_functions import bsnm_zoneStim_paired
    
    stages = vbf.LoadStages(rootFolder+bsnm+'/'+bsnm)
    
    if getAniID(bsnm) in ['mrr11', 'mrr12']:
        inds = np.array([i for i, d in enumerate(stages['desen']) if 'condtest' in d.lower()])
        stages = stages.iloc[inds]
        desen = np.array([['OFF', 'ON'][len(d.split('_')) > 1] for d in stages['desen'].ravel()])
        seshLabs = np.array([b.split('_')[-1] for b in stages['filebase'].ravel()])
    else:
        desen = np.array([d.split(' ')[-1] for d in stages['desen']])
        seshLabs = np.array([b.split('_')[-1] for b in stages['filebase']])

    seshInds = np.flatnonzero(desen == 'ON')
    stimYs = []
    for seshLab_ in seshLabs[seshInds]:
        b = rootFolder+bsnm+'/'+bsnm+'_'+seshLab_
        trk = load_EEGtracking(b)
        pulses = np.array(LoadIntervals(b, '.laser_pulse'))
        inds = np.array([findNearestInd(t, trk['t'].ravel()) for t in pulses[:,0]])
        stimYs.append(trk['y'].iloc[inds].ravel())
    stimYs = np.concatenate(stimYs)
    zonePaired = bsnm_zoneStim_paired[bsnm]
    top, bottom = [np.nanpercentile(trk['y'], p) for p in [perc, 100-perc]]
    if zonePaired == 1:
        yBound = np.percentile(stimYs, 100-perc)
        yBound_ = bottom + (top - yBound) 
        def in_stim(y):
            return y > yBound
        def in_null(y):
            return y < yBound_
    elif zonePaired == 2:
        yBound = np.percentile(stimYs, perc)
        yBound_ = top - (yBound - bottom)
        def in_stim(y):
            return y < yBound
        def in_null(y):
            return y > yBound_

    b = rootFolder+bsnm+'/'+bsnm+'_'+seshLab
    trk = load_EEGtracking(b)
    zoneVisits = {}
    for f, zone in zip([in_null, in_stim], ['null', 'stim']):
        zoneVisits[zone] = np.row_stack(
            [[int(trk['t'].iloc[st]), int(trk['t'].iloc[en])] \
             for st, en in get_boutWindows(f(trk['y'].ravel()))])
    if return_everything:
        zone_inFuncs = {'null' : in_null, 
                        'stim' : in_stim}
        zoneBounds = {'null' : yBound_, 
                      'stim' : yBound}
        
        return zoneVisits, trk, zone_inFuncs, zoneBounds
    
    if returnTrk:
        return zoneVisits, trk
    else:
        return zoneVisits
        
        
        
def get_zoneVisits_inWindow(zoneVisits_, seshSt, seshEn):
    seshTimes = np.arange(seshSt, seshEn)
    zoneVisits = {}
    for zone in zoneVisits_:
        zoneVisits[zone] = []
        vStatus = np.column_stack([np.in1d(visitsi, seshTimes) for visitsi in zoneVisits_[zone].T])
        for vStatusi, visit in zip(vStatus, zoneVisits_[zone]):
            if not vStatusi.sum():
                continue
            elif vStatusi.sum() == 2:
                zoneVisits[zone].append(visit)
            elif np.array_equal(vStatus, [False, True]):
                zoneVisits[zone].append(seshSt, visits[1])
            elif np.array_equal(vStatus, [True, False]):
                zoneVisits[zone].append(visits[0], seshEn)
        if len(zoneVisits[zone]):
            zoneVisits[zone] = np.row_stack(zoneVisits[zone])
        else:
            zoneVisits[zone] = np.array([])
    return zoneVisits


def get_ledBoundFunc(trk):
    xMin, xMax = [f(trk['x']) for f in [np.nanmin, np.nanmax]]
    xLen = xMax-xMin

    def getZoneBoundaryCat(t, trk, prop=0.35):
        trki = findNearestInd(t, trk['t'].ravel())
        x = trk.iloc[trki]['x']
        sideChecks = []
        sideChecks.append(xMin <= x < xMin+xLen*prop)
        sideChecks.append(xMax >= x > xMax-xLen*prop)
        if any(sideChecks):
            cat = 'side'
        else:
            cat = 'centre'
        return cat

    return getZoneBoundaryCat

def get_kVisitTimes_cpp(trk):
    print('!!! depreciated - use ccw.get_zoneVisits_cpp')
    return None
    

def get_zoneVisits_cpp(trk, zones=['sal', 'coc']):
    zoneVisits = {}
    for zone in zones:
        boutTimes = getBoutTimes(trk[zone])
        if trk[zone].iloc[0]:
            boutTimes = np.row_stack([[0, np.flatnonzero(trk[zone] == False)[0]-1], boutTimes])
        zoneVisits[zone] = np.row_stack([[int(trk['t'].iloc[st]), int(trk['t'].iloc[en])] for st, en in boutTimes])
    return zoneVisits


def get_behData(bsnm, testStage, rootFolder, runSpeed=True):
    #
    seshLab = getSeshLab4testStage(bsnm, testStage, rootFolder)
    seshLen = loadLFPs(bsnm, seshLab, rootFolder).shape[1]
    if testStage == 'iext':
        seshSt = 0
        seshEn = int(20*60*eegHz)
    else:
        seshSt = seshLen-int(20*60*eegHz)
        seshEn = seshLen-1
    ledTimes = getLEDtimes(bsnm, seshLab, rootFolder)
    visits, trk = getVisitEEGtimestamps(bsnm, seshLab, rootFolder, returnTrk=True)
    if runSpeed:
        speed = sig.resample(np.nan_to_num(trk['speed']), seshLen)
    else:
        speed = None
    #
    return ledTimes, visits, trk, speed, seshSt, seshEn, seshLen



def findNearest(val, array, returnVal=True):
    '''returnVal=True to return value rather than index'''
    array = np.asarray(array)
    ind = (np.abs(array - val)).argmin()
    if returnVal is False:
        return ind
    #
    return array[ind]

def findNearestInd(val, array):
    array = np.array(array)
    d = np.abs(array - val)
    np.nan_to_num(d, False, np.nanmax(d))
    ind = d.argmin()
    return ind

def find_smallest_positive_index(arr):
    arr = np.array(arr)
    positive_indices = np.where(arr > 0)[0]
    
    if len(positive_indices) == 0:
        return -1  # Return -1 if no positive numbers found
    
    min_positive_index = np.argmin(arr[positive_indices])
    return positive_indices[min_positive_index]

def find_smallest_positive_indices(arr):
    positive_mask = arr > 0
    positive_indices = np.argmax(positive_mask, axis=0)
    positive_indices[~positive_mask.any(axis=0)] = -1
    return positive_indices

def rolling_mean(x, N):
    cumsum = np.nancumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)


def get_rollingY_nD(arr, N, axis=0, f=np.mean):
    '''
    rolling {f} for a ND array
    returns:
    rollingX, rollingY
    '''
    if axis != 0:
        raise ValueError('ccw to code')
    bins = get_slidingBins(0, arr.shape[axis], 1, N)
    rollingX = bins.mean(axis=1)
    rollingY = np.array([f(arr[st:en], axis=0) for st, en in bins])
    return rollingX, rollingY


def effSize(x, y):
    '''x-y'''
    return (np.mean(x)-np.mean(y))/np.std(np.concatenate([x, y]))






def saveAsPickle(name, obj):
    with open(name+'.pkl', 'wb') as handle:
        pk.dump(obj, handle, pk.HIGHEST_PROTOCOL)



def saveAsJSON(obj, name):
    with open(name, 'w') as fileOut:
        json.dump(obj, fileOut)

def getLedPaired(bsnm, path='/Dupret_Lab/analysis/ccwilliams_analysis/figure_1/aniInfo/ledsPaired'):
    #
    from ccw_coc_closedLoop_functions import aniID_mor_dayBsnms, aniIDs2use_daGlu, aniIDDayBsnms
    if bsnm in np.concatenate([[aniID_mor_dayBsnms[aniID][k] for k in aniID_mor_dayBsnms[aniID]] for aniID in aniID_mor_dayBsnms]):
        path = '/Dupret_Lab/analysis/ccwilliams_analysis/figure_1/aniInfo/ledsPaired_morphine'
    
    if bsnm in np.concatenate([[aniIDDayBsnms[aniID][k] for k in aniIDDayBsnms[aniID]] for aniID in aniIDs2use_daGlu]):
        for aniID in aniIDs2use_daGlu:
            if bsnm in [aniIDDayBsnms[aniID][k] for k in aniIDDayBsnms[aniID]]:
                break
    else:
        aniID = getAniID(bsnm)    

    info = pd.read_csv(path, sep=' ', header=None)
    anis = np.asarray(info.loc[:, 0])
    leds = np.asarray(info.loc[:, 1])
    #
    for ind, ani in enumerate(anis):
        if ani == aniID:
            ledPairedOut = leds[ind]
    try:
        return ledPairedOut
    except:
        print('warning, jj LED paired assigned for ', Basename)
        return 0



def convert20KHzToEEG(x, toEEG=True):
    '''ccw edit on 19/10/18'''
    #return np.dot(x, (1250./20000.))
    x = np.array(x, dtype=int)
    if toEEG:
        return np.array(np.multiply(1250., np.divide(x, 20000)), dtype=int)
    else:
        return np.array(np.multiply(20000., np.divide(x, 1250.)), dtype=int)



def getBinEdgesOn(edge1, edge2, binMin, binMax, nBinsWithin):
    #
    edgeMin, edgeMax = [func([edge1, edge2]) for func in [np.min, np.max]]
    histEdges = np.linspace(edgeMin, edgeMax, nBinsWithin)
    histInt = float(np.diff(histEdges)[0])
    #
    edgeLims = []
    for edgeID, edge in zip(['min', 'max'], [binMin, binMax]):
        if edgeID == 'min':
            nEdges2add = np.ceil(np.abs(histEdges[0]-minVal) / histInt)+1
            edgeLims.append(histEdges[0]-nEdges2add*histInt)
        if edgeID == 'max':
            nEdges2add = np.ceil(np.abs(histEdges[-1]-maxVal) / histInt)+1
            edgeLims.append(histEdges[-1]+nEdges2add*histInt)
            #
    histEdges = np.arange(edgeLims[0], edgeLims[-1], histInt)
    return histEdges







def try_IO(call, tries=100):
    # for call() to retry:
    # ccw.try_IO(lambda: call())
    #
    assert tries > 0
    error = None
    result = None
    while tries:
        try:
            result = call()
        except IOError as e:
            error = e
            tries -= 1
        else:
            break
    if not tries:
        raise error
    return result


def getSpikesToTrigDiffs(spikes, triggers, maxLag):
    spikeDiffsOut = []
    for spikeTime in spikes:
        spikeDiffs = spikeTime-triggers
        minDiff = spikeDiffs[np.argmin(np.abs(spikeDiffs))]
        if abs(minDiff) < maxLag:
            spikeDiffsOut.append(minDiff)
    return spikeDiffsOut




def getEEGtracking(b):
    trk = vbf.load_tracking(b)
    #eegLen = eegHz*(np.size(trk,0)/39.)
    eegLen = np.size(trk,0)*32
    #trkTime = np.linspace(0, eegLen, np.size(trk,0))
    trk['t'] = np.linspace(0, eegLen, np.size(trk,0))
    trkSelInds = getNotNA(trk['x'], returnInds=True)
    trk = trk.loc[trkSelInds, :]
    trk = trk.reset_index(drop=True)
    trk['speed'] = np.concatenate([[0], getInstSpeed(trk)])
    return trk




def LoadIntervals(b,ext):
    ## Read in the file, store as <DataFrame> and return
    interval = pd.read_csv(b+ext, sep='\s+', header=None,names=['start','end'])
    interval['start'] = (interval['start']*(1250./20000)).astype(int)
    interval['end'] = (interval['end']*(1250./20000)).astype(int)
    return interval

def get_sucQuinLEDtimes(b):
    quiStr = '.LEDstrip_pulse'
    sucStr = '.LEDpanel_pulse'
    
    kTimes = {}
    for k, s in zip(['q', 's'], [quiStr, sucStr]):
        kTimes[k] = np.array(LoadIntervals(b, s))
    return kTimes

def get_sucQuinToneTimes(b):
    kTimes = {}
    for k, s in zip(['q', 's'], ['.toneQuinine_pulse', '.toneSucrose_pulse']):
        kTimes[k] = np.array(LoadIntervals(b, s))
    return kTimes




def getInstSpeed(trk, boxLenCm=46.):
    #
    pxPerCm = (np.max(trk['y']-np.min(trk['y'])))/boxLenCm
    tInt_s = (trk['t'].ravel()[-1] - trk['t'].ravel()[0])/(np.size(trk,0)*eegHz)
    def getDistCm(dx, dy):
        return np.sqrt(np.power(dx,2)+np.power(dx,2))/pxPerCm
    #
    def getSpeeds(coords):
        return [getDistCm(np.diff([coords[i, 0], coords[i+1, 0]])[0], np.diff([coords[i, 1], coords[i+1, 1]])[0])\
                /tInt_s for i in range(np.size(coords,0)-1)]
    #
    coords = np.column_stack([trk['x'], trk['y']])
    #
    return getSpeeds(coords)



def speedMap(trk, pixSq=30, smooth=[True, 0.5]):
    #
    trk['x'] = trk['x'] - np.min(trk['x'])
    trk['y'] = trk['y'] - np.min(trk['y'])
    #
    trk['speed'] = np.concatenate([[0], getInstSpeed(trk)])
    #
    trkRange = np.max([np.max(trk['x']), np.max(trk['y'])])*1.1
    pixBins = np.linspace(0, trkRange, pixSq+1)
    pixMap = np.zeros([pixSq]*2)
    pixCount = np.copy(pixMap)
    #
    xBins = [int(np.floor_divide(co, trkRange/float(pixSq))) for co in trk['x']] 
    yBins = [int(np.floor_divide(co, trkRange/float(pixSq))) for co in trk['y']] 
    trk['xPix'] = xBins
    trk['yPix'] = yBins
    #
    #
    for i, speed in enumerate(trk['speed']):
        pixMap[trk['xPix'][i], trk['yPix'][i]] += speed
        pixCount[trk['xPix'][i], trk['yPix'][i]] += 1
    #
    if smooth[0]:
        pixMap /= pixCount
        pixMap, _ = vbf.MatrixGaussianSmooth(rotate(pixMap, nRots=3), smooth[1])
    else:
        pixMap = rotate(pixMap/pixCount, nRots=3)
    #
    return pixMap


def nanCor(a, b, method='pearson'):
    return pd.DataFrame(np.column_stack([a, b])).corr(method=method).iloc[0,1]





'''
def getSpeedMod(strengths, trk, speedMax=20., speedInt=3., sampleWindow_ms=1000./512.):
    #
    nDims=np.size(strengths,0)
    sampleWindow_eeg = int(np.ceil(eegHz*sampleWindow_ms/1000.))
    #
    eegBins = getBins(0, np.size(strengths,1), sampleWindow_eeg)
    speedBins = getBins(0, speedMax, speedInt)
    #
    speedBinStrengths = [[] for i in range(np.size(speedBins,0))]
    for binSt, binEn in eegBins:
        trkSt, trkEn = [findNearest(t, trk['t'].ravel(), False) for t in [binSt, binEn]]
        speedb = np.nanmean(trk['speed'].ravel()[trkSt:trkEn])
        if speedb <= np.max(speedBins) and speedb != 0:
            speedBinInd = [i for i, b in enumerate(speedBins) if b[0] < speedb <= b[1]][0]
            speedBinStrengths[speedBinInd].append(
                [np.mean(strengths[i, binSt:binEn]) for i in range(nDims)])
    #
    speedBinStrengths = \
    [np.column_stack(speedBinStrengthsi) for speedBinStrengthsi in speedBinStrengths]
    speedLabs = [str(int(binSt))+'-'+str(int(binEn)) for binSt, binEn in speedBins]
    #
    return speedBinStrengths, speedLabs
'''

def getSpeedMod(strengths, trk, speedBins=None, speedMax=20., speedInt=3.):
    '''returns strengths[binInd][strengthInd, :]'''
    #
    if speedBins is None:
        speedBins = getBins(0, speedMax, speedInt)
    speedBinStrengths = []
    for binSt, binEn in speedBins:
        speedi = (binSt < trk['speed']) & (trk['speed'] <= binEn).ravel()
        eegInds = [int(round(ti)) for ti in trk['t'][speedi]]
        speedBinStrengths.append(np.column_stack([strengths[:, i] for i in eegInds]))
    speedLabs = [str(int(binSt))+'-'+str(int(binEn)) for binSt, binEn in speedBins]
    #
    return speedBinStrengths, speedLabs


def getSpeedMod2(strengths, speeds, speedBins=None, speedMax=20., speedInt=3.):
    '''returns strengths[binInd][strengthInd, :]
    speeds should be a 1D array, len() == np.size(strengths,1)
    '''
    #
    if speedBins is None:
        speedBins = getBins(0, speedMax, speedInt)
        #speedBins = getBins(0, speedMax, speedInt)
    speedBinStrengths = []
    for binSt, binEn in speedBins:
        speedInds = [i for i, speedi in enumerate(speeds) if (binSt < speedi) and (speedi <= binEn)]  #(binSt < speeds) & (speeds <= binEn)
        #eegInds = [int(round(ti)) for ti in trk['t'][speedi]]
        speedBinStrengths.append(np.column_stack([strengths[:, i] for i in speedInds]))
    speedLabs = [str(int(binSt))+'-'+str(int(binEn)) for binSt, binEn in speedBins]
    #
    return speedBinStrengths, speedLabs




def get_CrossCor(a, b, maxShift=250, shiftInt=1, samplingRate=1250., method='pearson'):
    '''shifts a rel 2 b

    RETURNS:
    displacements_ms, crossCor

    '''
    #
    traceLen = len(a)
    displacements = np.arange(-maxShift, maxShift+shiftInt, shiftInt)
    crossCor = []
    displacements_ms = []
    for d in displacements:
        if np.sign(d) == -1:
            a_ = np.concatenate([a, [np.nan for i in range(abs(d))]])
            b_ = np.concatenate([[np.nan for i in range(abs(d))], b])
        elif np.sign(d) == 1:
            a_ = np.concatenate([[np.nan for i in range(abs(d))], a])
            b_ = np.concatenate([b, [np.nan for i in range(abs(d))]])
        else:
            a_ = a
            b_ = b
        #
        crossCor.append(nanCor(a_, b_, method=method))
        displacements_ms.append(1000.*(d/samplingRate))
    #
    displacements_ms = np.array(displacements_ms)
    crossCor = np.array(crossCor)
    return displacements_ms, crossCor


def normVec(v):
    ed0 = eucDist([0]*len(v), v)
    return np.multiply(v, (1./ed0))


def nx_matrix_to_Graph(g):
    ''' weights of unconnected nodes and the diagonal should be set to -1 '''
    G = nx.Graph()
    nNodes = len(g)
    for i in range(nNodes):
        for j in range(i + 1, nNodes):
            if g[i, j] > 0:
                G.add_edge(i, j, weight=g[i, j])
    return G

def nx_plot_graph_(G, pos=None, labels=None):
    #pos = nx.spring_layout(G, seed=42)
    if pos is None:
        pos = nx.circular_layout(G)
    nx.draw(G, pos, with_labels=True, labels=labels, node_size=500, node_color='skyblue', font_size=10, font_weight='bold')
    nx.draw_networkx_edges(G, pos, width=1.0, alpha=0.5)
    edge_labels = nx.get_edge_attributes(G, "weight")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=14)
    
    

def getCrossDot(a, b, maxShift=250, shiftInt=1, samplingRate=1250.):
    #
    traceLen = len(a)
    displacements = np.arange(-maxShift, maxShift+shiftInt, shiftInt)
    crossDot = []
    displacements_ms = []
    for d in displacements:
        if np.sign(d) == -1:
            a_ = np.concatenate([a, [np.nan for i in range(abs(d))]])
            b_ = np.concatenate([[np.nan for i in range(abs(d))], b])
        elif np.sign(d) == 1:
            a_ = np.concatenate([[np.nan for i in range(abs(d))], a])
            b_ = np.concatenate([b, [np.nan for i in range(abs(d))]])
        else:
            a_ = a
            b_ = b
        #
        crossDot.append(np.nanmean(np.multiply(a_, b_)))
        displacements_ms.append(1000.*(d/samplingRate))
    #
    return displacements_ms, crossDot




def getInstAmpFreq2phase(instAmp, instFreq, phase, maxF=120, minF=14, nPhaseBins=32, logdf=0.5, normAmps=True):
    '''
    returns:
    phaseAmp, phaseCens, freqs
    '''
    smps = ~np.isnan(phase)
    phase0 = phase[smps]
    phasetg,freqs = myemd.getTG(nPh=32, minF=minF, maxF=maxF, logdf=logdf)
    if normAmps:
        phaseAmp, phaseCens, freqs = myemd.getTG(phase0, stats.zscore(instAmp[smps]), instFreq[smps], 
                                                 nPh=nPhaseBins, minF=minF, maxF=maxF, logdf=logdf)
    else:
        phaseAmp, phaseCens, freqs = myemd.getTG(phase0, instAmp[smps], instFreq[smps], 
                                                 nPh=nPhaseBins, minF=minF, maxF=maxF, logdf=logdf)
    return phaseAmp, phaseCens, freqs



'''
def binSpikes(spikes, binLen_ms=25, binFunc=np.sum):
    binLen_eeg = binLen_ms*eegHz/1000
    nDim = np.ndim(spikes) 
    if nDim == 1:
        binSt = np.arange(0, len(spikes), binLen_eeg, dtype=int)
        spikes = [spikes]
    else:
        binSt = np.arange(0, np.size(spikes,1), binLen_eeg, dtype=int)
    binEn = binSt[1:]
    bins = np.column_stack([binSt[:-1], binEn])
    binnedSpikes = []
    for spTrain in spikes:
        binnedTrain = []
        for binSt, binEn in bins:
            binnedTrain.append(binFunc(spTrain[binSt:binEn]))
        binnedSpikes.append(binnedTrain)
    if nDim == 1:
        return binnedSpikes[0]
    else:
        return np.row_stack(binnedSpikes)
'''
def get_psd(lfp, maxFreq=300, pointsPerHz=4, eegHz=eegHz):
    ''' returns freqAx, psd '''
    psdIndMax = maxFreq*pointsPerHz
    psd = sig.welch(lfp, fs=eegHz, nperseg=int(eegHz)*pointsPerHz)
    #
    return psd[0][0:psdIndMax], psd[1][0:psdIndMax]


def hilopass(x, f, sr, order=4):
    f = (2./sr)*f
    b,a = sig.butter(order, f, 'lowpass')
    x = sig.filtfilt(b, a, x)
    return x



def binSpikes(spikes, binSize_ms=25, binFunc=np.sum, spkHz=eegHz, returnCutLen=False):
    binLen_eeg = binSize_ms*spkHz/1000
    nDim = np.ndim(spikes) 
    if nDim == 1:
        binEdges = np.arange(0, len(spikes), binLen_eeg, dtype=int)
        spikes = [spikes]
    else:
        binEdges = np.arange(0, np.size(spikes,1), binLen_eeg, dtype=int)
    bins = np.column_stack([binEdges[:-1], binEdges[1:]])
    #
    binnedSpikes = []
    for spTrain in spikes:
        binnedTrain = []
        for binSt, binEn in bins:
            binnedTrain.append(binFunc(spTrain[binSt:binEn]))
        binnedSpikes.append(binnedTrain)
    if nDim == 1:
        binnedSpikes = np.array(binnedSpikes[0])
    else:
        binnedSpikes = np.row_stack(binnedSpikes)
    #
    if returnCutLen:
        cutLen = int(np.size(spikes,1) - float(bins[-1,1]))
        return binnedSpikes, cutLen
    else:
        return binnedSpikes






def cens2edges(cens):
    d=np.diff(cens)/2.
    edges = np.add(cens, d[0])
    edges = np.concatenate([[cens[0]-d[0]], edges])
    return edges



def edges2cens(edges):
    cens = np.convolve(edges,[.5,.5],'same')
    cens = cens[1::]
    return cens



def hist(obs, bins_, weights=None, norm=False):

    counts, edges = np.histogram(obs, bins_, weights=weights)
    cens = edges2cens(edges)
    if norm:
        counts = np.divide(counts, counts.sum(dtype=float))

    return cens, counts, edges

def get_edges_0split(edgeMax, nHistBins):
    edgeMax = np.abs(edgeMax)
    edges = np.concatenate([-np.linspace(0, edgeMax, nHistBins)[::-1][:-1], np.linspace(0, edgeMax, nHistBins)])
    return edges

def fix_data_for_de(data, y, f=np.concatenate):
    for k in data:
        data[k] = f(data[k])
    data[y] = data[y].astype(float)
    data = pd.DataFrame(data)
    return data

def sort_psth(psth, trigAx, win4cor=(-300, 300), zScore=True):
    if zScore:
        psth = stats.zscore(psth, axis=1)
    def get_cellRespScores():
        cellScores=[]
        for i, psthi in enumerate(psth):
            tSt, tEn = [findNearestInd(t, trigAx) for t in win4cor]
            cellScores.append(np.corrcoef(psthi[tSt:tEn], psth[:,tSt:tEn].mean(axis=0))[0,1])
        return np.array(cellScores)

    cellScores = get_cellRespScores()
    sortInds = np.argsort(np.argsort(-cellScores))
    psth_sort = np.zeros_like(psth)
    for i, sortInd in enumerate(sortInds):
        psth_sort[sortInd, :] = psth[i, :]
    return psth_sort




def get_seshIndsSleep(bsnm, rootFolder, sleep=False):
    '''use for not ccw desen...'''
    if sleep:
        seshInds = np.array([i for i, desen in enumerate(vbf.LoadStages(rootFolder+bsnm+'/'+bsnm)['desen']) \
                             if 's ' in desen])
    else:
        seshInds = np.array([i for i, desen in enumerate(vbf.LoadStages(rootFolder+bsnm+'/'+bsnm)['desen']) \
                             if 's ' not in desen])
    #
    return seshInds



def getPostSeshLabInd(stageID, stages):
    #
    if stageID in ['post', 'rec']:
        si = [i for i, desen in enumerate (stages['desen'].ravel()) if 'condtestledloc' in desen.lower()][0]
    if stageID in ['iext', 'ext']:
        si = [i for i, desen in enumerate (stages['desen'].ravel()) if 'exttestledloc' in desen.lower()][0]
    if stageID == 'ren':
        inds = np.where(['phase_laser2' in desen for desen in stages['desen'].ravel()])[0]
        if len(inds):
            si = inds[0]
        else:
            si = [i for i, desen in enumerate (stages['desen'].ravel()) if 'condtestledloc' in desen.lower()][-1]
    elif stageID == 'bl':
        try:
            si = [i for i, desen in enumerate (stages['desen'].ravel()) if 'fam' in desen.lower()][0]
        except:
            print('warning: no sesh found for: '+stageID)
            si = None
    #
    return si



def getCondSeshLabInd(stageID, stages):
    #
    try:
        si = [i for i, desen in enumerate (stages['desen'].ravel()) if stageID.lower() in desen.lower()][0]
    except:
        print('warning: no sesh found for: '+stageID)
        si = None
    return si




def getPreSeshLabInd(stageID, stages):
    #
    if stageID in ['il', 'pre']:
        try:
            si = [i for i, desen in enumerate (stages['desen'].ravel()) if 'condtestledloc' in desen.lower()][0]
        except:
            print('warning: no sesh found for: '+stageID)
            si = None
    #
    return si




def getSeshLab4testStage(bsnm, testStage, rootFolder, mntStr=None, mnfs=False):
    '''for merged data'''

    tetRegs, tetLabs, seshLabs, tetUseInds = get_bsnmInfo(bsnm, rootFolder, multi_site=False)

    stages = vbf.LoadStages(rootFolder+bsnm+'/'+bsnm)
    if testStage in ['rec', 'iext', 'ext', 'ren']:
        seshLab = seshLabs[getPostSeshLabInd(testStage, stages)]
    elif testStage in ['il', 'pre']:
        seshLab = seshLabs[getPreSeshLabInd(testStage, stages)]
    elif testStage in ['sal', 'coc']:
        seshLab = seshLabs[getCondSeshLabInd(testStage, stages)]
    else:
        print('test stage not recognised')
        seshLab = None
    return seshLab

def get_seshLab4testStage_cpp(bsnm, testStage):
    rootFolder = get_rootFolder4bsnm(bsnm)
    desen = vbf.LoadStages(rootFolder+bsnm+'/'+bsnm)['desen'].ravel()
    tetRegs, tetLabs, seshLabs, tetUseInds = get_bsnmInfo(bsnm, rootFolder)
    if testStage == 'pre':
        seshLab = seshLabs[[i for i, d in enumerate(desen) if 'test' in d.lower()][0]]
    else:
        if testStage == 'rec':
            seshi = [i for i, d in enumerate(desen) if 'condtest' in d.lower()][0]
        elif testStage == 'ext':
            seshi = [i for i, d in enumerate(desen) if 'exttest' in d.lower()][0]
        else:
            seshi = [i for i, d in enumerate(desen) if 'condtest' in d.lower()][-1]
        seshLab = seshLabs[seshi]
    return seshLab
    
    

def get_inStage(testStage):
    if testStage == 'iext':
        def in_stage(en, seshLen):
            return en < 20*60*eegHz
    else:
        def in_stage(en, seshLen):
            return en>(seshLen-20*60*eegHz)
    return in_stage




def getSeshInd(bsnm, rootFolder):
    '''to get a ~baseline sesh for any animal (merged data)'''
    stages = vbf.LoadStages(rootFolder+bsnm+'/'+bsnm)
    if bsnm in ['mhb17-171222', 'mrr02-180503', 'mrr04-180520', 'mccw03-180219', \
                'mccw05-180525', 'mrr03-180602', 'mrr05-181125']:
        seshIndOut = getPostSeshLabInd('rec', stages)
    else:
        seshIndOut = [i for i, desen in enumerate(stages['desen'].ravel()) \
                      if desen[:3] != 's b'][0]
    #
    return seshIndOut



def getInstAmp(trace):
    analyticSignal = sig.hilbert(trace, axis=0)
    instAmp = np.abs(analyticSignal)
    return instAmp


def getInstFreq(trace, smooth=True, eegHz=1250.):
    '''Based on Vitor defaults '''
    analyticSignal = sig.hilbert(trace, axis=0)
    iphase = np.unwrap(np.angle(analyticSignal), axis=0)
    orig_dim = iphase.ndim
    if iphase.ndim == 2:
        iphase = iphase[:, :, None]
    #
    if smooth:
        sig.medfilt(iphase, 5)
    if orig_dim == 2:
        iphase = iphase[:, :, 0]
    #
    iphase = iphase + np.pi/2
    #
    iphase = np.gradient(iphase, axis=0)
    instFreq = iphase / (2.0*np.pi) * eegHz
    return instFreq








def detectCellAssemblies(cellSpikesb, return_all=False):
    '''Vitor method
    if return_all:
        return eigVals, eigVecs, lambdaMax, assemblies
    return assemblies
    '''
    #
    cellSpikesb = stats.zscore(cellSpikesb, axis=1)
    corMat = np.corrcoef(cellSpikesb)
    eigVals, eigVecs = np.linalg.eigh(corMat)
    eigVals = eigVals[::-1]
    eigVecs = eigVecs[::-1]
    #
    lambdaMax = np.power((1 + (np.sqrt(np.size(corMat,0)/float(np.size(cellSpikesb,1))))), 2)
    #
    #eigVecsUse = np.array([ eigVecs[ind] for ind, eVal in enumerate(eigVals) if eVal > lambdaMax ]).T
    #projMat = np.dot(corMat, eigVecsUse)
    #model = FastICA(n_components=np.size(eigVecsUse, 0))
    #assemblies = model.fit(projMat).transform(projMat).T
    nICs = np.sum(eigVals > lambdaMax)
    if nICs == 0:
        print('WARNING: no significant patterns detected')
        assemblies = None
    else:
        model = FastICA(n_components=nICs)
        assemblies = model.fit(cellSpikesb).transform(cellSpikesb).T
        assemblies = np.row_stack([flipSign(ass) for ass in assemblies])
    #
    if return_all:
        return eigVals, eigVecs, lambdaMax, assemblies
    return assemblies

def get_assemblyStrength(ass, cellSpikesb):
    op = np.outer(ass, ass)
    np.fill_diagonal(op, 0)
    strength = np.array([np.dot(np.dot(op, act_t), act_t.T) for act_t in cellSpikesb.T])
    return strength


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

    #x = np.concatenate(x)
    #y = np.concatenate(y)

    return x, y

def plot_histFill(cens, Ys, lw=2, color='k', fill=[True, 'white', 0.5]):
    x, y = get_xy_4_histFill(cens, Ys)
    plt.plot(x, y, color=color, zorder=3, lw=lw)
    if fill[0]:
        plt.fill_between(x, np.zeros_like(x), y, color=fill[1], alpha=fill[2], zorder=3)





def optimiseSimMat(patternsX, patternsY, method='euclidean'):
    '''
    returns:
    simMat, rowInds, colInds
    '''

    if method == 'cosine':
        from sklearn.metrics.pairwise import cosine_similarity as getsim
        inputs = {'X': patternsX, 'Y': patternsY}
        simMat0 = getsim(**inputs)
    elif method == 'pearson':
        from sklearn.metrics import pairwise_distances
        simMat0 = pairwise_distances(patternsX, patternsY, metric='correlation')
    elif method == 'euclidean':
        simMat0 = np.zeros([patternsX.shape[0], patternsY.shape[0]])
        for x, vecx in enumerate(patternsX):
            for y, vecy in enumerate(patternsY):
                simMat0[x, y] = eucCloseness(vecx, vecy)
    else:
        print(method +' for similarity has not been implemented yet.')
        return
    #
    #
    #
    def fillmissingidxs(ind,n):
        missing = list(set(np.arange(n))-set(ind))
        ind = np.array(list(ind)+missing)
        return ind

    import scipy.optimize as optimize
    rowind,colind = optimize.linear_sum_assignment(simMat0)

    rowInds = fillmissingidxs(rowind,np.size(simMat0,0))
    colInds = fillmissingidxs(colind,np.size(simMat0,1))

    simMat = np.zeros_like(simMat0)

    for rowi, coli in zip(rowInds, colInds):
        simMat[:, coli] = simMat0[rowi]
    #
    return simMat, rowInds, colInds







####################
### IMF PLOTTING ###
####################

def plotIMFs(emd, IMFcols=None, Title=None, traceLen_s=1., eegHz=1250., figSize=2):
    nImfs, recLen = np.shape(emd)
    if IMFcols is None:
        IMFcols = sb.color_palette('husl', nImfs)
    #
    tracesLen = int(traceLen_s*eegHz)
    traceStartInd = rdm.randint(0, (recLen-tracesLen))
    plt.figure(figsize=(10*figSize, nImfs*figSize))
    for imfInd, imfTrace in enumerate(emd):
        seshi = seshLabs[0]
        imfSample = emd[imfInd, traceStartInd:(traceStartInd+tracesLen)]
        plt.subplot(nImfs, 1, imfInd+1)
        plt.plot(imfSample, color=IMFcols[imfInd])
        #if imfInd == 0:
        #    plt.suptitle = Title
        plt.axis('off')


def plotImfPsds(emd, IMFcols=None, Title=None, pointsPerHz=2, maxFreq=100):
    nImfs, recLen = np.shape(emd)
    if IMFcols is None:
        IMFcols = sb.color_palette('husl', nImfs)
    psdIndMax = maxFreq*pointsPerHz
    def zScore(vec):
        SD = np.std(vec)
        M = np.mean(vec)
        zOut = []
        for i in vec:
            zOut.append((i-M)/SD)
        return zOut
    for imfInd, imfTrace in enumerate(emd):
        psd = sig.welch(imfTrace, fs=eegHz, nperseg=int(eegHz)*pointsPerHz)
        fMax = str(int(psd[0][np.argmax(psd[1])]))+' Hz'
        plt.plot(psd[0][0:psdIndMax], zScore(psd[1])[0:psdIndMax], color=IMFcols[imfInd], 
                 label=fMax)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Power (z-score)')
        #fig, ax = plt.plot(psd[0][0:psdIndMax], psd[1][0:psdIndMax], color=imfCols[imfInd], label=imfFreqLabs[imfInd])
        #ax.set_facecolor('red')
    plt.legend(title='IMF frequencies', frameon=True, facecolor='white')

def plot_title(title, loc='left', fontweight='bold', color='gray', fontsize=12):
    plt.title(title, loc=loc, fontweight=fontweight, color=color, fontsize=fontsize)

def plot_xlim(x, f=0.5):
    d = (x[1]-x[0])*f
    plt.xlim(x[0]-d, x[-1]+d)
    
def plot_vec(vec, cols=None, yshift=0., xshift=0., alpha=1, vertical=True, lw=None, s=None, zorder=3): #, ms=None):
    #
    if cols is None:
        cols = np.repeat('gray', len(vec))
    if vertical:
        for i, col in enumerate(cols):
            plt.hlines(y=i+yshift, xmin=0.+xshift, xmax=vec[i]+xshift, color=col, alpha=alpha, lw=lw, zorder=zorder)
            plt.scatter(vec[i]+xshift, i+yshift, color=col, alpha=alpha, s=s, zorder=zorder) 
            #plt.plot(vec[i]+xshift, i+yshift, 'o', color=col, alpha=alpha, ms=ms, edgecolors=col)
        plt.vlines(x=xshift, ymin=-0.5+yshift, ymax=len(cols)-0.5+yshift, color='gray', lw=lw, alpha=alpha, zorder=zorder-1)
    else:
        for i, col in enumerate(cols):
            plt.vlines(x=i+xshift, ymin=0.+yshift, ymax=vec[i]+yshift, color=col, alpha=alpha, lw=lw, zorder=zorder)
            plt.scatter(i+xshift, vec[i]+yshift, color=col, alpha=alpha, s=s, zorder=zorder) #plt.plot(i+xshift, vec[i]+yshift, 'o', color=col, alpha=alpha, ms=ms)

        plt.hlines(y=yshift, xmin=-0.5+xshift, xmax=len(cols)-0.5+xshift, color='gray', lw=lw, alpha=alpha, zorder=zorder-1)



def get_minMaxFreqs4imfTit(imfTit, region=None):
    if imfTit in ['slow-gamma']:
        if region is None:
            minFreq, maxFreq = None, None
        else:
            if region == 'vta':
                minFreq, maxFreq = (20, 35)
            else:
                minFreq, maxFreq = (25, 40)          
    elif imfTit == 'beta':
        minFreq, maxFreq = (12, 25)
    elif imfTit == 'theta':
        minFreq, maxFreq = (4, 12)
    elif imfTit == '4Hz':
        minFreq, maxFreq = (3.5, 5.5)
    else:
        minFreq, maxFreq = None, None
    return minFreq, maxFreq

def get_minMaxFreqs_propCyHz(x, sr, percThresh=.66, nBins=100):
    cycles0 = getCycles(x, 0, 100, 't')[0]
    cyHzs = 1./(np.subtract(cycles0[:,1], cycles0[:,0])/sr)
    counts, cens, _ = vbf.hist(cyHzs, np.linspace(np.percentile(cyHzs, 1), np.percentile(cyHzs, 99), nBins))
    thresh = counts.sum()*percThresh
    inds = np.array([counts.argmax()])
    for jj in range(len(counts)):
        toL, toR = [inds[i] for i in [0, -1]]

        if counts[toL] > counts[toR] and toL != 0:
            inds = np.concatenate([[toL-1], inds])

        else:
            if toR != len(counts)-1:
                inds = np.concatenate([inds, [toR+1]])

        if counts[inds].sum() > thresh:
            break

    minFreq, maxFreq = [inds[i] for i in [0, -1]]
    return minFreq, maxFreq



def fix_sbHeatmap(pad=0):
    b, t = plt.ylim()
    b += (0.5+pad)
    t -= (0.5+pad)
    plt.ylim(b, t)
    #
    l, r = plt.xlim()
    r += (pad)
    l -= (pad)
    plt.xlim(l, r)


def runICA(dims, nICs, max_iter=500):
    from sklearn.decomposition import FastICA
    model = FastICA(n_components=nICs, max_iter=max_iter)
    proj = model.fit(dims.T).transform(dims.T)
    return proj.T

def get_boutWindows(boutStatus):
    boutStartInds = list(np.array(locateMotif([False, True], boutStatus))+1)
    boutEndInds = list(np.array(locateMotif([True, False], boutStatus)))
    if all([len(boutStartInds), len(boutEndInds)]): #techincall might want option to incliude start/end
        if boutEndInds[0] < boutStartInds[0]:
            boutEndInds = boutEndInds[1:]
        if boutStartInds[-1] > boutEndInds[-1]:
            boutEndInds.append(len(boutStatus)-1)
        bouts = np.column_stack([boutStartInds, boutEndInds])
    else:
        bouts = []
    #
    return bouts

def getBurstInfo(amp, percThresh=75): 
    thresh = np.percentile(amp, percThresh)
    aboveThresh = amp>=thresh
    bursts = getBoutWindows(aboveThresh)
    burstDurs = np.array([1000.*((ben-bst)/eegHz) for bst, ben in bursts])
    #
    return bursts, burstDurs






############
## GRAPHS ##
############

def getTraceAv(traceList):
    if len(traceList) != 0:
        sampleLen = len(traceList[0])
        traceAv = np.zeros(sampleLen)
        for trace in traceList:
            traceAv = traceAv+np.asarray(trace)
        traceAv = traceAv/len(traceList)
    if len(traceList) == 0:
        print('No traces in input') # traceAv = np.zeros(traceLength)
    return traceAv




def getSpecAv(traceList, freqs, eegHz=1250.):
    # this bad as loose on smaller freqs - better to get spec for whole rec then use inds
    # to cut and average spec
    if len(traceList) != 0:
        avSpecs = np.zeros((len(freqs), len(traceList[0])))
        for trace in traceList:
            spec = vbf.WvSpectrogram(trace, eegHz, freqs)[0]
            avSpecs = avSpecs+spec
        avSpecs = avSpecs/len(traceList)
    if len(traceList) == 0:
        print('No traces in input') # avSpecs = np.zeros((len(freqs), traceLen))
    #
    return avSpecs




def trigTraceAndSpecAv(Signal2Trig, Trigger, TrigAvLength, SpecFreqs=np.arange(10, 125, 5), SamplingRate=1250., Z_score=False, returnSpec=True):
    '''Returns triggered trace av and the triggered spectrogram Average, along with freq/time axis. Currently a bug with z scoring. 
Currently Signal2Trig must be 1D'''
    #
    if Z_score is True:
        Signal2Trig = zScore(Signal2Trig)
    #
    TimeSamples = (np.arange(TrigAvLength)-(TrigAvLength/2.)).astype(int)
    Trigger = Trigger[Trigger>(1+np.ceil(TrigAvLength/2))]
    Trigger = Trigger[Trigger<(len(Signal2Trig)-1-np.ceil(TrigAvLength/2))] # chnage here if ndim>1
    #
    trigTraceAv = np.zeros((len(TimeSamples)))
    trigSpecAv = np.zeros((len(SpecFreqs), len(trigTraceAv)),dtype=complex)
    #
    nTriggers = len(Trigger)
    #
    for triggerIndi, triggeri in enumerate(Trigger):
        samplesi = Trigger[triggerIndi]+TimeSamples
        #
        trigTracei = Signal2Trig[samplesi]
        trigTraceAv = trigTraceAv + trigTracei
        #
        if returnSpec:
            trigSpeci = vbf.WvSpectrogram(trigTracei, SamplingRate, SpecFreqs)[0]
            trigSpecAv = trigSpecAv + trigSpeci
        #
    trigTraceAv = trigTraceAv/nTriggers
    if returnSpec:
        trigSpecAv = trigSpecAv/nTriggers
    timeAxis = TimeSamples*1000/SamplingRate
    #
    if returnSpec:
        return trigTraceAv, trigSpecAv, timeAxis, SpecFreqs
    else:
        return trigTraceAv, timeAxis



def plotSEpoly(Xs, Ys, SEs, col):
        above = np.array([y+se for y, se in zip(Ys, SEs)])
        below = np.array([m-se for m, se in zip(Ys, SEs)])
        poly_x = np.concatenate([Xs, np.flip(Xs, 0)])
        poly_y = np.concatenate([above, np.flip(below, 0)])
        plt.fill(poly_x, poly_y, color=col)





def get_imfMainFreqs(imfs):
    mainFreqs = []
    for imf in imfs.T:
        psd = sig.welch(imf[:125000], fs=eegHz, nperseg=int(eegHz*1))
        mainFreqs.append(psd[0][np.argmax(psd[1])])
    return np.array(mainFreqs)

'''
def getIMFstrengths(imfTraces, fMax=120, fMin=2, fStep=2):
    #
    specFreqs=np.arange(fMin, fMax+fStep, fStep)
    imfStrengthsOut = np.zeros(np.shape(imfTraces))
    for imfi, imfTracei in enumerate(imfTraces):
        speci = vbf.WvSpectrogram(imfTracei, eegHz, specFreqs)[0]
        psdi = sig.welch(imfTracei[:125000], fs=eegHz, nperseg=int(eegHz*(1./fStep)))
        maxFi = [i for i, f in enumerate(psdi[0]) if f == fMax][0]
        minFi = [i for i, f in enumerate(psdi[0]) if f == fMin][0]
        psdVec = psdi[1][minFi:(maxFi+1)]
        imfStrengthsOut[imfi, :] = \
        stats.zscore([np.dot(psdVec, specit) for specit in speci.T])
    return imfStrengthsOut
'''



def getImfPowStrengths(imfTraces, fMax=120, fMin=2, fStep=2, usePSD=False):
    #
    specFreqs=np.arange(fMin, fMax+fStep, fStep)
    imfStrengthsOut = np.zeros(np.shape(imfTraces))
    for imfi, imfTracei in enumerate(imfTraces):
        speci = vbf.WvSpectrogram(imfTracei, eegHz, specFreqs)[0]
        if usePSD:
            psdi = sig.welch(imfTracei[:125000], fs=eegHz, nperseg=int(eegHz*(1./fStep)))
            maxFi = [i for i, f in enumerate(psdi[0]) if f == fMax][0]
            minFi = [i for i, f in enumerate(psdi[0]) if f == fMin][0]
            psdVec = psdi[1][minFi:(maxFi+1)]
            imfStrengthsOut[imfi, :] = \
            stats.zscore([np.dot(psdVec, specit) for specit in speci.T])
        else:
            imfStrengthsOut[imfi, :] = stats.zscore([np.sum(spect) for spect in speci.T])
    return imfStrengthsOut



## think ok, maybe bug...
def getBinEdgesOn(edge1, edge2, binMin, binMax, nBinsWithin):
    #
    edgeMin, edgeMax = [func([edge1, edge2]) for func in [np.min, np.max]]
    histEdges = np.linspace(edgeMin, edgeMax, nBinsWithin)
    histInt = float(np.diff(histEdges)[0])
    #
    edgeLims = []
    for edgeID, edge in zip(['min', 'max'], [binMin, binMax]):
        if edgeID == 'min':
            nEdges2add = np.ceil(np.divide(np.abs(histEdges[0]-binMin), histInt))+1
            edgeLims.append(histEdges[0]-nEdges2add*histInt)
        if edgeID == 'max':
            nEdges2add = np.ceil(np.divide(np.abs(histEdges[-1]-binMax), histInt))+1
            edgeLims.append(histEdges[-1]+nEdges2add*histInt)
            #
    histEdges = np.arange(edgeLims[0], edgeLims[-1], histInt)
    return histEdges




def getBins(binMin, binMax, binInt):
    binEdges = np.arange(binMin, binMax+binInt, binInt)
    binsOut = np.column_stack([binEdges[:-1], binEdges[1:]])
    return binsOut



def getTracesSpecs(bsnm, seshLab, rootFolder, specFreqs=np.arange(5, 90, 10), returnSpec=False, eegHz=1250.,
                   divMean=True, subBL=False, concatenateSpec=True):
    #
    par = vbf.LoadPar(rootFolder+bsnm+'/'+bsnm)
    base = rootFolder+bsnm+'/'+bsnm+'_'+seshLab
    lfpInds, _ = vbf.GetLFPis(base, ChAreas='all')
    lfp = vbf.MapLFPs(base+'.eeg', par['nch'])
    tetUseInds = getAniTetUseInds(bsnm)
    #
    tracesOut = []
    for tetInd in tetUseInds:
        tracesOut.append(lfp[lfpInds[tetInd], :])
    tracesOut = np.row_stack(tracesOut)
    #
    if returnSpec:
        specOut = []
        for trace in tracesOut:
            regSpec = vbf.WvSpectrogram(trace, eegHz, specFreqs)[0]
            if divMean:
                regSpec = regSpec/np.mean(regSpec)
            if subBL:
                regSpec = specSubBL(regSpec)
            specOut.append(regSpec)
        #
        if concatenateSpec:
            specOut = np.row_stack(specOut)
        return tracesOut, specOut
    else:
        return tracesOut



def get_1stInBurst(spikeTimes, minDiff_ms=220, sr=eegHz):
    isis_ms = 1000*(np.diff(spikeTimes)/sr)
    return spikeTimes[1:][np.where(isis_ms>minDiff_ms)[0]]


def specSubBL(spec):
    #
    bl = [np.mean(powf) for powf in spec]
    #specOut = np.column_stack([np.subtract(spect, bl) for spect in specj.T])
    #return np.column_stack([np.subtract(spect, bl) for spect in spec.T]) #specOut
    #
    return np.row_stack([np.subtract(specf, bl[fi]) for fi, specf in enumerate(spec)])




def triggeredAverage(Signal2Trig, Trigger, TrigAvLength, trigSum=False, SamplingRate=1250.):
    #
    shape = np.shape(Signal2Trig)
    TimeSamples = (np.arange(TrigAvLength)-(TrigAvLength/2.)).astype(int)
    #
    Trigger = Trigger[Trigger>(1+np.ceil(TrigAvLength/2))]
    if np.ndim(Signal2Trig)>1:
        Trigger = Trigger[Trigger<(shape[1]-1-np.ceil(TrigAvLength/2))]
        TriggAverage = np.zeros((shape[0],len(TimeSamples)))
    else:
        Trigger = Trigger[Trigger<(len(Signal2Trig)-1-np.ceil(TrigAvLength/2))]
        TriggAverage = np.zeros((len(TimeSamples)))

    nTriggers = len(Trigger)
    #
    for triggeri in range(nTriggers):
        samples = Trigger[triggeri]+TimeSamples

        if np.ndim(Signal2Trig)>1:
            TriggAverage = TriggAverage+Signal2Trig[:,samples]
        else:
            TriggAverage = TriggAverage+Signal2Trig[samples]
    if not trigSum:
        TriggAverage = TriggAverage/nTriggers
    #
    TimeAxis = TimeSamples*1000/SamplingRate
    #
    return TimeAxis,TriggAverage

def triggeredMatrix(Signal2Trig,Trigger,TrigAvLength,SamplingRate=1250.):

    shape = np.shape(Signal2Trig)

    TimeSamples = (np.arange(TrigAvLength)-(TrigAvLength/2.)).astype(int)

    Trigger = Trigger[Trigger>(1+np.ceil(TrigAvLength/2))]
    if np.ndim(Signal2Trig)>1:
        Trigger = Trigger[Trigger<(shape[1]-1-np.ceil(TrigAvLength/2))]
        TriggMat = np.zeros((len(Trigger),shape[0],len(TimeSamples)))
    else:
        Trigger = Trigger[Trigger<(len(Signal2Trig)-1-np.ceil(TrigAvLength/2))]
        TriggMat = np.zeros((len(Trigger),len(TimeSamples)))

    nTriggers = len(Trigger)

    for triggeri in range(nTriggers):
        samples = Trigger[triggeri]+TimeSamples

        if np.ndim(Signal2Trig)>1:
            TriggMat[triggeri] = Signal2Trig[:,samples]
        else:
            TriggMat[triggeri] = Signal2Trig[samples]

    TimeAxis = TimeSamples*1000/SamplingRate

    return TimeAxis, TriggMat


def get_trigAx(trigLen, sr):
    trigAx = triggeredAverage(np.zeros(trigLen*10), np.array([trigLen+3]), trigLen, SamplingRate=sr)[0]
    return trigAx


def getSpikePhaseCoherence(spikePhases, minPhases=[False, 200, np.nan]):
    if minPhases[0] and len(spikePhases) < minPhases[1]:
        return minPhases[2]
    return np.abs(np.mean(np.exp(1j*spikePhases)))



# default for plotBarcode()
gridProps_default={'col' : 'k', 
                   'lw' : 2, 
                   'colBord' : 'k',
                   'lwBord' : 2}



def plot_barcodeTit(sti, stateCols, fontsize=12, loc='left', fontweight='bold'):
    plt.title('Barcode '+str(sti+1), color=stateCols[sti], fontsize=fontsize, loc=loc, fontweight=fontweight)



def switch_barcode_regOrder(state, regImfUseInds, intDimInds, regions_from=regions_old, regions_to=regions[::-1]):
    #
    import ccw_imfBarcode_functions as barcode
    intVec = barcode.barcode2vec(state, intDimInds)

    dimInfo_from = []
    for region in regions_from:
        for imfi in regImfUseInds[region]:
            dimInfo_from.append({'region' : region, 
                                 'imfi' : imfi})


    def get_toi_4_fromi(fromi):
        region0, imfi0 = [dimInfo_from[fromi][k] for k in ['region', 'imfi']]
        toi = 0
        for region in regions_to:
            if region != region0:
                toi += len(regImfUseInds[region])
            else:
                toi += np.where(regImfUseInds[region0] == imfi0)[0][0]
                break
        return toi

    intDimInds_ = np.row_stack([[get_toi_4_fromi(fromi) for fromi in [xi, yi]] for xi, yi in intDimInds])

    state_ = barcode.vec2barcode(intVec, intDimInds_)

    return state_, regions_to



def plot_barcode(state, regImfUseInds, intDimInds=None, switchRegOrder=True, regFontsize=10, regGrid=True, cmap='RdBu_r', 
                setvRange=[False, None, None], cbar=False, addLabs=True, xDimLabsRot=None, yDimLabsRot=None, mask_tri=False,
                useDimLabs=False, dimLabsRight=True, dimLabsTop=True, dimFontsize=10, regions_=regions_old, gridProps=None):

    if switchRegOrder:
        state, regions_ = switch_barcode_regOrder(state, regImfUseInds, intDimInds)
    #    
    if mask_tri:
        state = np.tril(state, k=-1)
    if gridProps is None:
        gridProps = gridProps_default
    else:
        gridProps_ = gridProps_default
        for k in gridProps:
            gridProps_[k] = gridProps[k]
        #

        del gridProps
        gridProps = gridProps_
    #fig, ax = plt.subplots()
    reg_nDims=np.array([len(regImfUseInds[region]) for region in regions_])

    nRegs = len(regImfUseInds)
    #
    regTits = []
    regCols = []
    dimLabs = []
    dimLabCols = []

    for region in regions_:
        if region in list(regImfUseInds.keys()):
            regTit, regCol = getTitColForReg(region)
            regTits.append(regTit)
            regCols.append(regCol)
            for imfi in regImfUseInds[region]:
                dimLabs.append(regTit+'_'+get_imfTit4imfi(imfi, region))
                dimLabCols.append(regCol)
    #
    #
    if setvRange[0]:
        vmin = setvRange[1]
        vmax = setvRange[2]
        #
    else:
        maxWeight = np.max(np.abs(state))
        vmin = -maxWeight
        vmax = maxWeight
    #
    sb.heatmap(state, cmap=cmap, square=True, vmin=vmin, vmax=vmax, cbar=cbar,
               xticklabels=False, yticklabels=False)
    #
    if addLabs:
        if useDimLabs:
            #
            labMids = np.array(range(len(dimLabs)))+0.5

            if dimLabsRight:
                if dimLabsTop:
                    plt.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False, 
                                    left=False, right=True, labelleft=False, labelright=True)
                else:
                    plt.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True, 
                                    left=False, right=True, labelleft=False, labelright=True)

            else:
                if dimLabsTop:
                    plt.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False, 
                                    left=True, right=False, labelleft=True, labelright=False)
                else:
                    plt.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True, 
                                    left=True, right=False, labelleft=True, labelright=False)

            if dimLabsTop:
                if xDimLabsRot is None:
                    xDimLabsRot=-90
            else:
                if xDimLabsRot is None:
                    xDimLabsRot = 90


            plt.xticks(labMids, dimLabs, rotation=xDimLabsRot, fontsize=dimFontsize)
            if yDimLabsRot is None:
                yDimLabsRot = 0
            plt.yticks(labMids, dimLabs, rotation=yDimLabsRot, fontsize=dimFontsize)
            for tickLab, col in zip(plt.gca().get_xticklabels(), dimLabCols):
                tickLab.set_color(col)
            for tickLab, col in zip(plt.gca().get_yticklabels(), dimLabCols):
                tickLab.set_color(col)
            #
            #
        else:
            midShifts = np.array([dimsPerReg/2. for dimsPerReg in reg_nDims])
            regMids = np.subtract(np.cumsum(reg_nDims), midShifts)
            #regMids = np.arange(midShift, dimsPerReg*nRegs, dimsPerReg)
            plt.xticks(regMids, regTits, fontsize=regFontsize, fontweight='bold')
            for tickLab, regCol in zip(plt.gca().get_xticklabels(), regCols):
                tickLab.set_color(regCol)
            #
            plt.yticks(regMids, regTits, fontsize=regFontsize, fontweight='bold')
            for tickLab, regCol in zip(plt.gca().get_yticklabels(), regCols):
                tickLab.set_color(regCol)
    #
    indMin = 0
    indMax = reg_nDims.sum()

    plt.xlim(indMin, indMax)
    plt.ylim(indMin, indMax)
    if regGrid:
        # plot inside grid
        regSeps = np.concatenate([[0], np.cumsum(reg_nDims)])
        for sep in regSeps:
            plt.hlines(y=sep, xmin=0, xmax=np.size(state,0), color=gridProps['col'], lw=gridProps['lw'])
            plt.vlines(x=sep, ymin=0, ymax=np.size(state,1), color=gridProps['col'], lw=gridProps['lw'])

        # plot Border
        '''plt.axhline(y=indMin, color=gridProps['colBord'],linewidth=gridProps['lwBord'])
        plt.axhline(y=state.shape[1], color=gridProps['colBord'],linewidth=gridProps['lwBord'])
        plt.axvline(x=indMin, color=gridProps['colBord'],linewidth=gridProps['lwBord'])
        plt.axvline(x=state.shape[0], color=gridProps['colBord'],lw=gridProps['lwBord'])'''
        #
        plt.hlines(y=[indMin, indMax], xmin=indMin, xmax=indMax, color=gridProps['colBord'], lw=gridProps['lwBord'])
        plt.vlines(x=[indMin, indMax], ymin=indMin, ymax=indMax, color=gridProps['colBord'], lw=gridProps['lwBord'])


def plotSqareBarcode(state, imfUseInds=[0], regFontsize=10, regGrid=True, cmap='RdBu_r', setvRange=[False, None, None], 
                cbar=False, addLabs=True, useDimLabs=False, regions=regions, gridProps=None):
    #
    #
    if gridProps is None:
        gridProps = gridProps_default
    else:
        gridProps_ = gridProps_default
        for k in gridProps:
            gridProps_[k] = gridProps[k]
        #
        del gridProps
        gridProps = gridProps_
    #fig, ax = plt.subplots()
    dimsPerReg=len(imfUseInds)
    nRegs = len(regions)
    #
    regTits = []
    regCols = []
    dimLabs = []
    dimLabCols = []
    for region in regions:
        regTit, regCol = getTitColForReg(region)
        regTits.append(regTit)
        regCols.append(regCol)
        for imfi in imfUseInds:
            #dimLabs.append(regTit+'_IMF'+str(imfi+1))
            dimLabs.append(regTit+'_'+getTit4imfi(imfi))
            dimLabCols.append(regCol)
    #
    #
    if setvRange[0]:
        vmin = setvRange[1]
        vmax = setvRange[2]
        #
    else:
        maxWeight = np.max(np.abs(state))
        vmin = -maxWeight
        vmax = maxWeight
    #
    sb.heatmap(state, cmap=cmap, square=True, vmin=vmin, vmax=vmax, cbar=cbar,
               xticklabels=False, yticklabels=False)
    #
    if addLabs:
        if useDimLabs:
            #
            labMids = np.array(range(len(dimLabs)))+0.5
            plt.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False, 
                            left=False, right=True, labelleft=False, labelright=True)
            plt.xticks(labMids, dimLabs, rotation=-90)
            plt.yticks(labMids, dimLabs, rotation=0)
            for tickLab, col in zip(plt.gca().get_xticklabels(), dimLabCols):
                tickLab.set_color(col)
            for tickLab, col in zip(plt.gca().get_yticklabels(), dimLabCols):
                tickLab.set_color(col)
            #
            #
        else:
            midShift = dimsPerReg/2.
            regMids = np.arange(midShift, dimsPerReg*nRegs, dimsPerReg)
            plt.xticks(regMids, regTits, fontsize=regFontsize, fontweight='bold')
            for tickLab, regCol in zip(plt.gca().get_xticklabels(), regCols):
                tickLab.set_color(regCol)
            #
            plt.yticks(regMids, regTits, fontsize=regFontsize, fontweight='bold')
            for tickLab, regCol in zip(plt.gca().get_yticklabels(), regCols):
                tickLab.set_color(regCol)
    #
    indMin = 0
    indMax = dimsPerReg*nRegs

    plt.xlim(indMin, dimsPerReg*nRegs)
    plt.ylim(indMin, dimsPerReg*nRegs)
    if regGrid:
        # plot inside grid
        regSeps = np.concatenate([[0], np.cumsum([dimsPerReg]*nRegs)])
        for sep in regSeps:
            plt.hlines(y=sep, xmin=0, xmax=np.size(state,0), color=gridProps['col'])
            plt.vlines(x=sep, ymin=0, ymax=np.size(state,1), color=gridProps['col'])

        # plot Border
        '''plt.axhline(y=indMin, color=gridProps['colBord'],linewidth=gridProps['lwBord'])
        plt.axhline(y=state.shape[1], color=gridProps['colBord'],linewidth=gridProps['lwBord'])
        plt.axvline(x=indMin, color=gridProps['colBord'],linewidth=gridProps['lwBord'])
        plt.axvline(x=state.shape[0], color=gridProps['colBord'],lw=gridProps['lwBord'])'''
        #
        plt.hlines(y=[indMin, indMax], xmin=indMin, xmax=indMax, color=gridProps['colBord'], lw=gridProps['lwBord'])
        plt.vlines(x=[indMin, indMax], ymin=indMin, ymax=indMax, color=gridProps['colBord'], lw=gridProps['lwBord'])







def getPhasePhaseCoh(phasex, phasey, nPhaseBins=60):
    #
    phaseBins_ = np.linspace(-np.pi, np.pi, nPhaseBins+1)
    phaseBins = np.column_stack([phaseBins_[:-1], phaseBins_[1:]])
    #
    #binDiffs.append(phaseBins_[1]-phaseBins_[0])
    phaseXY = np.zeros([nPhaseBins]*2)
    #
    cohs = []
    for bi, phaseBin in enumerate(phaseBins):
        phasexInds = np.where(np.logical_and(phasex>phaseBin[0], phasex<=phaseBin[1], ))
        #
        cohs.append(getSpikePhaseCoherence(phasey[phasexInds]))
    #
    return np.mean(cohs)


def plot_binPhaseMat(binPhaseMat, phaseCens4plot, waveCol='w', lw=2, vmin=None, vmax=None, cmap=None, 
                     yticklabels=[], xticklabels=[], cbar=False):
    #
    sb.heatmap(binPhaseMat, vmin=vmin, vmax=vmax, cmap=cmap, yticklabels=yticklabels, xticklabels=xticklabels, cbar=cbar)
    fix_sbHeatmap()
    plt.yticks(rotation=0, fontweight='bold')
    ax2 = plt.twiny()
    plotWave(phaseCens4plot, col=waveCol, lw=lw)
    ax2.set_xticks(xticklabels) #np.arange(-180, (nCy+1)*360), 180)
    ax2.grid(False)



def nanCor(a, b, method='pearson'):
    return pd.DataFrame(np.column_stack([a, b])).corr(method=method).iloc[0,1]

def plotWave(x, toDeg=True, col='gray', ls='--', lw=1, alpha=1):
    if toDeg:
        phase_rad = [convRadDeg(d,False) for d in x]
    else:
        phase_rad = x
    wave = [np.cos(r) for r in phase_rad]
    #
    ax2 = plt.twinx()
    ax2.plot(x, wave, color=col, ls=ls, lw=lw, alpha=alpha)
    ax2.grid(False)
    ax2.set_yticks([])



def getPhaseSpikes(phase, spikes, cycles=None, normCell=True, nPhaseBins=32, endCut=1250):
    #
    if cycles is None:
        endInd = len(phase)-endCut
        cyPhase = phase[:endInd]
        cySpikes = spikes[:, :endInd]
    else:
        cyPhase = []
        cySpikes = []
        for cySt, cyEn in cycles:
            if cyEn < np.size(spikes,1)-endCut:
                cyPhase.append(phase[cySt:cyEn])
                cySpikes.append(spikes[:, cySt:cyEn])
        #
        cyPhase = np.concatenate(cyPhase)
        cySpikes = np.column_stack(cySpikes)
    #
    _, phaseSpikes, phaseCentres = \
    vbf.ModulationIndex(cyPhase, cySpikes, nPhaseBins)
    #
    if normCell:
        for celli, phaseMod in enumerate(phaseSpikes):
            phaseSpikes[celli, :] = np.divide(phaseSpikes[celli, :], 
                                              float(np.sum(phaseSpikes[celli, :])))
    #
    return phaseSpikes, phaseCentres


def getRegCellInds(reg, spikeInfo, cellType='all'):
    #
    if cellType == 'all':
        cellIndsOut = [i for i, cType in enumerate(spikeInfo['des'].ravel()) if \
                       cType[1:] == reg]
    else:
        cellIndsOut = [i for i, cType in enumerate(spikeInfo['des'].ravel()) if \
                       cType == cellType+reg]
    return cellIndsOut


def run_hilbert(X):
    def next_power_of_2(n):
        '''Return next power of 2 greater than or equal to n'''
        n -= 1                 # short for "n = n - 1"
        shift = 1
        while (n+1) & n:       # the operator "&" is a bitwise operator: it compares every bit of (n+1) and n, and returns those bits that are present in both
            n |= n >> shift    
            shift <<= 1        # the operator "<<" means "left bitwise shift": this comes down to "shift = shift**2"
        return n + 1

    hilbert = sig.hilbert(X, next_power_of_2(len(X)))[range(len(X))]
    return hilbert


def get_cycles6(X, f_minmax=None, edgeID='azc', sr=sr, prop_psd=0.9, freq_filter=False):
    '''
    edgeID : str
        'azc' : ascending zero-crossings
        'dzc' : descending zero-crossings
        't' : troughs
        'p' : peaks

    RETURNS:
    cycles6, cyFreqs, cyChain
    '''
    
    edgeIDs = np.array(['p', 'dzc', 't', 'azc'])
    # sort cycle stages so first is the edgeID to detect start of cycle
    i = np.where(edgeIDs[::-1] == edgeID)[0]
    if not len(i):
        print('if edgeID is target phase, ccw to implement, otherwise edgeID entry error')
    i = i[0]
    edgeIDs_ = np.roll(edgeIDs, i+1)


    # get inds of cycle stages
    edgeInds = {}
    for edgeID in edgeIDs:
        if edgeID in ['p', 't']:
            s = [-1, 1][edgeID == 'p']
            edgeInds[edgeID] = vbf.detect_peaks(X*s, mph=0)
        else:
            if edgeID == 'azc':
                motif = [-1, 1]
            elif edgeID == 'dzc':
                motif = [1, -1]
            edgeInds[edgeID] = locateMotif(motif, np.sign(X))    


    cycles6 = []
    for i, i_next in zip(edgeInds[edgeIDs_[0]], edgeInds[edgeIDs_[0]][1:70]):
        candidate = np.repeat(-1, 6)
        candidate[0] = i
        candidate[4] = i_next

        for ei, edgeID in enumerate(edgeIDs_):
            if ei:
                try:
                    laterInds = np.where(np.sign(edgeInds[edgeID] - i) == 1)[0]
                    if edgeID == edgeIDs_[1]:
                        cand_is = [ei, ei+4]
                        cy_is = [laterInds[j] for j in range(2)]
                    else:
                        cand_is = [ei]
                        cy_is = [laterInds[0]]

                    for cand_i, cy_i in zip(cand_is, cy_is):
                        candidate[cand_i] = edgeInds[edgeID][cy_i]
                except:
                    continue

            if np.array_equal(candidate, sorted(candidate)):
                cycles6.append(candidate)

    cycles6 = np.row_stack(cycles6)
    cyFreqs = 1./(np.subtract(cycles6[:, -2], cycles6[:, 0])/sr)

    if freq_filter:
        if f_minmax is None:
            f, p = get_psd2(X, sr)
            st, en = get_inds4propAUC(p, prop_psd)
            minFreq, maxFreq = f[st], f[en]
        else:
            minFreq, maxFreq = f_minmax
        
        keepInds = np.where(np.digitize(cyFreqs, [minFreq, maxFreq])==1)[0]
        cycles6 = cycles6[keepInds]
    
    cyChain = get_cyChain(np.column_stack([cycles6[:, 0], cycles6[:, -2]]))
    
    return cycles6, cyFreqs, cyChain

def cycles6_2_cycles(cycles6):
    cycles = np.column_stack([cycles6[:, 0], cycles6[:, -2]])
    return cycles


def getCycles(trace, minFreq, maxFreq, edgeID='t', samplingRate=1250.):
    '''
    edgeID : str
        'azc' : use ascending zero-crossings
        'dzc' : use descending zero-crossings
        't' : use troughs
        'p' : use peaks

    RETURNS:
    cycles, troughs, cycleInfo
    '''
    print('depreciated. use: get_cycles6')
    #
    phase = np.angle(vbf.runHilbert(trace)) # from -pi - pi
    #
    if edgeID=='t':
        cyEdges = vbf.detect_peaks(-phase)
    elif edgeID=='p':
        cyEdges = vbf.detect_peaks(phase)
    elif edgeID=='azc':
        cyEdges = locateMotif([-1., 1.], np.sign(phase))
    elif edgeID=='dzc':
        cyEdges = locateMotif([1., -1.], np.sign(phase))
    else:
        print('error for cyle detect edgeID')
    cyclesOut = []
    troughsOut = []
    cycleInfoOut = []
    for ti, cyEdge in enumerate(cyEdges[:-1]):
        cycle = [cyEdge, cyEdges[ti+1]]
        cycDur_ms = 1000.*(cycle[1]-cycle[0])/(samplingRate)
        cycFreq = 1000./cycDur_ms
        #
        if minFreq <= cycFreq <= maxFreq:
            cyclesOut.append(cycle)
            troughsOut.append(
                np.argmin(phase[cycle[0]:cycle[1]])+cycle[0])
            cycleInfoOut.append({'dur' : cycDur_ms, 'freq' : cycFreq})
    cyclesOut = np.row_stack(cyclesOut)
    troughsOut = np.array(troughsOut)

    #
    return cyclesOut, troughsOut, cycleInfoOut

def count_cycles(X, nPhaseBins=60):
    phaseBinInds = np.digitize(get_instPhase(X), np.linspace(0, 2*np.pi, nPhaseBins+1))-1
    phi_start = phaseBinInds[0]
    if phi_start == phaseBinInds.max():
        phi_next = 0
    else:
        phi_next = phi_start+1
    inds = locateMotif([phi_start, phi_next], phaseBinInds)
    nCy = len(inds)
    if len(inds):
        nCy += len(np.unique(phaseBinInds[inds[-1]:phaseBinInds[-1]]))/float(nPhaseBins)
    else:
        nCy += len(np.unique(phaseBinInds[0:phaseBinInds[-1]]))/float(nPhaseBins)
    return nCy

def get_cyChain(cycles):
    """
    input : cycles ([nCy x 2] matrix)
    returns a 1D array of length cycles, each element denoting the number of cycles that cycle was in a chain of """
    cyPrevAmp = np.concatenate([[False], [cycles[i,0]==cycles[i-1,1] for i in np.arange(1, cycles.shape[0])]])
    cyChain = np.ones(cycles.shape[0])
    for bSt, bEn in get_boutWindows(cyPrevAmp):
        chainInds = np.arange(bSt, bEn+1)
        for i in chainInds:
            cyChain[i] = len(chainInds)
    return cyChain

def get_win4chainBout(chainBout, cycles):
    chainSt, chainEn = chainBout
    st = cycles[chainSt,0]
    en = cycles[chainEn,1]
    return st, en



def get_imfBouts(trace, amp, minFreq, maxFreq, cyAmpPercThresh=50, minChain=10, cyAmpMode='imf', sr=eegHz):
    """ get bout start and end times for a trace. 
    1. detect cycles
    2. filter cycles where amp is above cyAmpPercThresh
    3. extract bouts

    INPUT (to sort)
    cyAmpPercThresh
        float/int  | only a chain if the cycle amp is above this percThresh
    minChain
        int | min cycles in a chain AFTER removing those below cyAmpPercThresh
    cyAmpMode
        'imf' or 'barcode' : ampltude signal to filter cycles
    """
    #
    cycles = getCycles(trace, minFreq, maxFreq, edgeID='t')[0] # sort this?
    cyAmps = np.array([amp[st:en].mean() for st, en in cycles])
    selInds = np.where(cyAmps > np.percentile(cyAmps, cyAmpPercThresh))[0]
    cycles = cycles[selInds, :]
    cyAmps = cyAmps[selInds]
    cyChain = get_cyChain(cycles)
    chainBouts = get_boutWindows(cyChain > minChain)

    xChainMets = {'max_amp' : [], 
                 'mean_amp' : [],
                 'med_amp' : [],
                 'max_amp_t' : []}
    for chainBout in chainBouts:
        st, en = get_win4chainBout(chainBout, cycles)
        for f, k in zip([np.max, np.mean, np.median], ['max_amp', 'mean_amp', 'med_amp']):
            xChainMets[k].append(f(amp[st:en]))

        xChainMets['max_amp_t'].append(st+amp[st:en].argmax())

    for k in xChainMets:
        xChainMets[k] = np.array(xChainMets[k])

    return cycles, chainBouts, xChainMets





def get_linearGradient1d(x, y):
    from sklearn import linear_model
    lmodel = linear_model.LinearRegression()
    inds = np.array([i for i, y_ in enumerate(y) if not np.isnan(y_)])
    lmodel.fit(np.column_stack([x[inds]]), np.column_stack([y[inds]]))
    m = lmodel.coef_[0][0]
    return m




def plot_coupXlab(coupReg, coupImfTit, phaseCens4plot, nCy, fontweight='bold', fontsize=12, 
                  lw=2, addWave=True, xdeg=90, addXlab=True):
    coupTit, coupCol = getTitColForReg(coupReg)
    plt.xticks(np.arange(0, 360*nCy+180, xdeg), fontsize=fontsize)
    if coupImfTit == '4Hz':
        coupImfTit = '4-Hz'
        
    if addXlab:
        plt.xlabel(coupTit+' '+coupImfTit+' phase (deg.)', color=coupCol, fontweight=fontweight, fontsize=fontsize)
    if addWave:
        plotWave(phaseCens4plot, nCy, col=coupCol, lw=lw)

def get_cyPhaseMod(x, phase, cycles, nPhaseBins, phaseMin=0, phaseMax=2*np.pi, binFunc=np.mean):
    #
    '''currently only works if x is 1-d array

    RETURNS:
    cyPhaseMod, phaseCens
    '''
    phaseBins_ = np.linspace(phaseMin, phaseMax, nPhaseBins+1)
    phaseBinInds = np.digitize(phase, phaseBins_)-1
    phaseCens = edges2cens(phaseBins_)
    #
    cyPhaseMod = []
    nDim = np.ndim(x)
    if nDim == 1:
        for cySt, cyEn in cycles:
            cyX = x[cySt:cyEn]
            cyPhaseInds = phaseBinInds[cySt:cyEn]
            phaseModi = np.repeat(np.nan, nPhaseBins)
            for pi in range(nPhaseBins):
                phaseModi[pi] = binFunc(cyX[np.where(cyPhaseInds==pi)[0]])
            #
            cyPhaseMod.append(phaseModi)
        #
        cyPhaseMod = np.row_stack(cyPhaseMod)
    elif nDim == 2:
        for cySt, cyEn in cycles:
            cyX = x[:, cySt:cyEn]
            cyPhaseInds = phaseBinInds[cySt:cyEn]
            phaseModi = np.row_stack([np.repeat(np.nan, nPhaseBins) for di in range(x.shape[0])])
            for pi in range(nPhaseBins):
                phaseModi[:, pi] = binFunc(cyX[:, np.where(cyPhaseInds==pi)[0]], axis=1)
            #
            cyPhaseMod.append(phaseModi)
    else:
        for cySt, cyEn in cycles:
            cyX = x.T[cySt:cyEn].T
            cyPhaseInds = phaseBinInds[cySt:cyEn]

            phaseModi = np.array([np.full_like(x.T[0], np.nan) for pi in range(nPhaseBins)])

            for pi in range(nPhaseBins):
                phaseModi[pi] = cyX.T[np.where(cyPhaseInds==pi)[0]].mean(axis=0)
            cyPhaseMod.append(phaseModi.T)
        #
        #
    cyPhaseMod = np.array(cyPhaseMod)
    #
    return cyPhaseMod, phaseCens


def get_zoneCyPhaseMod(strengths, coupPhase, coupCycles, visits, nPhaseBins, seshLen=None, zone_zScore=False,
                       phaseMin=0., phaseMax=2*np.pi):
    #
    zCyPhaseMod = {}
    cyTimes = np.mean(coupCycles, axis=1, dtype=int)
    for zone in zones:
        if seshLen is None:
            zTimes = np.concatenate([np.arange(en, ex) for en, ex in visits[zone]])
        else:
            zTimes = np.concatenate([np.arange(en, ex) for en, ex in \
                                     visits[zone] if en > seshLen-int(20*60*eegHz)])
        if zone_zScore:
            if np.ndim(strengths) == 3:
                strengths4mod = np.zeros_like(strengths)
                for x in range(strengths.shape[0]):
                    for y in range(strengths.shape[1]):
                        m = np.mean(strengths[x, y, zTimes])
                        sd = np.std(strengths[x, y, zTimes])
                        strengths4mod[x, y, :] = (strengths4mod[x, y, :]-m)/sd
            elif np.ndim(strengths) == 2:
                m = np.mean(strengths[:, zTimes], axis=1)
                sd = np.std(strengths[:, zTimes], axis=1)
                strengths4mod = np.row_stack([(strengthi-mi)/sdi for strengthi, mi, sdi in zip(strengths, m, sd)])
        else:
            strengths4mod = strengths

        zCyInds = np.where(np.in1d(cyTimes, zTimes))[0]
        zCycles = coupCycles[zCyInds, :]
        #
        zCyPhaseMod[zone] = get_cyPhaseMod(strengths4mod, coupPhase, zCycles, nPhaseBins, 
                                           phaseMin=phaseMin, phaseMax=phaseMax)[0]
    return zCyPhaseMod

def get_xPhaseMod(phase, x, nPhaseBins=32, f=np.mean, fe=stats.sem):
    
    '''
    RETURNS:
    phaseCens, xPhaseMod_M, xPhaseMod_E
    '''
    edges = np.linspace(0, 2*np.pi, nPhaseBins+1)

    phaseBinInds = np.digitize(phase, edges)-1

    ndim = np.ndim(x) 
    if ndim == 1:
        xPhaseMod_M = np.zeros(nPhaseBins)
        xPhaseMod_E = np.zeros(nPhaseBins)
    elif ndim == 2:
        xPhaseMod_M = np.zeros([x.shape[0], nPhaseBins])
        xPhaseMod_E = np.zeros([x.shape[0], nPhaseBins])
    else:
        raise TypeError('x must be 1- or 2-D')

    for bi in range(nPhaseBins):
        inds = np.flatnonzero(phaseBinInds == bi)
        if ndim == 1:
            m, se = [f_(x[inds]) for f_ in [f, fe]]
            xPhaseMod_M[bi] = m
            xPhaseMod_E[bi] = se
        else:
            m, se = [f_(x[:, inds], axis=1) for f_ in [f, fe]]
            xPhaseMod_M[:, bi] = m
            xPhaseMod_E[:, bi] = se
    phaseCens = edges2cens(edges)
    
    return phaseCens, xPhaseMod_M, xPhaseMod_E



def getBoutTimes(boutStatus):
        #
        boutStartInds = list(np.array(locateMotif([False, True], boutStatus))+1)
        boutEndInds = list(np.array(locateMotif([True, False], boutStatus)))

        if not all([len(boutStartInds), len(boutEndInds)]):
            if len(boutStartInds):
                bouts = np.column_stack([boutStartInds[0], boutStatus.shape[0]-1])
            elif len(boutEndInds):
                bouts = np.column_stack([0, boutEndInds[0]])
            else:
                if boutStatus.sum()==len(boutStatus):
                    bouts = np.column_stack([0, (len(boutStatus)-1)])
                else:
                    bouts = None
        else:
            if all([len(boutStartInds)==1, len(boutEndInds)==1, boutEndInds[0] < boutStartInds[0]]):
                bouts = None
            else:
                if boutEndInds[0] < boutStartInds[0]:
                    boutEndInds = boutEndInds[1:]
                if boutStartInds[-1] > boutEndInds[-1]:
                    #if eegTrack is not None:
                        #boutEndInds.append(np.size(eegTrack,0)-1)
                    #else:
                    boutEndInds.append(boutStatus.shape[0]-1)
                #
                #if eegTrack is None:
                bouts = np.column_stack([boutStartInds, boutEndInds])
                #else:
                    #boutStartTimes = [ int(round(eegTrack['t'].iloc[ind])) for ind in boutStartInds ]
                    #boutEndTimes = [ int(round(eegTrack['t'].iloc[ind])) for ind in boutEndInds ]
            #

        return bouts

def view_image(filename):
    try:
        from wand.image import Image as WImage
        img = WImage(filename=filename)
        return img
    except:
        print('must be in directory where filename exists')



def get_nTrans4ideal(nICs, kLabs0, classBsnms):

    bsnms = np.setdiff1d(np.unique(classBsnms), ['all'])

    N_TRANSITIONS = 0 # will count the number of transitions u need to make for 'ideal' clustering

    k_avStateis0 = [] # find the original cluster indices assigned to the AV-STATES by kmeans
    for sti in range(nICs):
        k_avStateis0.append(np.intersect1d(np.where(kLabs0==sti)[0], np.where(classBsnms=='all')[0]))

    ''' first find the all possible transfer combinations - 
    only will move avstates if they are not in an ideal kluster  '''

    stis_2trans=[] # the n x [ki, stj]; where ki is the original cluster it was assigned and stj the state index corresponsind to classXXXX
    empty_kiis=[] # the klusters which can accept an avState
    k_avStateis = {}
    for ki in range(nICs):
        if len(k_avStateis0[ki])==1:
            k_avStateis[ki] = k_avStateis0[ki][0]
        elif len(k_avStateis0[ki])>1:
            #k_avStateis[ki] = k_avStateis0[ki][0]
            for stj in k_avStateis0[ki]:
                stis_2trans.append([ki, stj]) # these are the states to be moved
        else:
            empty_kiis.append(ki)
    try:
        stis_2trans = np.row_stack(stis_2trans)
    except:
        return N_TRANSITIONS
    #
    #
    transDestinationCombos = list(itertools.permutations(empty_kiis, len(empty_kiis)))
    ovPop_ks = np.unique(stis_2trans[:,0])
    ovPop_assignInds=[]
    for ki in ovPop_ks:
        ovPop_assignInds.append(np.where(stis_2trans[:,0]==ki)[0])
    #
    #
    n2reorder = len(np.concatenate(ovPop_assignInds))
    transPerms_ = list(itertools.permutations(np.arange(n2reorder), n2reorder))
    transPerms = [] # all the combos of leaving 1 in clu - then for each of these
                    # generate additional assign combos
    for perm in transPerms_:
        permOK = [] # make sure 1 of the ovStates remains in its original assigned cluster
        for ovi, ovPop_assignIndsi in enumerate(ovPop_assignInds):
            permOK.append([perm[i] in ovPop_assignIndsi for i in ovPop_assignIndsi])
        if all(np.concatenate(permOK)):
            transPerms.append(perm)
    #
    #
    transPerm = transPerms[0] # <-------- only need to consider one perm for basic method
    transPerm_clu=[]
    for ovPop_assignIndsi in ovPop_assignInds:
        cluPerm = []
        for i in ovPop_assignIndsi:
            cluPerm.append(transPerm[i])
        #
        transPerm_clu.append(cluPerm)

    #
    inds2trans = []
    for ki, transPerm_clui in enumerate(transPerm_clu):
        inds2trans.append(np.setdiff1d(ovPop_assignInds[ki], [transPerm_clui[0]]))

    inds2trans = np.concatenate(inds2trans)
    #
    kLabs = np.copy(kLabs0)
    for stii2trans, k2send in zip(inds2trans, empty_kiis):
        classInd = stis_2trans[stii2trans, 1]
        kLabs[classInd] = k2send
        #kLabs
        #
        #kTransCount=0 <-- if diff permutations did make a difference, would need to count each here
        for ki in range(nICs):
            kInds = np.where(kLabs==ki)[0]
            np.setdiff1d(classBsnms[kInds], ['all'])
            #
            bsnmsInKlu = np.setdiff1d(classBsnms[kInds], ['all'])
            N_TRANSITIONS += np.sum([abs(len(np.where(bsnm==bsnmsInKlu)[0])-1) for bsnm in bsnms])
    #
    return N_TRANSITIONS




def send_email2( receivers, subject, body, filenames=[], text_type='plain', images=[]):

        from email.mime.multipart import MIMEMultipart
        from email.utils import COMMASPACE,formatdate
        from email.mime.text import MIMEText
        import smtplib

        message = MIMEMultipart()
        message['From'] = 'spikedata0@gmail.com'
        message['To'] = COMMASPACE.join(receivers)
        message['Date'] = formatdate(localtime=True)
        message['Subject'] = subject
        message.attach(MIMEText(body, text_type, 'utf-8'))
        if filenames:
                from email.mime.base import MIMEBase
                from email import encoders
                import ntpath
                for filename in filenames:
                        part = MIMEBase('application', 'base64')
                        part.set_payload(open(filename, "rb").read())
                        encoders.encode_base64(part)
                        part.add_header('Content-Disposition', 'attachment; filename=%s' % ntpath.basename(filename))
                        message.attach(part)
        if images:
                from email.MIMEImage import MIMEImage
                for i in images:
                        fp = open(i[0], 'rb')
                        msgImage = MIMEImage(fp.read(), _subtype="jpeg+svg+png")
                        fp.close()
                        msgImage.add_header('Content-ID', i[1])
                        message.attach(msgImage)

        server = smtplib.SMTP('smtp.gmail.com', 587)
        server.ehlo()
        server.starttls()
        server.login('spikedata0@gmail.com'  , 'SouIff4927')
        server.sendmail('spikedata0@gmail.com'  , receivers, message.as_string())
        server.quit() 


import ccw_tailorEMD_functions as te
