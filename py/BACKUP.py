
w_lda, w_psd = 7, 5
wTot = w_lda + w_psd
h = 4
hTot = nChannels * h

facecolor = None

plt.figure(figsize=(wTot, hTot))
grid = plt.GridSpec(hTot, wTot, hspace=3.5, wspace=0.5)

reload_functions()

currH = 0
for channel in channels:

    timeofday_colors = {'morning' : 'gray', 'evening' : chColors[channel]}

    currW = 0
    obs_, obsInfo_ = get_obs_4_channel(channel)

    lda = LDA()
    lda.fit(obs_, obsInfo_['timeofday'])

    plt.subplot(grid[currH:(currH+h), currW:(currW+w_lda)], facecolor=facecolor)
    flo.plot_title(channel, color=chColors[channel], fontsize=16)

    obProj_ = np.array([np.dot(lda.coef_, ob) for ob in obs_])
    rdmY = np.array([np.random.random() for _ in obProj_])

    for timeofday in timeofdays:
        inds = np.flatnonzero(obsInfo_['timeofday'] == timeofday)
        plt.plot(obProj_[inds], rdmY[inds], 'o', color=timeofday_colors[timeofday], label=timeofday)
        plt.ylim(-3, 3)
    plt.legend()
    plt.xlabel('Projection (a.u.))')
    currW += w_lda

    plt.subplot(grid[currH:(currH+h), currW:(currW+w_lda)], facecolor=facecolor)
    lw = 2
    for timeofday in timeofdays:
        color = timeofday_colors[timeofday]
        inds = np.flatnonzero(obsInfo_['timeofday'] == timeofday)
        m, se = flo.get_nanMSE4obs(obs_[inds, :])
        plt.plot(freqAx, m, color=color, lw=lw)
        plt.fill_between(freqAx, m-se, m+se, color=color, alpha=0.5)
    plt.ylabel('Power ()', color=timeofday_colors['evening'])
    plt.xscale('log')
    plt.yscale('log')
    ax2 = plt.twinx()
    yMax = np.abs(lda.coef_[0]).max()*1.8
    ax2.plot(freqAx, lda.coef_[0], color='k')
    ax2.set_ylim(-yMax, yMax)
    ax2.set_ylabel('Coef.')

    currH += h

flo.save_image('jjLDA_PSD')
