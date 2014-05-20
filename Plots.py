

figure = plt.Figure()							#Figure needed to plot cbar
axes = figure.add_subplot(111)			#For the heatmap
axes.cla()										#Clear if anything was left. I don't know why I need this, but it improves functionality
xticks = np.linspace(min(X), max(X), 4)
yticks = np.linspace(min(Y), max(Y), 4)
cticks = np.linspace(min(,cmax,numctics)

		axes.set_xlim([xmin, xmax])
		axes.set_ylim([ymin, ymax])

		axes.xaxis.set_ticks(xticks)
		axes.yaxis.set_ticks(yticks)

		axes.tick_params(axis='x',which='both',bottom='on',top='off',direction='out',width=1.5)
		axes.tick_params(axis='y',which='both',left='on',right='off',direction='out',width=1.5)

		axes.set_xticklabels(xticks,fontweight='bold',fontsize=fontsize)
		axes.set_yticklabels(yticks,fontweight='bold',fontsize=fontsize)

		axes.set_xlabel("X-plane / %i nm" % int(aval), fontweight='bold',fontsize=fontsize)
		axes.set_ylabel("Y-plane / %i nm" % int(aval), fontweight='bold',fontsize=fontsize)

		t = axes.set_title(r'$\mathbf{|\vec{E}|^2/ |\vec E_0|^2}$', fontweight='bold',fontsize=fontsize+1)
		t.set_y(1.04) 							#Offset the title so it doesn't overlap with the plot
		im = axes.imshow(data,cmap=color,origin='lower',alpha=0.9,extent=[xmin,xmax,ymin,ymax])
		im.set_clim(cmin,cmax) 		#Set the colorbar limits
		cbar = figure.colorbar(im, ticks = cticks)
		cbar.ax.set_yticklabels(cticks,fontweight='bold',fontsize=fontsize)

		canvas.draw() #Refreshed the canvas so that plot appears

		#Need to reset these to get the plot to diplay and resize properly
		sizer = wx.BoxSizer(wx.HORIZONTAL)
		sizer.Add(canvas, 1, wx.ALL | wx.EXPAND | wx.ALIGN_CENTER)
		SetSizer(sizer)
		Fit()
		GetParent().ResetSash()
		cidpress = canvas.mpl_connect('button_press_event',OnPress)
		cidrelease = canvas.mpl_connect('button_release_event',OnRelease)
