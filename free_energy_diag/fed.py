#| - IMPORT Modules


#__|

class FreeEnergyDiagram:

	def __init__(self, reaction):
		self.reaction = reaction
		fig = pylab.figure()
		fig.set_facecolor("white")
		figlegend = pylab.figure(figsize=(3,2))
		ax = fig.add_subplot(111)


	# # def plot_ideal_
    #
    #
	# #| - Incorporated &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	# # # plt.figure().patch.set_facecolor("white")
	# # fig = pylab.figure()
	# # fig.set_facecolor("white")
	# # figlegend = pylab.figure(figsize=(3,2))
	# # ax = fig.add_subplot(111)
	# #__|
    #
    #
	# #| - Additional Plotting                                                    2222
	# if plot_ideal_steps == True:
	# 	if reaction == "2e":
	# 		# plt.plot(xdata, ideal_ydata_2e, color = grey)
	# 		ax.plot(xdata, ideal_ydata_2e, color = blue, label = "dop")
    #
	# 	if reaction == "4e":
	# 		# plt.plot(xdata, ideal_ydata_4e, color = grey)
	# 		ax.plot(xdata, ideal_ydata_4e, color = blue, label = "dop")
    #
	# #__|                                                                        2222
    #
    #
	# #| - Plotting System Data                                                   2222
	# data = hib_cu_b_yc
	# ax.plot(xdata, data['ydata_eq'], color = blue, label = "dop")
    #
	# data = hib_cu_m_yc
	# ax.plot(xdata, data['ydata_eq'], color = red, label = "goople")
    #
	# data = hib_ni_b
	# ax.plot(xdata, data['ydata_eq'], color = magenta, label = "goople")
    #
	# # plt.plot(xdata,data["ydata_nobias"],
	# # color = blue, label="hib_ni_b|OP: "+str(round(data["overpotential"],2)))
    #
	# # data = hib_cu_b_yc
	# # plt.plot(xdata,data["ydata_nobias"],
	# # color = red, label="hib_cu_b_yc|OP: "+str(round(data["overpotential"],2)))
	# #
	# #
	# # data = hib_cu_m_yc
	# # plt.plot(xdata,data["ydata_nobias"],
	# # color = green, label="hib_cu_m_yc|OP: "+str(round(data["overpotential"],2)))
	# #
	# # data = hib_ni_b
	# # plt.plot(xdata,data["ydata_nobias"],
	# # color = blue, label="hib_ni_b|OP: "+str(round(data["overpotential"],2)))
    #
	# #__|                                                                        2222
    #
    #
	# #| - Format Plot                                                            2222
    #
	# #| - TMP
	# # plt.legend(fontsize=16, loc="best").draggable()
	# # plt.ylabel("Gibbs Free Energy [eV]",fontsize=18)
	# # plt.xlabel("Reaction Coordinate", fontsize=18)
	# #
	# # if reaction 	== "2e":
	# # 	rxn_coord_lab = ["$O_2$","*OOH", "$H_2O_2$"]
	# # 	plt.xticks(np.array(range(3))+0.5, rxn_coord_lab,rotation=35,fontsize=18)
	# # elif reaction 	== "4e":
	# # 	rxn_coord_lab = ["$O_2$","*OOH", "*O + $H_2O$","*OH + $H_2O$","$2H_2O$"]
	# # 	plt.xticks(np.array(range(5))+0.5, rxn_coord_lab,rotation=35,fontsize=18)
	# #
	# # plt.tick_params(labelsize=18)
	# # plt.tight_layout()
	# #__|
    #
	# handles, labels = ax.get_legend_handles_labels()
	# figlegend.legend(handles, labels)
    #
	# # ax.set_title('foo')
	# ax.set_xlabel("Reaction Coordinate", fontsize=18)
	# ax.set_ylabel("Gibbs Free Energy [eV]",fontsize=18)
	# ax.tick_params(labelsize=20, bottom="off", top="off", right="off")
    #
	# if reaction 	== "2e":
	# 	rxn_coord_lab = ["$O_2$","*OOH", "$H_2O_2$"]
	# 	ax.set_xticks(np.array(range(3))+0.5)
	# 	ax.set_xticklabels(rxn_coord_lab,rotation=35,fontsize=18)
	# elif reaction 	== "4e":
	# 	rxn_coord_lab = ["$O_2$","*OOH", "*O + $H_2O$","*OH + $H_2O$","$2H_2O$"]
	# 	ax.set_xticks(np.array(range(5))+0.5)
	# 	ax.set_xticklabels(rxn_coord_lab,rotation=35,fontsize=18)
	# #__|                                                                        2222
    #
    #
	# #| - Output Figure (Saving and Displaying)                                  2222
	# # plt.savefig("fig_fed.svg", format="svg",spi=1200,transparent=False)
	# # plt.show()
    #
	# # fig.show()
	# # figlegend.show()
    #
	# fig.savefig('figure.svg')
	# figlegend.savefig('legend.svg')
    #
	# #__|                                                                        2222
