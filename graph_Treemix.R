
#Create graph in Treemix

#OUTPUT tree and residual graph

###################################################################

args <- commandArgs(trailingOnly = TRUE)

folder=args[1]
ord_pop=args[2]

#folder="/home/marco/pool_first/Treemix_new2/inter_no_mig/t_1"
#ord_pop="/home/marco/pool_first/Treemix_new2/end_ord"

source("plotting_funcs.R")
png(width=1000,height=800,file=paste(folder,"_tree.png",sep=""))
plot_tree(folder)
dev.off()

png(width=1000,height=800,file=paste(folder,"_resi.png",sep=""))
plot_resid(folder,ord_pop)
dev.off()



#############################################################