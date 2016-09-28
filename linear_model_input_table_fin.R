################################################################################################################################
######
##create input for linear model



pops_sel=c("p07C", "p07D", "p07E", "p07F", "p07G", "p07H", "p07J", "p07K", "p07L", "p07M", "p07N", "p07O", "p07P", "p07Q", "p07R", "p11A", "p11AA", "p11AB", "p11AC", "p11AE", "p11AG", "p11AH", "p11AJ", "p11B", "p11C", "p11D", "p11E", "p11F", "p11G", "p11H", "p11J", "p11K", "p11L_F0", "p11M", "p11N", "p11O", "p11P", "p11Q", "p11R", "p11S", "p11T", "p11U", "p11V", "p11W", "p11X", "p11Y", "p11Z", "p14A", "p14B", "p14C", "p14D", "p14E")

nome_out="/home/marco/Dropbox/dottorato_Neuchatel/npstat_output_fin/tab_lm_1_0_fin"


tab_tot=read.table("/home/marco/Dropbox/dottorato_Neuchatel/npstat_output_fin/sum_tot_varscan_1_0")
tab_inter=read.table("/home/marco/Dropbox/dottorato_Neuchatel/npstat_output_fin/sum_intergenic_1000bp_varscan_1_0")
tab_intro=read.table("/home/marco/Dropbox/dottorato_Neuchatel/npstat_output_fin/sum_intron_varscan_1_0")
tab_cds=read.table("/home/marco/Dropbox/dottorato_Neuchatel/npstat_output_fin/sum_CDS_varscan_1_0")


tree_in="/home/marco/pool_first/Treemix_tot_05_2/tot_05_m7/t_49.treeout"



tab_tot=tab_tot[pops_sel,]
tab_inter=tab_inter[pops_sel,]
tab_intro=tab_intro[pops_sel,]
tab_cds=tab_cds[pops_sel,]

tab_ms=read.table("/home/marco/Dropbox/dottorato_Neuchatel/npstat_output_fin/SummaryStats_160824",header = T)
rownames(tab_ms)=paste("p",tab_ms[,1],sep="")
tab_ms=tab_ms[pops_sel,]


tab52=read.table("/home/marco/Dropbox/dottorato_Neuchatel/npstat_output_fin/Alyrata_P52",header = TRUE,sep="\t")
tab52=tab52[pops_sel,]





##calculate distance (like Pip)
library(PBSmapping)

# Function: Convert degrees to radians: http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
deg2rad <- function(deg) return(deg*pi/180)
#
# Function: Calculates the geodesic distance between two points specified by radian latitude/longitude using
# Vincenty inverse formula for ellipsoids (vif): http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
gcd.vif <- function(long1, lat1, long2, lat2) {
  
  # WGS-84 ellipsoid parameters
  a <- 6378137         # length of major axis of the ellipsoid (radius at equator)
  b <- 6356752.314245  # ength of minor axis of the ellipsoid (radius at the poles)
  f <- 1/298.257223563 # flattening of the ellipsoid
  
  L <- long2-long1 # difference in longitude
  U1 <- atan((1-f) * tan(lat1)) # reduced latitude
  U2 <- atan((1-f) * tan(lat2)) # reduced latitude
  sinU1 <- sin(U1)
  cosU1 <- cos(U1)
  sinU2 <- sin(U2)
  cosU2 <- cos(U2)
  
  cosSqAlpha <- NULL
  sinSigma <- NULL
  cosSigma <- NULL
  cos2SigmaM <- NULL
  sigma <- NULL
  
  lambda <- L
  lambdaP <- 0
  iterLimit <- 100
  while (abs(lambda-lambdaP) > 1e-12 & iterLimit>0) {
    sinLambda <- sin(lambda)
    cosLambda <- cos(lambda)
    sinSigma <- sqrt( (cosU2*sinLambda) * (cosU2*sinLambda) +
                        (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda) )
    if (sinSigma==0) return(0)  # Co-incident points
    cosSigma <- sinU1*sinU2 + cosU1*cosU2*cosLambda
    sigma <- atan2(sinSigma, cosSigma)
    sinAlpha <- cosU1 * cosU2 * sinLambda / sinSigma
    cosSqAlpha <- 1 - sinAlpha*sinAlpha
    cos2SigmaM <- cosSigma - 2*sinU1*sinU2/cosSqAlpha
    if (is.na(cos2SigmaM)) cos2SigmaM <- 0  # Equatorial line: cosSqAlpha=0
    C <- f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
    lambdaP <- lambda
    lambda <- L + (1-C) * f * sinAlpha *
      (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))
    iterLimit <- iterLimit - 1
  }
  if (iterLimit==0) return(NA)  # formula failed to converge
  uSq <- cosSqAlpha * (a*a - b*b) / (b*b)
  A <- 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
  B <- uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
  deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM^2) -
                                             B/6*cos2SigmaM*(-3+4*sinSigma^2)*(-3+4*cos2SigmaM^2)))
  s <- b*A*(sigma-deltaSigma) / 1000
  
  return(s) # Distance in km
}

# Calculates the geodesic distance between two points specified by degrees (DD) latitude/longitude using
# Vincenty inverse formula for ellipsoids (vif)
gcd <- function(long1, lat1, long2, lat2) {
  
  # Convert degrees to radians
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)
  
  return(list( dist = gcd.vif(long1, lat1, long2, lat2)) )
}


##get parent node
getParent<-function(tree,node){
  
  pare<-tree$edge[which(tree$edge[,2]==node),1]
  return(pare)
}



eastsubset <- subset(tab52, GeneticCluster=="E", select=c(x, y))
easthull <- chull(x=eastsubset$x, y=eastsubset$y)
easthullpoints <- eastsubset[easthull,]
easthullpoints$PID <- 1
easthullpoints$POS <- c(1,2,3,4,5,6,7)
names(easthullpoints) [names(easthullpoints)=="x"] <- "X"
names(easthullpoints) [names(easthullpoints)=="y"] <- "Y"
easthullcentroidPre <- as.PolySet(easthullpoints, projection = "LL", zone = NULL)
easthullcentroid <- calcCentroid(easthullcentroidPre, rollup=1)
#easthullcentroid # 1 -78.00514 39.44167

westsubset <- subset(tab52, GeneticCluster=="W", select=c(x, y))
westhull <- chull(x=westsubset$x, y=westsubset$y)
westhullpoints <- westsubset[westhull,]
westhullpoints$PID <- 2
westhullpoints$POS <- c(1,2,3,4,5,6,7,8)
names(westhullpoints) [names(westhullpoints)=="x"] <- "X"
names(westhullpoints) [names(westhullpoints)=="y"] <- "Y"
westhullcentroidPre <- as.PolySet(westhullpoints, projection = "LL", zone = NULL)
westhullcentroid <- calcCentroid(westhullcentroidPre, rollup=1)
#westhullcentroid # -88.15296 44.1265

# For each population, calculate distance from the two centroids of polygon hulls
centroids <- rbind(easthullcentroid, westhullcentroid) 

gcd_centroid=rep(NA,dim(tab52)[1])


for(i in 1:dim(tab52)[1]) {
  if(tab52$GeneticCluster[i]=="E") {
    gcd_centroid[i]=unlist(gcd(tab52$x[i], tab52$y[i], centroids[1,"X"], centroids[1,"Y"]))
  } else {
    gcd_centroid[i]=unlist(gcd(tab52$x[i], tab52$y[i], centroids[2,"X"], centroids[2,"Y"]))
  }

  
}


###################################################
###Adimixture
Admix=rep("no",dim(tab52)[1])
names(Admix)=rownames(tab52)
Admix[c("p11H","p11G","p11X","p14B","p14D","p14E","p07H")]="yes"

AdmixWE=rep("no",dim(tab52)[1])
names(AdmixWE)=rownames(tab52)
AdmixWE[c("p14B","p14D","p14E","p07H")]="yes"



###################################################################################################################################
Cluster3=as.character(tab52[,"GeneticCluster"])
names(Cluster3)=rownames(tab52)
#Cluster3[c("p11A","p11B","p11G")]="SW"
Cluster3[c("p11A","p11B")]="SW"

#ice sheet

Ice=rep("low",dim(tab52)[1])
names(Ice)=rownames(tab52)
bla=names(Ice)[which(tab52$GeneticCluster=="W")]
pop_old=c("p11X","p11H","p11G","p11A","p11B","p11W")
Ice[setdiff(bla,pop_old)]="high"
Ice[c("p11V","p11U","p07G","p07F","p07D", "p07E")]="high"

########################################################

######################
# geographic distance based on topology
library(ape)
require(phytools)



###West

##calculation coordinate each Cluster separate
aa=read.tree(paste(tree_in,"_noout",sep=""))

tab=cbind(tab52[,c("x","y")],Cluster3)
colnames(tab)=c("lon","lat","Cluster")
tab=tab[aa$tip.label,]

clu_sel="W"
pop_rem=setdiff(rownames(tab),rownames(tab)[which(tab$Cluster==clu_sel)])

aa=drop.tip(aa, pop_rem)
aa$node.labels=(length(aa$tip.label)+1):(((length(aa$tip.label)+1)+Nnode(aa,internal.only = TRUE))-1)

plot.phylo(aa,show.node.label = T)
anc_node=30

coo=tab[which(tab$Cluster==clu_sel),c("lat","lon")]
rownames(coo)=aa$tip.label
xx<-phylo.to.map(aa,coo,plot=FALSE,xlim=c(min(tab$lon)-0.5,max(tab$lon)+0.5), ylim=c(min(tab$lat)-0.5,max(tab$lat)+0.5))
#plot(xx,type="direct")
bb=phylomorphospace(aa,xx$coords[,2:1])
dev.off()

pro=c(rownames(tab)[which(tab$Cluster==clu_sel)],as.character(aa$node.labels))
tab_coo=cbind(bb$xx,bb$yy)
rownames(tab_coo)=pro

tab_coo2=tab_coo[anc_node:dim(tab_coo)[1],]
rownames(tab_coo2)=75:102 

###East

aa=read.tree(paste(tree_in,"_noout",sep=""))

clu_sel="E"
pop_rem=setdiff(rownames(tab),rownames(tab)[which(tab$Cluster==clu_sel)])

aa=drop.tip(aa, pop_rem)
aa$node.labels=(length(aa$tip.label)+1):(((length(aa$tip.label)+1)+Nnode(aa,internal.only = TRUE))-1)
plot.phylo(aa,show.node.label = TRUE)
anc_node=22
coo=tab[which(tab$Cluster==clu_sel),c("lat","lon")]
rownames(coo)=aa$tip.label
xx<-phylo.to.map(aa,coo,plot=FALSE,xlim=c(min(tab$lon)-0.5,max(tab$lon)+0.5), ylim=c(min(tab$lat)-0.5,max(tab$lat)+0.5))

bb=phylomorphospace(aa,xx$coords[,2:1])
dev.off()

pro=c(rownames(tab)[which(tab$Cluster==clu_sel)],as.character(aa$node.labels))
tab_coo=cbind(bb$xx,bb$yy)
rownames(tab_coo)=pro


tab_coo1=tab_coo[anc_node:dim(tab_coo)[1],]
rownames(tab_coo1)=55:74

colnames(tab_coo1)=colnames(tab_coo2)=c("lon","lat")
tab_coo_f=rbind(tab_coo1,tab_coo2)

write.table(tab_coo_f,file=paste(nome_out,"_nodes",sep=""))

##################################
##calculation distance from each pops

tab_coo=rbind(tab[,1:2],tab_coo_f)

aa=read.tree(paste(tree_in,"_noout",sep=""))
aa$node.labels=53:((53+Nnode(aa,internal.only = TRUE))-1)
plot.phylo(aa,show.node.label = TRUE)

##distance from 77
anc_node=77
clu_sel="W"
seli=setdiff(rownames(tab)[which(tab$Cluster==clu_sel)],c("p11G","p11X","p11H"))
gcd_tree_anc=rep(0,length(seli))
names(gcd_tree_anc)=seli


for(i in 1:length(gcd_tree_anc)) {
  pap=getParent(aa,which(aa$tip.label==names(gcd_tree_anc)[i]))
  gcd_tree_anc[names(gcd_tree_anc)[i]]=unlist(gcd(tab_coo[names(gcd_tree_anc)[i],1], tab_coo[names(gcd_tree_anc)[i],2], tab_coo[as.character(pap),1], tab_coo[as.character(pap),2]))
  pap_pre=pap
  if(pap==anc_node) {
    flag=FALSE
  } else {
    flag=TRUE
  }
  
  while(flag) {
    if(pap==anc_node) {
      flag=FALSE
    } else {
      pap=getParent(aa,pap)  
      dista=unlist(gcd(tab_coo[as.character(pap_pre),1], tab_coo[as.character(pap_pre),2], tab_coo[as.character(pap),1], tab_coo[as.character(pap),2]))
      if(length(dista)==0) { #per rami a tre
        dista=0
      }
      gcd_tree_anc[names(gcd_tree_anc)[i]]=gcd_tree_anc[names(gcd_tree_anc)[i]]+dista
      pap_pre=pap
    }

  }
}

gcd_tree_anc_W=gcd_tree_anc


##distance from 56
anc_node=56
clu_sel="E"
seli=setdiff(rownames(tab)[which(tab$Cluster==clu_sel)],c("p11O"))
gcd_tree_anc=rep(0,length(seli))
names(gcd_tree_anc)=seli


for(i in 1:length(gcd_tree_anc)) {
  pap=getParent(aa,which(aa$tip.label==names(gcd_tree_anc)[i]))
  gcd_tree_anc[names(gcd_tree_anc)[i]]=unlist(gcd(tab_coo[names(gcd_tree_anc)[i],1], tab_coo[names(gcd_tree_anc)[i],2], tab_coo[as.character(pap),1], tab_coo[as.character(pap),2]))
  pap_pre=pap
  if(pap==anc_node) {
    flag=FALSE
  } else {
    flag=TRUE
  }
  
  while(flag) {
    if(pap==anc_node) {
      flag=FALSE
    } else {
      pap=getParent(aa,pap)  
      dista=unlist(gcd(tab_coo[as.character(pap_pre),1], tab_coo[as.character(pap_pre),2], tab_coo[as.character(pap),1], tab_coo[as.character(pap),2]))
      if(length(dista)==0) { #per rami a tre
        dista=0
      }
      gcd_tree_anc[names(gcd_tree_anc)[i]]=gcd_tree_anc[names(gcd_tree_anc)[i]]+dista
      pap_pre=pap
    }
    
  }
}

#sort(gcd_tree_anc)
gcd_tree_anc_E=gcd_tree_anc

gcd_tree_anc=rep(NA,dim(tab)[1])
names(gcd_tree_anc)=rownames(tab)
gcd_tree_anc[names(gcd_tree_anc_E)]=gcd_tree_anc_E
gcd_tree_anc[names(gcd_tree_anc_W)]=gcd_tree_anc_W

gcd_tree_anc=gcd_tree_anc[pops_sel]

###################################################################################################################################


gcd_anc=rep(0,52)
names(gcd_anc)=pops_sel

aa=read.tree(paste(tree_in,"_noout",sep=""))
aa$node.labels=53:((53+Nnode(aa,internal.only = TRUE))-1)
plot.phylo(aa,show.node.label = T)


tab=tab52[,c("x","y","GeneticCluster")]
colnames(tab)=c("lon","lat","Cluster")
##West
sel=rownames(tab[which(tab$Cluster=="W"),])
anc="77"

for(i in 1:length(sel)) {
  gcd_anc[sel[i]]=unlist(gcd(tab_coo_f[anc,1], tab_coo_f[anc,2], tab[sel[i],1], tab[sel[i],2]))
}

##East
sel=rownames(tab[which(tab$Cluster=="E"),])

anc="56"
for(i in 1:length(sel)) {
  gcd_anc[sel[i]]=unlist(gcd(tab_coo_f[anc,1], tab_coo_f[anc,2], tab[sel[i],1], tab[sel[i],2]))
}

gcd_anc=gcd_anc[pops_sel]


gcd_mix_anc=gcd_tree_anc
gcd_mix_anc[which(is.na(gcd_mix_anc))]=gcd_anc[which(is.na(gcd_mix_anc))]

#####
#create table

tab_sel=as.data.frame(cbind(
  
  gcd_anc=gcd_anc,
  gcd_mix_anc=gcd_mix_anc,
  gcd_centroid=gcd_centroid,
  gcd_tree_anc=gcd_tree_anc,
  
  tot_theta_0.5=tab_tot$Wat_w_0.5,
  tot_pi_0.5=tab_tot$Pi_w_0.5,
  tot_Tajima_D_0.5=tab_tot$Tajima_D_w_0.5,
  tot_FayWu_H_0.5=tab_tot$FayWu_H_w_0.5,
  
  inter_theta_0.5=tab_inter$Wat_w_0.5,
  inter_pi_0.5=tab_inter$Pi_w_0.5,
  inter_Tajima_D_0.5=tab_inter$Tajima_D_w_0.5,
  inter_FayWu_H_0.5=tab_inter$FayWu_H_w_0.5,
  
  intro_theta_0.5=tab_intro$Wat_w_0.5,
  intro_pi_0.5=tab_intro$Pi_w_0.5,
  intro_Tajima_D_0.5=tab_intro$Tajima_D_w_0.5,
  intro_FayWu_H_0.5=tab_intro$FayWu_H_w_0.5,
  
  cds_theta_0.5=tab_cds$Wat_w_0.5,
  cds_pi_0.5=tab_cds$Pi_w_0.5,
  cds_Tajima_D_0.5=tab_cds$Tajima_D_w_0.5,
  cds_FayWu_H_0.5=tab_cds$FayWu_H_w_0.5,
  
  Dos_0.5=tab_cds$Dos_w_0.5,
  
  lat=tab52$y,
  lon=tab52$x,
  Cluster=as.character(tab52$GeneticCluster),
  Cluster3=Cluster3,

  Mat=as.character(tab52$MatingSystem),
  Sub=as.character(tab52$Substrate),
  Disturbance=as.character(tab52$Disturbance..1.3.),
  Canopy=tab52$CanopyCover....,
  Admix=Admix,
  AdmixWE=AdmixWE,
  Ice=Ice,
  
  A=tab_ms$A,
  He=tab_ms$He,
  Cens=tab52$PopSizeBoltedPlants,
  Cens_log=log10(tab52$PopSizeBoltedPlants),
  Dens=tab52$DensityBolted..m.2.,
  Dens_log=log10(tab52$DensityBolted..m.2.)))


rownames(tab_sel)=rownames(tab_cds)
write.table(tab_sel,file=nome_out,quote=FALSE,sep="\t")


