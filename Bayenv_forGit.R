###########Script to analyse output from bayenv env assoc runs######
####author: Mitra Menon
####date: 12-July-2018####
##modified: 16-Jan-2020##
################NOTE:#########################
#The script is made to fit the folder and file name structure used for my outputs from bayenv. Please modify the paths as seen fit. It also assumes that the bayenv output for climate and soil are stored in two seperate directories. If you combined them all into one file, you will only need to run one of the script from PART-1 and specify the directory name in PART-2
###
###################################################################

library(plyr)
library(data.table)

#there should be three replicate runs for each env varible, obtained from 3 runs in bayenv
#These files are not provided.
climFiles<-list.files("clim/",pattern = "Bay",full.names = T)
soilFiles<-list.files("Soil/",pattern = "Bay",full.names = T)

##PART1#################################################################################
################# generate outliers shared between rho and BF for each predictor variable and write them out (This needs to be done for each replicate file)
#######################################################################################

######FOR SOIL VARIABLES################
for (s in 1:length(soilFiles)){
  
  Pop98<-read.table("98PopIDs",sep="\t")
  run<-paste0("run",s)
  
  cat("Reading file",soilFiles[s],"\n")
    BFsoil <- fread(soilFiles[s],sep="\t",data.table=F)
    top<-round(nrow(BFsoil)*0.01) #the cutoff used to declare SNPs as significant
    
    soil<-read.table("SoilGrid1km.txt",header=T,sep="\t")
    soil98<-soil[soil$Site%in%Pop98$V1, ]
    
    colnames(BFsoil)[1]<-"ID"
    colnames(BFsoil)[2:31]<-rep(colnames(soil98[c(4:13)]),each=3) #triplicates due to three estimates from bayenv
    rownames(BFsoil)<-BFsoil$ID
    BFsoil<-BFsoil[ ,-c(1)]
    
    #Generate summary file for each variable using the cutoff
    output<-as.data.frame(matrix(nrow=1,ncol=5))
    colnames(output)<-c("bf_median","bf_min","bf_max","rho_avg","overlap")
    
    a<-1 #start at col1
    b<-3 #3 due to the three estimates for each loci per variable given by bayenv
    
    while(b<=ncol(BFsoil)){
      
      cat("working on",colnames(BFsoil[a]),"\n")
      
      #first look at bf
      var<-BFsoil[ ,a:b]
      var_bf_soil<-var[order(var[ ,1],decreasing=T), ]
      var_bf_soil<-var_bf_soil[1:top,]
      bf_min<-min(var_bf_soil[ ,1])
      bf_max<-max(var_bf_soil[ ,1])
      OUT<-c(median(var_bf_soil[ ,1]),bf_min,bf_max)
      
      #now look at rho
      var[ ,2]<-abs(var[ ,2])
      var_rho<-var[order(var[ ,2],decreasing=T), ]
      var_rho<-var_rho[1:top, ]
      OUT<-c(OUT,mean(var_rho[ ,2]))
      
      #shared outliers between rho and bf
      common_var<-var_bf_soil[rownames(var_bf_soil)%in%rownames(var_rho), ]
      common_var$ID<-rownames(common_var)
      
      OUT<-c(OUT,nrow(common_var))
      ID<-colnames(var_bf_soil[a])
      
      outDIR<-paste0("Soil/",run)
      dir.create(outDIR)
      write.table(common_var,file=paste0(outDIR,"/overlaptop_",ID,".txt"),sep="\t",row.names = F,quote=F)
      
      output<-rbind(output,OUT)
      
      a<-b+1
      b<-b+3 #moves to the next variable (skip col 2 and col 3)
    }
    
    output<-output[-1,]
    output$predictor<-colnames(BFsoil)[seq(1,ncol(BFsoil),3)]
    summLOC<-paste0(outDIR,"/SummaryTop1_SOIL_scaled.txt")
    write.table(output,file=summLOC,sep="\t",row.names = F,quote=F)
}
  
############NOW FOR CLIMATE DATASET############################

for (c in 1:length(climFiles)){
  
  Pop98<-read.table("98PopIDs",sep="\t")
  
  run<-paste0("run",c)
  
  cat("Reading file",climFiles[c],"\n")
  
  bf_clim<-fread(climFiles[c],sep="\t",data.table=F)
  str(bf_clim)
  clim<-read.table("periphery_climate1981-2010.txt",header=T,sep="\t")
  clim98<-clim[clim$Site%in%Pop98$V1, ]
  
  colnames(bf_clim)[1]<-"ID"
  colnames(bf_clim)[2:244]<-rep(colnames(clim98[c(2:82)]),each=3) 
  rownames(bf_clim)<-bf_clim$ID
  bf_clim<-bf_clim[ ,-c(1)]
  
  top<-round(nrow(bf_clim)*0.01) #cutoff
  
  output<-as.data.frame(matrix(nrow=1,ncol=5))
  colnames(output)<-c("bf_median","bf_min","bf_max","rho_avg","overlap")
  
  a<-1 #start at col1
  b<-3 #3 due to the three estimates for each loci per variable given by bayenv
  
  while(b<=ncol(bf_clim)){
    
    cat("working on",colnames(bf_clim[a]),"\n")
    
    var<-bf_clim[ ,a:b]
    var_bf_clim<-var[order(var[ ,1],decreasing=T), ]
    var_bf_clim<-var_bf_clim[1:top,]
    bf_min<-min(var_bf_clim[ ,1])
    bf_max<-max(var_bf_clim[ ,1])
    OUT<-c(median(var_bf_clim[ ,1]),bf_min,bf_max)
    
    var[ ,2]<-abs(var[ ,2])
    var_rho<-var[order(var[ ,2],decreasing=T), ]
    var_rho<-var_rho[1:top, ]
    OUT<-c(OUT,mean(var_rho[ ,2]))
    
    common_var<-var_bf_clim[rownames(var_bf_clim)%in%rownames(var_rho), ]
    common_var$ID<-rownames(common_var)
    
    OUT<-c(OUT,nrow(common_var))
    ID<-colnames(bf_clim[a])
    
    outDIR<-paste0("clim/",run)
    dir.create(outDIR)

    write.table(common_var,file=paste0(outDIR,"/overlaptop_",ID,".txt"),sep="\t",row.names = F,quote=F)
    
    output<-rbind(output,OUT)
    
    a<-b+1
    b<-b+3 
  }
  
  output<-output[-1,]
  output$predictor<-colnames(bf_clim)[seq(1,ncol(bf_clim),3)]
  summLOC<-paste0(outDIR,"/SummaryTop1_climNA_scaled.txt")
  write.table(output,file=summLOC,sep="\t",row.names = F,quote=F)
  
}
  

############################################################################################
##PART2###############
#####COMBINING OUTLIERS FROM ALL THREE RUNS TO GET ONLY SHARED OUTLIERS AS THE MOST STRINGENT SET
########################################################################################

step2JOIN<-function(type){
  
  ######type is either Soil or climate### 
  #(needs to exactly match the folder names that holds output from part1 for soil and clim)
  #######################################
  
  DIR<-paste0("./",type)
  varNames<-list.files(paste0(DIR,"/run1"),pattern="overlap") #JUST OBTAINING ALL FILE NAMES
  
  commonOverlap<-vector("list",length(varNames))
  
  for (i in 1:length(varNames)){
    
    print(varNames[i])
    
    
    df1<-read.table(paste0(DIR,"/run1/",varNames[i]),sep="\t",header=T)
    df2<-read.table(paste0(DIR,"/run2/",varNames[i]),sep="\t",header=T)
    df3<-read.table(paste0(DIR,"/run3/",varNames[i]),sep="\t",header=T)
    
    
    top1Stringet<-join_all(list(df1,df2,df3),by="ID",match="all",type = "inner")
    convDIR<-paste0(DIR,"/convergence")
    dir.create(convDIR)
    write.table(top1Stringet,file=paste0(convDIR,"/Overlap1_3Chains",varNames[i]),sep="\t",
                row.names = F,quote=F)
    
  }
  
}

  



