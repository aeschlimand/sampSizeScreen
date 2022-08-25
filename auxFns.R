 # This file contains helper f'ns for the Shiny app 'SSMDSampleSize'
 getSSMDSingle=function(x, muNow, nMi, dfPerFOVb, XbarDMSOb, 
   varDMSO, sdModelb, sdNow, nDMSOb){
   sampNow=(sdNow/sdModelb)*rt(nMi, dfPerFOVb)+muNow
   Xbar=mean(sampNow)
   varMi=var(sampNow)
   df=(varMi/nMi+varDMSO/nDMSOb)^2/((((varMi/nMi)^2)/(nMi-1))+
     (((varDMSO/nDMSOb)^2)/(nDMSOb-1)))
   rawDiff=Xbar-XbarDMSOb
   betaHat=rawDiff/(sqrt(varMi + varDMSO))
   t=rawDiff/(sqrt(varMi/nMi + varDMSO/nDMSOb))
   return(c(betaHat, Xbar, t, df))  
 }

 getLevels=function(x) switch(x, "01 ExtWeak", "02 VeryWeak", "03 Weak", "04 FairWeak",
   "05 FairMod", "06 Mod", "07 FairStrong", "08 Strong", "09 VeryStrong",
   "10 ExtStrong")
 # getLevels2 assigns midpoints of the intervals
 getLevels2=function(x){
   if(x==0.375){
     y="VeryWeak=0.375"
   }else if(x==0.625){
     y="Weak=0.625"
   }else if(x==0.875){
     y="FairWeak=0.875"
   }else if(x==1.14){
     y="FairMod=1.14"
   }else if(x==1.46){
     y="Mod=1.46"
   }else if(x==1.82){
     y="FairStrong=1.82"
   }else if(x==2.5){
     y="Strong=2.5"
   }else if(x==4){
     y="VeryStrong=4"
   }else if(x==7.5){
     y="ExtStrong=7.5"
   }else y="StrangeSSMD"
   return(y)
 }
 getGraphCentered=function(SSMD.Now, sd.Now, metricNow, dfPerFOV, mu.DMSO, sdDMSO, colsNowA, compClass){
   sdModel=sqrt(dfPerFOV /(dfPerFOV-2)) # Variance of a t-distribution is df/(df-2)
   # Read in values for generating SDs and dist'n of DMSO ----
   cat('hereAA\n')
   # Values are read in, start crunching ----
   # Get mu (mean of the FOV measurements of the current metric) for this compound. ----
   mu.Now=SSMD.Now*sqrt(sdDMSO^2+sd.Now^2)+mu.DMSO
   seqq=seq(0.0001, 0.9999, length=1e4)
   tsNow=qt(seqq, dfPerFOV)
   pdfNow=dt(tsNow, dfPerFOV)
   FOVsNow=mu.Now+(sd.Now/sdModel)*tsNow
   DMSOFOVsNow=mu.DMSO+(sdDMSO/sdModel)*tsNow
   bigDat=data.frame(explanVar=c(FOVsNow, DMSOFOVsNow), compClass=c(rep(compClass, 1e4), rep('DMSO', 1e4)))
   if(metricNow=="MPDC"){
     xlabNow=expression(paste("Centered MPDC (", mu, m^2, "/sec )"))
   }else{
     xlabNow=expression(paste("Centered Q3JL (", mu, m, ")"))
   }
   p <- bigDat %>%
     ggplot( aes(x=explanVar, fill=compClass, y=..density..)) +
     geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins=30) +
     scale_fill_manual(values=colsNowA[c(3,8)]) +
     # theme_ipsum() +
     labs(fill="")+ xlab(xlabNow)+ylab("Density")
   p
   # par(mar=c(5, 5, 1, 1))
   # if(metricNow=="MPDC"){
   #   (FOVsNow, pdfNow, pch=16, main='', col=colsNowA[3],
   #        xlab=expression(paste("Centered MPDC (", mu, m^2, "/sec )")), 
   #        ylab='PDF', cex.lab=1.2, cex.axis=1.2)
   # }else if(metricNow=="Q3JL"){
   #   plot(FOVsNow, pdfNow, pch=16, main='', col=colsNowA[3],
   #        xlab=expression(paste("Centered Q3JL (", mu, m, ")")), 
   #        ylab='PDF', cex.lab=1.2, cex.axis=1.2)
   # }
   # # abline(v=mu.DMSO, col=colsNowA[8], lwd=4)
   # text(mu.DMSO, 0.2*min(pdfNow)+0.8*max(pdfNow), "DMSO", cex=1.5, col=colsNowA[8], pos=2)
   # points(FOVsNow, pdfNow, pch=16, col=colsNowA[3])
   # points(DMSOFOVsNow, pdfNow, pch=16, col=colsNowA[3])
 }
 
 # getLevels(10)
 # getScreenResults=function(x){
 #   # Available variables are pi0, nDMSO, SSMD.t1, nCompounds, and nFOVs.t1
 #   # How many active compounds do we have?
 #   samp.DMSO=rt(nDMSO, dfPerFOV)
 #   meanDMSO=mean(samp.DMSO)
 #   varDMSO=var(samp.DMSO)
 #   vecc=1:nActives
 #   SSMDActives=sort(sapply(vecc, getSSMDSingle, mu2Use=mu.t1))
 #   SSMDNonActives=sort(sapply(vecc, getSSMDSingle, mu2Use=0))
 #   loLim95Actives=SSMDActives[fivePercentActives]
 #   hiLim95Actives=SSMDActives[ninetyFivePercentActives]
 #   loLim95NonActives=SSMDNonActives[fivePercentNonActives]
 #   hiLim95NonActives=SSMDNonActives[ninetyFivePercentNonActives]
 #   medianActives=median(SSMDActives)
 #   medianNonActives=median(SSMDNonActives)
 #   return(c(loLim95Actives, medianActives, hiLim95Actives, loLim95NonActives,
 #     medianNonActives, hiLim95NonActives))
 # }
 
 
 
 
