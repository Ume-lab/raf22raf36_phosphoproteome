#This code were used for analysing of raf22raf36 phosphoproteome data
#cewated on : 2020/Sep/30
#Author : Sotaro Katagiri

setwd('  ') 
csv = read.csv("forRforR.csv", header = T)
peps = levels(csv$peptide) 
pepn = length(peps) 

countUP <- function(foldchange, p) {   
  if (p >= 0.05) return(0)
  if (foldchange > 2) return(1)
  return(0)
}

Result = data.frame(                
  pepseq = numeric(pepn),
  WT0 = numeric(pepn),
  WT15 = numeric(pepn),
  WT30 = numeric(pepn),
  WT90 = numeric(pepn),
  DKO0 = numeric(pepn),
  DKO15 = numeric(pepn),
  DKO30 = numeric(pepn),
  DKO90 = numeric(pepn),
  fc_WT15_WT0 = numeric(pepn),
  fc_WT30_WT0 = numeric(pepn),
  fc_WT90_WT0 = numeric(pepn),
  fc_DKO0_WT0 = numeric(pepn),
  fc_DKO15_WT15 = numeric(pepn),
  fc_DKO30_WT30 = numeric(pepn),
  fc_DKO90_WT90 = numeric(pepn),
  p_WT15_WT0 = numeric(pepn),
  p_WT30_WT0 = numeric(pepn),
  p_WT90_WT0 = numeric(pepn),
  p_DKO0_WT0 = numeric(pepn),
  p_DKO15_WT15 = numeric(pepn),
  p_DKO30_WT30 = numeric(pepn),
  p_DKO90_WT90 = numeric(pepn),
  ABA_up = numeric(pepn),
  ABA_down = numeric(pepn),
  DKO_up = numeric(pepn),
  DKO_down = numeric(pepn)
)




for(i in 1:pepn){      
  p = peps[i]
  a = csv[csv$peptide == p,]  
  Result[i,"pepseq"] = p
  
  Result[i,"WT0"] = mean(a[a$sample == "WT0","area"])    
  Result[i,"WT15"] = mean(a[a$sample == "WT15","area"])
  Result[i,"WT30"] = mean(a[a$sample == "WT30","area"])
  Result[i,"WT90"] = mean(a[a$sample == "WT90","area"])
  Result[i,"DKO0"] = mean(a[a$sample == "DKO0","area"])
  Result[i,"DKO15"] = mean(a[a$sample == "DKO15","area"])
  Result[i,"DKO30"] = mean(a[a$sample == "DKO30","area"])
  Result[i,"DKO90"] = mean(a[a$sample == "DKO90","area"])
  
  Result[i,"fc_WT15_WT0"] = Result[i,"WT15"] / Result[i,"WT0"]   
  Result[i,"fc_WT30_WT0"] = Result[i,"WT30"] / Result[i,"WT0"]
  Result[i,"fc_WT90_WT0"] = Result[i,"WT90"] / Result[i,"WT0"]
  Result[i,"fc_DKO0_WT0"] = Result[i,"DKO0"] / Result[i,"WT0"]
  Result[i,"fc_DKO15_WT15"] = Result[i,"DKO15"] / Result[i,"WT15"]
  Result[i,"fc_DKO30_WT30"] = Result[i,"DKO30"] / Result[i,"WT30"]
  Result[i,"fc_DKO90_WT90"] = Result[i,"DKO90"] / Result[i,"WT90"]
  
  Result[i,"p_WT15_WT0"] = t.test(a[a$sample == "WT15","area"], a[a$sample == "WT0","area"],var.equal=T)$p.value      
  Result[i,"p_WT30_WT0"] = t.test(a[a$sample == "WT30","area"], a[a$sample == "WT0","area"],var.equal=T)$p.value
  Result[i,"p_WT90_WT0"] = t.test(a[a$sample == "WT90","area"], a[a$sample == "WT0","area"],var.equal=T)$p.value
  Result[i,"p_DKO0_WT0"]  = t.test(a[a$sample == "DKO0","area"], a[a$sample == "WT0","area"],var.equal=T)$p.value
  Result[i,"p_DKO15_WT15"] = t.test(a[a$sample == "DKO15","area"], a[a$sample == "WT15","area"],var.equal=T)$p.value
  Result[i,"p_DKO30_WT30"] = t.test(a[a$sample == "DKO30","area"], a[a$sample == "WT30","area"],var.equal=T)$p.value
  Result[i,"p_DKO90_WT90"] = t.test(a[a$sample == "DKO90","area"], a[a$sample == "WT90","area"],var.equal=T)$p.value
  
  
  
  
  Result[i,"ABA_up"] = (countUP(Result[i,"fc_WT15_WT0"], Result[i,"p_WT15_WT0"]) 
  + countUP(Result[i,"fc_WT30_WT0"],Result[i,"p_WT30_WT0"]) 
  + countUP(Result[i,"fc_WT90_WT0"],Result[i,"p_WT90_WT0"]))
  
  Result[i,"ABA_down"] = (countUP(1/ Result[i,"fc_WT15_WT0"], Result[i,"p_WT15_WT0"])
   + countUP(1/ Result[i,"fc_WT30_WT0"],Result[i,"p_WT30_WT0"]) 
   + countUP(1/ Result[i,"fc_WT90_WT0"],Result[i,"p_WT90_WT0"]))
  
  Result[i,"DKO_up"] = (countUP(Result[i,"fc_DKO0_WT0"], Result[i,"p_DKO0_WT0"]) 
  + countUP(Result[i,"fc_DKO15_WT15"],Result[i,"p_DKO15_WT15"]) 
  + countUP(Result[i,"fc_DKO30_WT30"],Result[i,"p_DKO30_WT30"]) 
  + countUP(Result[i,"fc_DKO90_WT90"],Result[i,"p_DKO90_WT90"]))
  
  Result[i,"DKO_down"] = (countUP(1/ Result[i,"fc_DKO0_WT0"], Result[i,"p_DKO0_WT0"])
  + countUP(1/ Result[i,"fc_DKO15_WT15"],Result[i,"p_DKO15_WT15"])
  + countUP(1/ Result[i,"fc_DKO30_WT30"],Result[i,"p_DKO30_WT30"])
  + countUP(1/ Result[i,"fc_DKO90_WT90"],Result[i,"p_DKO90_WT90"]))
}

write.csv(Result, "rsf2236kaiseki001.csv",row.names=FALSE)


#version 1 written by sotaro katagiri 2020/apr/23
#R version 3.6.3
