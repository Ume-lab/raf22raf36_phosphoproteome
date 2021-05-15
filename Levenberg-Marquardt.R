#This code were used for normalise of raf22raf36 phosphoproteome
#cewated on : 2020/Sep/30
#Author : Sotaro Katagiri

library(minpack.lm) 
library(snowfall) 

Data = read.csv('norm.csv', header = T, row.name = 1) 
Data <- matrix(as.matrix(Data), nrow(Data), ncol(Data)) 
samplenum = 24
parStart = numeric(samplenum * 3) 
parStart = parStart + c(1, 1, 1) 


peparea = function(X, A) {               
  B = A * X
  I = numeric(samplenum)  
  for (i in 1:samplenum) {
    I[i] = B[3 * i - 2] + B[3 * i - 1] + B[3 * i] 
  }
  res = 0
  for(i in 1:samplenum) {
    for(j in i+1:samplenum) {
      a = 0
      b = 0
      if(I[i] > 1) a = log(I[i])
      if((j <= samplenum) && (I[j] > 1)) b = log(I[j])
      if((a > 0) && (b > 0)) res = res + (a - b) ^ 2 
    }
  }
  return(sqrt(res)) 
}


sfInit(parallel = TRUE, cpus = 12) 
sfExportAll()

resid = function(par, data) {
  res = sfApply(data, 1, peparea, par) 
  return(res)
}

LM = nls.lm(par = parStart, lower = numeric(samplenum * 3) + 0.1, fn = resid, data = Data, control=nls.lm.control(maxiter=128,nprint=1))

sfRemoveAll()
sfStop()

LM
Data = read.csv('norm.csv',header = T)
for(row in 1:dim(Data)[1]) { 
  for(col in 2:dim(Data)[2]) {
    Data[row,col] <- Data[row,col] * LM$par[col - 1]
  }
}
write.csv(Data,
  "normarised.csv",
  row.names = FALSE
)
