#This code were used for normalise of raf22raf36 phosphoproteome
#cewated on : 2020/Sep/30
#Author : Sotaro Katagiri

library(minpack.lm) #非線形最小二乗法を行うのに必要なパッケージ
library(snowfall) # 並列処理に必要なパッケージ

Data = read.csv('norm.csv', header = T, row.name = 1) #データフレーム型で読み込む
Data <- matrix(as.matrix(Data), nrow(Data), ncol(Data)) #行列型に変換
samplenum = 24
parStart = numeric(samplenum * 3) # 係数の初期値。72はファイル数
parStart = parStart + c(1, 1, 1)  #初期値は前回行った結果をつかうとうまく計算できる。


peparea = function(X, A) {               #ペプチドごとの面積 Aが係数ベクトル,Xが単一ペプチドの実測値を収めたベクトル
  B = A * X
  I = numeric(samplenum)  #サンプルごとの定量値をしまう空箱 戻り値にする。
  for (i in 1:samplenum) {
    I[i] = B[3 * i - 2] + B[3 * i - 1] + B[3 * i] #３ファイルを足す
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
  return(sqrt(res)) # nls.lm関数で初手二乗されるため、ここでは平方根を取っておく。
}


sfInit(parallel = TRUE, cpus = 12) #パソコンのスレッド数
sfExportAll()

resid = function(par, data) {
  res = sfApply(data, 1, peparea, par)  #上の式をすべての行（ペプチド）に対し計算し、ベクトルとして受け取る。
  return(res)
}

LM = nls.lm(par = parStart, lower = numeric(samplenum * 3) + 0.1, fn = resid, data = Data, control=nls.lm.control(maxiter=128,nprint=1))

sfRemoveAll()
sfStop()

LM
Data = read.csv('norm.csv',header = T)
for(row in 1:dim(Data)[1]) { #pepmodseqの列を除くため[1,2]スタート
  for(col in 2:dim(Data)[2]) {
    Data[row,col] <- Data[row,col] * LM$par[col - 1]
  }
}
write.csv(Data,
  "normarised.csv",
  row.names = FALSE
)