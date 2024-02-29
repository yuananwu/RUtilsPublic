library("raters")
library('irr')

# 判断df中是否有缺失值,是的话就返回true
isDfContainNa = function(df){
  uniqueEle = unique(unlist(df))
  if( TRUE %in% is.na(uniqueEle)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# 将多rater对n类的评分数据, 在每个样本层面,转化为n类的频数
oriDataToLevelsFreq = function(df){
  raterNum = ncol(df)
  sampleNum = nrow(df)
  level =  sort(unique(unlist(df)))
  levelNum = length(level)
  if (isDfContainNa(df)){
    stop("数据组中存在缺失值, 请检查数据")
  }
  if (levelNum == 1){
    stop("数据组的枚举类型只有一种, 请检查数据")
  }
  freqVec = vector(mode="numeric",length=levelNum*sampleNum)
  i = 1
  while(i <= sampleNum){
    j = 1
    while(j <= raterNum){
      for (eleIndex in 1: levelNum){
        if( df[i,j] == level[eleIndex] ){
          freqIndex = levelNum * i - (eleIndex-1)
          freqVec[freqIndex] = freqVec[freqIndex] + 1
        }
      }
      j = j+1
    }
    i = i+1
  }
  result = matrix(freqVec,ncol = levelNum, byrow = T)
  return(result)
}

# 计算flesis kappa值
calculateFlesisKappaFor2Levels = function(dfFreq){
  # Monte Carlo检验方法用来计算Fleiss Kappa的95CI,适用于小数据集
  #检验方法选择参考: https://cran.r-project.org/web/packages/raters/raters.pdf
  f_k_rsl = concordance(dfFreq,test="MC",B=100,alpha = 0.05)
  kappaValue = round(f_k_rsl$Fleiss[1],3)
  kappaValue95L = round(f_k_rsl$Fleiss[2],3)
  kappaValue95U = round(f_k_rsl$Fleiss[3],3)
  result = paste(kappaValue, ' (', kappaValue95L, '-', kappaValue95U, ')', sep = '')
  cat('fleiss kappa:', result)
  return(c('fleiss kappa', result))
}

#这是多raters&多levels的计算方法
calculateFlesisKappaForMultiLevels = function(df){
  set.seed(42)
  kappa = wlin.conc(df,test="MC",B=25)
  kappaValue = round(kappa[1],3)
  kappaValue95L = round(kappa[3],3)
  kappaValue95U = round(kappa[4],3)
  result = paste(kappaValue, ' (', kappaValue95L, '-', kappaValue95U, ')', sep = '')
  cat('weighted kappa:', result)
  return(c('weighted kappa', result))
}
#连续变量计算icc
calculateIcc = function(df){
  for (i in names(df)){
    df[,i] = as.numeric(df[,i])
  }
  icc = icc(df, model="twoway", type="agreement")
  iccValue = round(icc$value,3)
  iccValue95L = round(icc$lbound,3)
  iccValue95U = round(icc$ubound,3)
  result = paste(iccValue, ' (', iccValue95L, '-', iccValue95U, ')', sep = '')
  cat('icc:', result)
  return(c('icc', result))
}
