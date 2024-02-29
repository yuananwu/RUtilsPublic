library("raters")

isDfContainNa = function(df){
  uniqueEle = unique(unlist(df))
  if( TRUE %in% is.na(uniqueEle)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# 将多位raters的评分数据, 转化为评分的频数
df2Freq = function(df, na.rm = F){
  uniqueFactor = sort(unique(unlist(df)))
  fNum = length(uniqueFactor)
  if (fNum == 1){
    stop('unique factor number equal 1')
  }
  if (isDfContainNa(df)){
    if(na.rm){
      df = na.omit(df)
      cat(names(df), '存在空值, 计算过程会自动将空值剔除掉, 剩余样本量: ', dim(df)[1], '\n')
    }else{
      stop('missing values exist')
    }
  }
  rst = NULL
  for (i in 1:nrow(df)){
    rstTmp = rep(0,fNum)
    for (j in 1:ncol(df)){
      for (idx in 1:fNum){
        if(df[i,j] == uniqueFactor[idx]){
          rstTmp[idx] = rstTmp[idx] + 1
        }
      }
    }
    rst = rbind(rst, rstTmp)
  }
  return(as.matrix.data.frame(rst))
}

# 当评分的level只有两个时, 调用该方法
getKappaFor2Levels = function(df){
  # Monte Carlo检验方法用来计算Fleiss Kappa的95CI,适用于小数据集
  #检验方法选择参考: https://cran.r-project.org/web/packages/raters/raters.pdf
  kp = concordance(df,test="MC",B=100,alpha = 0.05)
  value = kp$Fleiss
  rst = paste(round(value[1], 3), ' (',round(value[2], 3), '-', round(value[3], 3), ')',   sep = '')
  return(rst)
}

# 当评分的level大于两个时, 调用该方法
getKappaForMultiLevels = function(df){
  set.seed(42)
  kp = wlin.conc(df,test="MC",B=30)
  rst = paste(round(kp[1], 3), ' (',round(kp[3], 3), '-', round(kp[4], 3), ')',   sep = '')
  return(rst)
}
