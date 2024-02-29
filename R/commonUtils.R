library('dplyr')
library('stringdist')

# 判断df中某列是否有缺失值,是的话就返回true
isNaContain = function(df, feature){
  if( TRUE %in% is.na(df[,feature])){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
# 创建文件夹
mkDir = function(pathStr){
  if (!dir.exists(pathStr)){
    dir.create(pathStr, recursive = T)
  }
}

# 从一个特征list中找到和fstr最为接近的一个值并返回, 相似度度量方法是Levenshtein距离
getSimiFeatureName = function(fStr, featureList){
  simi = unlist(lapply(featureList, function(x) stringdist(a = fStr, b =x, method = 'lv')))
  return(featureList[which.min(simi)])
}

#保存plot的代码段, 传进来的是一个plot实例
savePlotInstance = function(pngName, width, height, plotInstance){
  png(
    filename = pngName, # 文件名称
    width = width,           # 宽
    height = height,          # 高
    units = "px",          # 单位
    bg = "white",          # 背景颜色
    res = 300)              # 分辨率
  plotInstance
  dev.off()
}
#将df中的某个连续变量按照某（几）个阈值离散化
discreateLvl = function(df, mdName, lvlThres){
  lvlNums = length(lvlThres) + 1
  newName = paste(mdName, '_lvl', lvlNums, sep = '')
  df[, newName] = lvlNums
  for (i in lvlThres){
    idx = which(df[, mdName] < i)
    df[idx, newName] = df[idx, newName] - 1
  }
  df[, newName] = as.factor(df[, newName])
  return(df)
}
singleP4Show = function(p, roundNum = 3){
  if (p < 0.001){
    rst = '<0.001'
  }else{
    pRound = round(p, roundNum)
    rst = str_pad(ifelse (p < 0.05 & pRound == 0.05
                          , str_sub(p, 1, roundNum + 2)
                          , pRound)
                  ,roundNum + 2, side = 'right', pad = '0')               
  }
  return(rst)
}
# 将data中的数据, 从oriMap的枚举中映射成新的map
reMap = function(vec, oriMap, newMap){
  if (length(oriMap) != length(newMap)){
    stop('新旧枚举映射需要长度一致')
  }
  vecNew = vec
  for (i in seq_len(length(oriMap))){
    idx = which(vec == oriMap[i])
    vecNew[idx] = newMap[i]
  }
  return(vecNew)
}

#小数位数不够补0
pointsZeroFill = function(num, pointsWidth = 1){
  roundNum = round(num, pointsWidth)
  if (is.na(roundNum)){
    return(NA)
  }
  # 1. 先检测四舍五入了之后是否是整数
  if(!str_detect(roundNum, '\\.')){
    return(paste(roundNum, rep(0, pointsWidth), sep = '.'))
  }else{#确认是小数就将小数位数给补齐
    sp = str_split(roundNum, '\\.')[[1]]
    return(paste(sp[1], str_pad(sp[2], pointsWidth, side = 'right', pad = '0'), sep = '.'))
    }
}
