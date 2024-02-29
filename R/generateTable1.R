
if (!require('car')) {
  install.packages('car')
  library('car')
} else {
  library('car')
}
if (!require('stringr')) {
  install.packages('stringr')
  library('stringr')
} else {
  library('stringr')
}

if (!require('tableone')) {
  install.packages('tableone')
  library('tableone')
} else {
  library('tableone')
}

# 连续变量
normStatusTest = function(strataVar, strataLvl, contVar){
  tryCatch({
    # 先对每一层的数据进行Shapiro—Wilk test， 看是否满足检验（p > 0.05）
    SWFlag = 0
    for (lvl in strataLvl){
      p = shapiro.test(contVar[which(strataVar == lvl)])$p.value
      if (p > 0.05){SWFlag = SWFlag + 1}
    }
    if (SWFlag == length(strataLvl)){ # 代表都通过了SW test
      #接下来进行方差齐性检验, P>0.05为方差齐
      p = leveneTest(contVar~factor(strataVar),data=as.data.frame(cbind(contVar, strataVar)))$`Pr(>F)`
      if (p[1] > 0.05){
        return('通过了方差齐性和正态性检测')
      }else{
        return('未通过方差齐性检测')
      }
    }else{
      return('未通过正态性检测')
    }
  },error = function(err) {
    return('出现错误')
  }
  )
}
compareGroupCont = function(strataVar, contVar){
  if (length(strataVar) != length(contVar)){
    stop('输入的分层变量和连续特征长度须一致')
  }
  strataVar = as.factor(strataVar)
  strataLvl = levels(strataVar)
  contVar   = as.numeric(contVar)
  
  rst = list()
  if(length(strataLvl) < 2){
    stop("分组枚举类型小于2个, 请检查输入数据")
  }else if (length(strataLvl) == 2){
    if( normStatusTest(strataVar, strataLvl, contVar) == '通过了方差齐性和正态性检测'){
      # t-检验
      p = t.test(contVar~factor(strataVar),data=as.data.frame(cbind(contVar, strataVar)),var.equal=T)$p.value
      rst[['检验方法']] = '正态_两组_t检验'
      rst[['p值']] = p
      rst[['p值解释']] = ifelse(p < 0.05, '组间存在显著差异', '组间不存在显著差异')
      
    }else{
      # Mann-Whitney U-Test (Wilcoxon秩和检验)
      p = wilcox.test(contVar~factor(strataVar),data=as.data.frame(cbind(contVar, strataVar)))$p.value
      rst[['检验方法']] = '非正态_两组_Mann-Whitney U检验'
      rst[['p值']] = p
      rst[['p值解释']] = ifelse(p < 0.05, '组间存在显著差异', '组间不存在显著差异')
      
    }
  }else{
    if(length(strataLvl) > 5){warning("分组枚举类型大于5个， 请注意")}
    if( normStatusTest(strataVar, strataLvl, contVar) == '通过了方差齐性和正态性检测'){
      # ANOVA（方差分析）
      aovRst = aov(contVar~factor(strataVar),data=as.data.frame(cbind(contVar, strataVar)),var.equal=T)
      p = summary(aovRst)[[1]]$`Pr(>F)`[1]
      rst[['检验方法']] = '正态_多组_ANOVA'
      rst[['p值']] = p
      rst[['p值解释']] = ifelse(p > 0.05, '各组总体均数全部相等', '各组总体均数不全相等')
    }else{
      # Kurskal-Wallis检验
      p = kruskal.test(contVar~factor(strataVar),data=as.data.frame(cbind(contVar, strataVar)))$p.value
      rst[['检验方法']] = '非正态_多组_Kurskal-Wallis检验'
      rst[['p值']] = p
      rst[['p值解释']] = ifelse(p > 0.05, '各组总体均数全部相等', '各组总体均数不全相等')
    }
  }
  return(rst)
}

#离散变量

chisqRoute = function(tb){
  #先计算理论频数
  n  = sum(tb)
  sr = rowSums(tb)
  sc = colSums(tb)
  e  = outer(sr, sc)/n # 将原始table转换为理论频数的table
  
  if(any(dim(tb) < 2)){
    return('不满足最低2x2的要求')
  }else if (all(dim(tb) == 2)){
    if(n < 40 | any(e < 1)){
      return('fisher')
    }else if(any(e < 5)){
      return('校正卡方')
    }else{
      return('卡方')
    }
  }else{
    # 获取RC列联表中的小于5的格子占比
    ratio = length(which(e < 5))/length(e)
    if (any(e < 1)){
      return('fisher')
    }else if (ratio > 0.2){
      return('fisher')
    }else{
      return('卡方')
    }
  }
}

compareGroupCat = function(strataVar, catVar){
  tryCatch({
    if (length(strataVar) != length(catVar)){
      stop('输入的分层变量和离散特征长度须一致')
    }
    strataVar = as.factor(strataVar)
    strataLvl = levels(strataVar)
    catVar   = as.factor(catVar)
    catLvl = levels(catVar)
    tb = table(strataVar, catVar)
    
    flag = chisqRoute(tb)
    rst = list()
    if(flag == 'fisher'){
      p = fisher.test(tb)$p.value
      rst[['检验方法']] = 'fisher精确检验'
      rst[['p值']] = p
      rst[['p值解释']] = ifelse(p < 0.05, '两个变量之间存在显著关联', '两个变量之间无关联')
    }else if(flag == '校正卡方'){
      p = chisq.test(tb, correct = T)$p.value
      rst[['检验方法']] = '校正卡方'
      rst[['p值']] = p
      rst[['p值解释']] = ifelse(p < 0.05, '两个变量之间存在显著关联', '两个变量之间无关联')
    }else if(flag == '卡方'){
      p = chisq.test(tb, correct = F)$p.value
      rst[['检验方法']] = '卡方'
      rst[['p值']] = p
      rst[['p值解释']] = ifelse(p < 0.05, '两个变量之间存在显著关联', '两个变量之间无关联')
    }else{
      NULL
    }
    return(rst)
  }, error = function(err) {
    rst = list()
    rst[['检验方法']] = '计算出错'
    rst[['p值']] = NA
    rst[['p值解释']] = '-'
    return(rst)
  }
  )
}

getNormVars = function(contVarsList, strataVar, df){
  strat = as.factor(df[, strataVar])
  normList = c()
  for (var in contVarsList) {
    rst = normStatusTest(strat,
                         levels(strat), 
                         as.numeric(df[, var]))
    # print(rst)
    if (rst == '通过了方差齐性和正态性检测'){ normList = c(normList, var) }
  }
  return(normList)
}
p4Show = function(pValueList, digitsNum = 3){
  pValueList = as.numeric(pValueList)
  
  rst = NULL
  for (p in pValueList){
    if (is.na(p)){
      rst = c(rst, NA)
    }else if (p < 0.001){
      rst = c(rst, '<.001')
    }else if (p == 1){
      rst = c(rst, '1')
    }else{
      rst = c(rst, 
              str_pad(str_replace(round(p, digitsNum), '0.', '.'), 1+digitsNum, 'right', '0'))
    }
  }
  return(rst)
}

getCmpPValueRst = function(df, contVarsList, catVarsList, strataVar, pRound = 3){
  rstAll = NULL
  for (var in catVarsList) {
    cmpRst = compareGroupCat(as.factor(df[, strataVar]), as.factor(df[, var]))
    rstAll = rbind(rstAll,
                   c(var, '离散', cmpRst$检验方法, cmpRst$p值, cmpRst$p值解释))
  }
  for (var in contVarsList) {
    cmpRst = compareGroupCont(as.factor(df[, strataVar]), as.numeric(df[, var]))
    rstAll = rbind(rstAll,
                   c(var, '连续', cmpRst$检验方法, cmpRst$p值, cmpRst$p值解释))
  }
  rstAll = as.data.frame(rstAll)
  colnames(rstAll) = c('变量', '变量属性', '检验方法', 'p值', 'p值解释')
  rstAll[, 'p值4Show'] = p4Show(rstAll[, 'p值'], pRound)
  return(rstAll)
}
generateCmpPValueList = function(df, allVarsList, catVarsList, strataVarName, pRound = 3){
  # 入参检查
  if(length(setdiff(allVarsList, names(df))) != 0){
    warning(paste('所有参数中包含数据中不存在的变量，需剔除掉', paste(setdiff(allVarsList, names(df)), collapse = ', ')))
  }
  allVarsList = intersect(allVarsList, names(df))
  if (strataVarName %in% allVarsList){
    warning(paste('分层变量出现在了需要统计的所有变量中， 自动剔除掉', strataVarName))
  }
  allVarsList = setdiff(allVarsList, strataVarName)
  if(length(setdiff(catVarsList, allVarsList)) != 0){
    warning(paste('分类变量中包含所有变量中不存在的变量，需剔除掉', paste(setdiff(catVarsList, allVarsList), collapse = ', ')))
  }
  catVarsList = intersect(catVarsList, allVarsList)
  contVarsList = setdiff(allVarsList, catVarsList)
  message(paste('未被指定的变量自动当做连续变量处理：', paste(contVarsList, collapse = ', ')))
  pValueRst = getCmpPValueRst(df, contVarsList, catVarsList, strataVarName, pRound)
  return(pValueRst)
}

getTable1 = function(df, allVarsList, catVarsList, strataVarName, pRound = 3){
  # 入参检查
  if(length(setdiff(allVarsList, names(df))) != 0){
    warning(paste('所有参数中包含数据中不存在的变量，需剔除掉', paste(setdiff(allVarsList, names(df)), collapse = ', ')))
  }
  allVarsList = intersect(allVarsList, names(df))
  if (strataVarName %in% allVarsList){
    warning(paste('分层变量出现在了需要统计的所有变量中， 自动剔除掉', strataVarName))
  }
  allVarsList = setdiff(allVarsList, strataVarName)
  if(length(setdiff(catVarsList, allVarsList)) != 0){
    warning(paste('分类变量中包含所有变量中不存在的变量，需剔除掉', paste(setdiff(catVarsList, allVarsList), collapse = ', ')))
  }
  catVarsList = intersect(catVarsList, allVarsList)
  contVarsList = setdiff(allVarsList, catVarsList)
  message(paste('未被指定的变量自动当做连续变量处理：', paste(contVarsList, collapse = ', ')))
  for (i in contVarsList) { # 要手动将num型进行数据转换不然tableone会自动将包含字符串的num当做factor
    df[, i] = as.numeric(df[, i])
  }
  
  tableOne <- CreateTableOne(vars = allVarsList,
                             factorVars = catVarsList,
                             strata = strataVarName,
                             data = df,
                             addOverall = T,
                             test = F)
  
  tableOne4Show = print(tableOne, 
                        nonnormal = setdiff(contVarsList, getNormVars(contVarsList, strataVarName, df)),
                        showAllLevels = T,
                        missing = T, 
                        printToggle = F
  ) 
  tbRowName = row.names(print(tableOne, explain = F, showAllLevels = T, printToggle = F))
  
  pValueRst = getCmpPValueRst(df, contVarsList, catVarsList, strataVarName, pRound = pRound)
  for (idx in 1 : length(tbRowName)) {
    rowN = tbRowName[idx]
    if(idx == 1 & rowN == 'n'){ # 初始化
      pList = c(rep('', 3))
    }else{
      if (rowN %in% pValueRst[, '变量']){
        pDetail = pValueRst[which(pValueRst[, '变量'] == rowN), ]
        pList = rbind(pList, c(pDetail[, 'p值4Show'], pDetail[, 'p值解释'], pDetail[, '检验方法']))
      }else{
        pList = rbind(pList, rep('', 3))
        
      }
    }
  }
  colnames(pList) = c('p值', 'p值解释', '检验方法')
  rstFinal = cbind(tableOne4Show, pList )
  # print(rstFinal)
  return(rstFinal)
}


