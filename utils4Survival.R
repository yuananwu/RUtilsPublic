library('survival')
library('survminer')
library('rms')
library('pec') ##验证模型
library('timeROC')
library('CsChange')
library('dplyr')
library('stringdist')


# 将模型系数转换为10分值, 并加权求和
getNuiScore = function(df, coxModel, maxPoints = 10, printDetail = FALSE){
  if (class(coxModel) != 'coxph'){
    stop("cox model's source doesn't match 'coxph'")
  }
  fIn = names(coxModel$assign)
  coef = coxModel$coefficients
  coefUni = round(coef/max(abs(coef))*10, 0)
  scoreAll = rep(0, dim(df)[1])
  for(f in fIn){
    if (class(df[,f]) == 'factor'){
      lvl = levels(df[,f])
      for (l in lvl[-c(1)]){
        colName = paste(f, l ,sep = '')
        score = coefUni[which(names(coefUni) == colName)]
        if(printDetail){cat('特征: ', f, '; lvl: ', l, '; 得分: ', score, '\n')}
        scoreAll = scoreAll + ifelse(df[, f] == l, score, 0)
      }
    }else if (class(df[,f]) == 'numeric'){
      score = coefUni[f]
      if(printDetail){cat('特征: ', f, '; 数值型; 得分: ', score, '\n')}
      scoreAll = scoreAll + score * df[,f]
    }
  }
  return(scoreAll)
}

# 将coxph模型的结果组合成一个vector
getCoxphCoef = function(coef = NULL, exec = T){
  if (exec == T){
    c = coef[1]
    se = coef[3]
    coef4Show = round(c, 4)
    coefLow95CI = round(c-1.96*se,4)
    coefUpper95CI = round(c+1.96*se,4)
    HR = round(exp(c),4)
    HRLow95CI = round(exp(c-1.96*se),4)
    HRUpper95CI = round(exp(c+1.96*se),4)
    pValue = round(coef[5],4)
    rst = c(coef4Show, coefLow95CI, coefUpper95CI, HR, HRLow95CI, HRUpper95CI, pValue)
  }else{
    rst = c('-', '-', '-', '-', '-', '-', '-')
  }
}
#拓展一下模型结果
getCoxphCoefExt = function(df, feature, label, fixedColList = NULL, ignoreNa = T){
  if (any(is.na(df[, feature]))){
    cat('\n', feature, '中存在缺失值', '\n')
    if(!ignoreNa){
      stop("特征中存在缺失值, 同时也设定了不能忽略缺失值, 请检查数据")
    }
  }
  dfCal = df[!is.na(df[, feature]),]#去除入组特征是空值的样本
  
  restDf = NULL
  if (length(unique(dfCal[,feature])) == 1){
    msg = '该特征只有一个枚举值,不进行计算'
    rst = c(feature, msg, '-', getCoxphCoef(exec = F))
  }else{
    coef = summary(coxph(as.formula(paste(label, paste(c(fixedColList, feature), collapse = '+'), sep = '~')), data = dfCal))$coefficients
    
    coefRstAll = NULL
    if(class(df[, feature]) == 'factor'){
      for (lvl in levels(df[, feature])[-c(1)]){
        rowname = paste(feature, lvl, sep = '')
        coefRst = c(feature, '正常', lvl, getCoxphCoef(coef = coef[rowname, ], exec = T))
        coefRstAll = rbind(coefRstAll, coefRst)
      }
    }else{
      coefRstAll = c(feature, '正常', 'numeric', getCoxphCoef(coef = coef[feature, ], exec = T))
    }
    return(coefRstAll)
  }
}
getSingleCoxPValueV2 = function(df, featureList, label, fixedColList = NULL, passStrategy = 'union', pValueThreshold = 0.1, ignoreNa = T, paralCoreNum = 1){
  
  if (paralCoreNum == 1){# 不开启并行计算
    rst = NULL
    for (feature in featureList){
      rst = rbind(rst, getCoxphCoefExt(df, feature, label, fixedColList, ignoreNa))
    }
  }else{
    if(!(paralCoreNum > 1 & paralCoreNum <= 8)){
      warning('并行计算core 数建议介于 2 到 8之间')
    }
    library(foreach)
    library(doParallel)
    registerDoParallel(cores = paralCoreNum) #注册并行计算核心数
    # 并行计算
    rstTmp = foreach(i = featureList) %dopar% {
      getCoxphCoefExt(df = df, feature = i, label, fixedColList, ignoreNa)
    }
    rst = NULL
    for(i in rstTmp){
      rst = rbind(rst, i)
    }
  }

  rst = as.data.frame(rst, row.names = F)
  colnames(rst) = c('feature', 'calculate_status', 'levels', 'coef', 'coef_95CI_lower', 'coef_95CI_upper', 'HR', 'HR_95CI_lower', 'HR_95CI_upper', 'p_value')
  
  rst[, 'p_status'] = ifelse(as.numeric(rst[, 'p_value']) < pValueThreshold, 1, 0)
  rst[which(is.na(rst[, 'p_status'])), 'p_status'] = 0 #填充缺失值为0, 也即不通过p值检验状态
  pGroup = rst %>% group_by(feature) %>% dplyr::summarise(lvlNum = dplyr::n(),pStatusNum = sum(p_status))
  if (passStrategy == 'union'){ #并集, 任一lvl通过就通过
    validList = as.vector(unlist(pGroup %>% filter(pStatusNum >= 1) %>% dplyr::select(feature)))
  }else if(passStrategy == 'intersection'){#交集, 所有lvl通过才通过
    validList = as.vector(unlist(pGroup %>% filter(pStatusNum == lvlNum) %>% dplyr::select(feature)))
  }else{
    stop('输入的passStrategy有误')
  }
  rst[, 'p_status'] = ifelse(rst[, 'feature'] %in% validList, 1, 0)
  return(rst)
}


# https://zyicu.cn/?p=13869
library('survivalROC')
getSurvivalCutoff = function(df, tm, status, variable, method = 'km', rocPredTm = 0){
  if (method == 'km'){
    cutoff = surv_cutpoint(df, time = tm, event = status, variables = c(variable))
    print(cutoff)
    return(cutoff[["cutpoint"]]['cutpoint'])
  }else if (method == 'roc'){
    if(rocPredTm == 0){
      stop("method 'roc' need to specify the predict time")
    }
    #确定DEPDC1的最佳cutoff值，构建生存ROC函数（survivalROC）
    rocM = survivalROC(Stime = df[, tm],    #生存时间=time
                     status = df[, status], #生存状态=event
                     marker = df[, variable], #需要分析的变量
                     predict.time = mean(df[, tm]),      #预测时间
                     method="KM")          #使用生存分析KM计算
    print(paste('auc: ', rocM$AUC, sep = ''))
    #约登指数最大对应的值为最佳cutoff
    cutoff = rocM$cut.values[which.max(rocM$TP-rocM$FP)]
    return(cutoff)
  }else{
    return("")
  }

}

simplePlotKm = function(df, label, fea = 1){
  fit <- survfit(as.formula(paste(label, fea, sep = '~')),data = df) 
  ggplot = ggsurvplot(fit, data = df, pval = T, risk.table = T)
  return(ggplot)
}


coxRst2Txt = function(rstList, txtPath){
  capture.output(rstList, file = txtPath)
}

# 根据data属性分别计算相关系数并最终返回一个matrix
getCorMat = function(df, cols){
  getCor = function(df, col1, col2){
    numFlagMid = c(is.numeric(df[, col1]), is.numeric(df[, col2]))
    numFlag = length(which(numFlagMid))
    if (numFlag == 2){
      method = 'pearson'
    }else if (numFlag == 1){
      method = 'kendall'
    }else if (numFlag == 0){
      method = 'spearman'
    }else{
      method = 'error'
    }
    return(cor(as.numeric(df[, col1]), as.numeric(df[, col2]), method = method))
  }
  nCol = length(cols)
  corMat = matrix(rep(0, nCol**2), nrow = nCol, dimnames = list(cols, cols))
  for (i in 1: (nCol - 1 )) {
    j = i + 1
    col1 = cols[i]
    while( j <= nCol){
      col2 = cols[j]
      corMat[col1, col2] = getCor(df, col1, col2)
      j = j + 1 
    }
  }
  return(corMat)
}
# 特征01化
scale01 = function(v){
  return((v-min(v))/(max(v)-min(v)))
}
#获取相关性剔除特征的结果
getCorFilterRst = function(df, colsIn, survLabel){
  corMat = getCorMat(df, colsIn)
  relatePairedIdx = which(corMat > 0.6, arr.ind = T)
  rstList = NULL
  for(pair in 1 : nrow(relatePairedIdx)){
    idx1 = relatePairedIdx[pair, 1]
    idx2 = relatePairedIdx[pair, 2]
    c1 = colsIn[idx1]
    c2 = colsIn[idx2]
    md = step(coxph(as.formula(paste(survLabel, '~', paste(c(c1, c2),collapse = '+'))), data = df), direction = 'backward', trace = F)
    ag = md$assign
    if(length(ag) == 1){
      method = '后向cox保留'
      remainCol = names(ag)
    }else{ #后向cox未能剔除掉其中任意一个
      coef = summary(md)$coefficients
      fn = names(which.min(coef[, 5]))#最小p值对应的特征名字
      #通过计算Levenshtein距离获取到正确的
      method = '后向cox未能剔除，取p值最小的特征保留'
      remainCol = getSimiFeatureName(fn, c(c1, c2))
    }
    rstList = rbind(rstList,
                    c(colsIn[idx1], colsIn[idx2], corMat[idx1, idx2], method, remainCol, setdiff(c(colsIn[idx1], colsIn[idx2]), remainCol)))
  }
  rstList = as.data.frame(rstList)
  colnames(rstList) = c('特征1', '特征2', '相似度', '特征保留方法', '保留特征', '剔除掉的特征')
  return(rstList)
}
# 共线性分析
getVifRst = function(df, cols, saveLogPath = NULL, survLabel){
  delCols = c()
  itor = 1
  vifExistFlag = T # 是否还存在共线性的特征，如果为F就不再执行迭代
  vifRstOutput = list()
  while (itor <= 100 & vifExistFlag) {
    mdCols = setdiff(cols, delCols)
    fit = rms::cph(
      data=df
      ,as.formula(paste(survLabel,'~',paste(mdCols,collapse = '+')))
    )
    vifV = rms::vif(fit)
    if (max(vifV) > 2){ 
      #有共线性特征存在
      vifRstOutput[[paste(itor, '轮迭代的vif值')]] = as.data.frame(vifV)
      #去掉最大vif对应的特征
      delCols = c(delCols, getSimiFeatureName(names(which.max(vifV)), mdCols))
    }else{     
      vifRstOutput[['====无特征的vif大于阈值2, 停止迭代, 现阶段各特征的vif值为：']] = as.data.frame(vifV)
      vifRstOutput[['共线性检测删掉的特征有: ']] = paste(delCols, collapse = ', ')
      vifExistFlag = F
    }
    itor = itor + 1
  }
  if (itor == 100){
    warning('迭代了100次结果中还是有共线性指标超过阈值的特征， 请手动检查')
  }
  if(is.null(saveLogPath)){
    print(vifRstOutput)
  }else{
    coxRst2Txt(vifRstOutput, saveLogPath)
  }
  return(delCols)
}
