# 判断df中某列是否有缺失值,是的话就返回true
isNaContain = function(df, feature){
  if( TRUE %in% is.na(df[,feature])){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
# 判断df中是否有缺失值,是的话就返回true
isDfContainNa = function(df){
  uniqueEle = unique(unlist(df))
  if( TRUE %in% is.na(uniqueEle)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
library('kknn')
fillNaWithKnn = function(df, featureList){
  # 将特征列表拆分为包含缺失值的特征组和不包含缺失值的特征组
  noNaFeature = featureList
  NaFeature = c()
  for(i in featureList){
    if(isNaContain(df,i) == TRUE){
      noNaFeature = setdiff(noNaFeature,i)
      NaFeature = c(NaFeature, i)
    }
  }

  #开始填充缺失值
  for (naFea in NaFeature){
    dfT = df[!is.na(df[naFea]),]
    dfV = df[ is.na(df[naFea]),]
    noNaFeatureTmp = noNaFeature
    for(i in noNaFeatureTmp){
      if(length(summary(as.factor(dfT[,i]))) == 1){
        cat('\n',i,'特征没有缺失值,但factor只有一个, 要过滤掉')
        noNaFeatureTmp = setdiff(noNaFeatureTmp,i)
      }
    }

    fo = paste(naFea,'~',paste(noNaFeatureTmp,collapse = '+'))
    knnModel = kknn(fo, dfT, dfV,
                    kernel = "triangular")
    # print(knnModel[["fitted.values"]])
    dfV[naFea] = knnModel[["fitted.values"]]
    df = rbind(dfT,dfV)
  }
  return(df)
}





library('pROC')
library('epiR')
getROCBestThreshold = function(df,prob,label){
  roc = roc(df[,label], df[,prob])
  bestThreshold = coords(roc,'best')[1,1]
  print(as.data.frame(coords(roc, "local maximas", ret=c("threshold", "sens", "spec", "ppv", "npv"),transpose = FALSE)),row.names = F)
  return(bestThreshold)
}
getTableSummary = function(df,feature,label){
  tableResult = table(df[,feature],df[,label])
  print
  #TP: True  Positive
  #FP: False Positive 假阳性: 金标准为负但检测为正
  #FN: False Negative 假阴性: 金标准为正但检测为负
  #TN: True  Negative
  #
  TP = tableResult[4]
  FP = tableResult[2]
  FN = tableResult[3]
  TN = tableResult[1]
  
  listTest = c(TP,FP,FN,TN)
  dat <- as.table(matrix(listTest, nrow = 2, byrow = TRUE))
  rval <- epi.tests(dat, conf.level = 0.95)
  # print(rval)
  # print(summary(rval))
  a = summary(rval)
  print(feature)
  cat('sen: ',round(a[3,2]*100,1),'% (',round(a[3,3]*100,1),'%-',round(a[3,4]*100,1),'%)', '\n' ,sep='')
  cat('spe: ',round(a[4,2]*100,1),'% (',round(a[4,3]*100,1),'%-',round(a[4,4]*100,1),'%)', '\n' ,sep='')
  cat('ppv: ',round(a[9,2]*100,1),'% (',round(a[9,3]*100,1),'%-',round(a[9,4]*100,1),'%)', '\n' ,sep='')
  cat('npv: ',round(a[10,2]*100,1),'% (',round(a[10,3]*100,1),'%-',round(a[10,4]*100,1),'%)', '\n' ,sep='')
  cat('acc: ',round(a[5,2]*100,1),'% (',round(a[5,3]*100,1),'%-',round(a[5,4]*100,1),'%)', '\n' ,sep='')
}
getModelResult = function(df,prob,label,bestThreshold = 0.5){
  roc1 = roc(df[,label], df[,prob])
  cat(paste('\nAUC: ',round(roc1$auc,4),' (',round(ci.auc(roc1)[1],4),'-',round(ci.auc(roc1)[3],4),')',sep = ''),'\n')
  
  cat('\n最优阈值:', bestThreshold, '\n')
  prob01 =  paste(prob,'_01',sep = '')
  df[,prob01] = ifelse(df[,prob] < bestThreshold,0,1)
  getTableSummary(df, prob01, label)
  return(df)
}

continuousDiscrete = function(df,feature,label,bestThreshold=0.5){
  roc = roc(df[,label], df[,feature])
  cat('\n',feature,'的最优阈值:', bestThreshold, '\n')
  prob01 =  paste(feature,'_01',sep = '')
  df[,prob01] = as.factor(ifelse(df[,feature] < bestThreshold,0,1))
  return(df)
}

getGlmCoefExt = function(coef){
  su = coef
  su = cbind(su,su[,1]-su[,2]*1.96) #coef 95CI Lower
  su = cbind(su,su[,1]+su[,2]*1.96) #coef 95CI Upper
  su = cbind(su,exp(su[,1])) #OR
  su = cbind(su,exp(su[,5])) #OR 95CI Lower
  su = cbind(su,exp(su[,6])) #OR 95CI Upper
  su = cbind(su,round(su[,4],2)) #P-value 保留两位
  su = cbind(su,round(su[,4],3)) #P-value 保留三位
  su = cbind(su,round(su[,1],2)) #coef 保留两位
  su = cbind(su,    paste(round(exp(su[,1]),1),'（',round(exp(su[,5]),1),'-',round(exp(su[,6]),1),'）',sep = 
                            '')) #OR 3in1
  su = cbind(row.names(su), su)
  colnames(su) = c('col', 'coef','se','zValue','pValue','coef_95CI_lower','coef_95CI_upper','OR'
                   ,'OR_95CI_lower','OR_95CI_upper', 'pValueIn2', 'pValueIn3', 'coefIn2', 'OR_3in1')
  return(as.data.frame(su,row.names = F))
}

# passStrategy = union并集(任一factor的p值小于阈值即通过)/intersection交集 (所有factor的p值小于阈值即通过)
getSingleLrPValue = function(df, feature, label, fixedColList = NULL, passStrategy = 'union', pValueThreshold = 0.1
                             , printDetail = F){
  if (TRUE %in% is.na(df[,feature])){
    meg = paste(feature, '中存在缺失值')
    stop(meg)
  }
  dfCal = df
  if (length(unique(dfCal[,feature])) == 1){
    cat('\n', feature, '该特征只有一个枚举值,不进行计算\n')
    return(FALSE)
    }else{
      fo = as.formula(ifelse(is.null(fixedColList)
                ,paste(label,'~',feature)
                ,paste(label,'~',paste(fixedColList,collapse = '+'),'+',feature)))
    
    glmModel = glm(
      data=dfCal,
      fo,
      family= binomial(link='logit'))
    coef = summary(glmModel)$coefficients
    coefExt = getGlmCoefExt(coef)
    if(printDetail){
      print(coefExt[-1,],row.names = F)
    }else{print(coefExt[-1,c(1,2,3,5)],row.names = F)
      }
    pValueList = as.numeric(coefExt[-1,'pValue'])
    if (passStrategy == 'union'){#并集
      if (TRUE %in% (pValueList < pValueThreshold)){
        return (TRUE)
      }else{return(FALSE)}
    }else if (passStrategy == 'intersection'){
      bool = 
      if(FALSE %in% (pValueList < pValueThreshold)){
        return(FALSE)
      }else {
        return(TRUE)
      }
    }else{
      print("输入的passStrategy有误")
      return(FALSE)
    }
  }
}

# 按照时间切分数据集
getDfLableRatio = function(df,label, minLabelInd){
  labelInd = unique(df[,label])
  if (length(labelInd) == 0){
    return(c(0,0,0))
  }
  rowsNum = dim(df)[1]
  rowsFilterNum = dim(df[which(df[,label] == minLabelInd),])[1]
  ratio = round(rowsFilterNum/rowsNum,2)
  result = c(rowsNum, rowsFilterNum, ratio)
  return(result)
}
splitDfByTmSummary = function(df,tm,format, label){
  tmpCol = 'tmForSplit'
  df[,tmpCol] = as.Date(df[,tm], format)
  tmList = unique(sort(df[,tmpCol]))
  minLabelInd = min(df[,label])
  cat("label计算的枚举类型为: ", minLabelInd, '\n')
  result = NULL
  
  for(i in tmList){
    cond = df[,tmpCol] <= i
    dfTrain = df[which(cond),]
    dfTest  = df[which(!cond),]
    tmListStr = unique(df[which(df[,tmpCol] == i),tm])
    resultTmp = c(
      tmListStr,i
      ,getDfLableRatio(df,label,minLabelInd)
      ,getDfLableRatio(dfTrain,label,minLabelInd)
      ,getDfLableRatio(dfTest,label,minLabelInd)
    )
    result = rbind(result,resultTmp)
    result = as.data.frame(result)
    colnames(result) = c('tm', 'tm2'
                         ,'dfOri_num','dfOri_label_num','dfOri_label_ratio'
                         ,'dfTrain_num','dfTrain_label_num','dfTrain_label_ratio'
                         ,'dfTest_num','dfTest_label_num','dfTest_label_ratio'
    )
    for(i in names(result)[-1]){
      result[,i] = as.numeric(result[,i])
    }
  }
  result['tainRatio'] = round(result$dfTrain_num/max(result$dfOri_num),2)
  
  return(result)
}

