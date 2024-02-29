


#输入 : df seed 是否分层,分层字段


getCrossValidationIndex = function(df, cvFold, randSeed, strataFlag, strataFeature = ''){
  set.seed(randSeed)
  oriCols = names(df)
  dfRowNum = dim(df)[1]
  joinFeature = 'uniqueIDFuncInner'
  df[joinFeature] = 1:dfRowNum
  cvResultColName = 'cvFlag'
  
  if (strataFlag == TRUE){ #分层抽样
    dfCVResult = as.data.frame(matrix(nrow=0,ncol=2)) #创建一个2列的空对象
    colnames(dfCVResult) = c(joinFeature, cvResultColName)
    for ( i in unlist(unique(df[strataFeature]))){
      dfTmp = getCVIndex(df[df[strataFeature] == i,]$uniqueIDFuncInner,cvFold)
      dfCVResult = rbind(dfCVResult,dfTmp)
    }
  }else{ #非分层抽样
    dfCVResult = getCVIndex(df[joinFeature],cvFold)
  }
  finalDf = merge(df,dfCVResult,by = joinFeature)
  cols = c(oriCols,cvResultColName)
  return(finalDf[,cols])
}

#获取cv切分每一次的数据量级
getFoldSize = function(dfRowNum,fold){
  if (dfRowNum == 0){
    warning("data's rows is empty")
  }
  if (dfRowNum <= fold){
    return(c(rep(1,dfRowNum),rep(0,fold-dfRowNum)))
  }else {
    partitionNum = round(dfRowNum/fold)
    sizeListTmp = rep(partitionNum,fold)
    plusNum = ifelse(sum(sizeListTmp) >= dfRowNum,-1,1)
    plusSize = abs(sum(sizeListTmp)-dfRowNum)
    sizeListFix = c(rep(plusNum,plusSize),rep(0,fold-plusSize))
    return(sizeListTmp + sizeListFix )
  }
}

#传入一个index列表和切分份数,返回唯一每个id对应的一个cv值
getCVIndex = function(index,fold){
  dfCV = data.frame(index,rep(0,length(index)))
  colnames(dfCV) = c("uniqueIDFuncInner", "cvFlag")
  for (i in 1: fold) {
    indexList = dfCV[dfCV$cvFlag == 0,]$uniqueIDFuncInner
    sampleIndexList = sample(x=indexList,size=getFoldSize(length(index),fold)[i],replace=F)
    dfCV[dfCV$uniqueIDFuncInner %in% sampleIndexList,]$cvFlag = i
  }
  return(dfCV)
}













