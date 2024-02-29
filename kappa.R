library("raters")
#将原始数据转化为频次数据
freq_stat_2levels = function(df){#数据转换,每一个样本转换为所有评分者结果的0,1频次集合
  freq_vec = vector(mode="numeric",length=2*nrow(df))
  i = 1
  while(i <= nrow(df)){
    j = 1
    while(j <= ncol(df)){
      if( df[i,j] == 0 ){
        freq_vec[2*i-1] = freq_vec[2*i-1]+1
      }else{
        freq_vec[2*i] = freq_vec[2*i]+1
      }
      j = j+1
    }
    i = i+1
  }
  return(matrix(freq_vec,ncol = 2,byrow = T))
}

#数据转换,每一个样本转换为所有评分者结果的0,1,2频次集合
freq_stat_3levels = function(df){
  freq_vec = vector(mode="numeric",length=3*nrow(df))
  i = 1
  while(i <= nrow(df)){
    j = 1
    while(j <= ncol(df)){
      if( df[i,j] == 3 ){
        freq_vec[3*i-2] = freq_vec[3*i-2]+1
      }else if( df[i,j] == 4 ){
        freq_vec[3*i-1] = freq_vec[3*i-1]+1
      }else {
        freq_vec[3*i  ] = freq_vec[3*i  ]+1
      }
      j = j+1
    }
    i = i+1
  }
  return(matrix(freq_vec,ncol = 3,byrow = T))
}


fleissKappa = function(df){
  # Monte Carlo检验方法用来计算Fleiss Kappa的95CI,适用于小数据集
  #检验方法选择参考: https://cran.r-project.org/web/packages/raters/raters.pdf
  f_k_rsl = concordance(df,test="MC",B=100,alpha = 0.05)
  print(paste(
    " Fleiss Kappa :", round(f_k_rsl$Fleiss[1],3),
    " (",           round(f_k_rsl$Fleiss[2],3),
    "-",               round(f_k_rsl$Fleiss[3],3), ") " ,sep = ""))
}

#这是多raters&多levels的计算方法
fleissKappaMultiLevels = function(df){
  set.seed(42)
  result_a = wlin.conc(df,test="MC",B=25)
  print(paste(
    " Kappa :", round(result_a[1],3),
    " (",           round(result_a[3],3),
    "-",               round(result_a[4],3), ") " ,sep = ""))
}
