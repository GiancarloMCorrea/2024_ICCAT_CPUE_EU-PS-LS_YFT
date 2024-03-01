influ_new = function(model, data=NULL, response=NULL, focus=NULL) {
  
  instance = list(
    model = model,
    data = data,
    response = response,
    focus = focus
  )
  
  instance = influ_addInit(instance)
  return(instance)
  
}

influ_addInit = function(obj) {
  
  # If no data was supplied....
  if(is.null(obj$data)){
    # See if the model has a data variable
    obj$data = obj$model$data
    if(is.null(obj$data)){
      # See if it is available as a model attribute
      obj$data = attr(obj$model,"data")
    }
    if(is.null(obj$data)){
      stop("You need to provide the data that was fitted to for this type of model")
    }
  }
  
  obj$terms = attr(obj$model$terms, "term.labels")
  if(is.null(obj$response)) obj$response = names(obj$model$model)[1]
  if(is.null(obj$focus)) obj$focus = obj$terms[1]
  
  obj$labels = list()
  for(term in obj$terms) obj$labels[term] = term
  
  obj$orders = list()
  for(term in obj$terms) obj$orders[term] = 'asis'
  
  #An alternative, but potentially buggy, way of getting response term
  #Obtain the response name. This seems the safest way, albeit long winded.
  #variablesList = attr(model$terms,"variables") #List of variables including the response
  #responseIndex = attr(model$terms,"response") #Index of the response variable in variablesList 
  #responseName = as.character(variablesList)[responseIndex+1]
  
  return(obj)
  
}

influ_coeffs = function(model, term){
  
  coeffs = coefficients(model)
  rows = substr(names(coeffs),1,nchar(term))==term
  out = c(0,coeffs[rows])
  return(out)
  
}


influ_ses = function(model, term){
  type = class(model)[1]
  if(type=='glm'|type=='negbin'){
    #The "cov.scaled" member of a glm object is the estimated covariance matrix of 
    #the estimated coefficients scaled by dispersion"    
    V = summary(model)$cov.scaled
  }
  else if(type=='survreg'){
    #The member "var" for a survreg object is the "the variance-covariance matrix for the parameters, 
    #including the log(scale) parameter(s)" 
    V = model$var
  }
  
  rows = substr(row.names(V),1,nchar(term))==term
  V = V[rows,rows]
  n = sum(rows)+1
  Q = matrix(-1/n, nrow=n, ncol=n-1)
  Q[-1,] = Q[-1,] + diag(rep(1,n-1))
  V0 = (Q%*%V) %*% t(Q)
  se = sqrt(diag(V0))
  return(se)
}


influ_effects = function(model, focus){
  
  coeffs = influ_coeffs(model, focus)
  out = exp(coeffs-mean(coeffs))
  return(out)
  
}

influ_addCalc = function(obj){
  #Get observed values
  observed = obj$model$model[,obj$response]
  if(class(observed)[1]=='Surv') observed = as.numeric(observed[,1])
  logged = substr(obj$response,1,4)=='log('
  #Create a data frame that is used to store various values for each level of focal term
  obj$indices = data.frame(level=levels(obj$model$model[,obj$focus]))
  #Add an unstandardised index
  if(logged | sum(observed<=0)==0){
    if(logged) log_observed = observed else log_observed = log(observed)
    # Calculate geometic mean
    obj$indices = merge(obj$indices, aggregate(list(unstan=log_observed),list(level=obj$model$model[,obj$focus]),mean))
    # Turn into a relative index
    obj$indices$unstan = with(obj$indices, exp(unstan-mean(unstan)))
  }else {
    # There are zeros (or even negative numbers) in the data so we cant calculate a geometric mean.
    # Use arithmetic mean instead
    obj$indices = merge(obj$indices,aggregate(list(unstan=observed),list(level=obj$model$model[,obj$focus]),mean))
    # Turn into a relative index
    obj$indices$unstan = with(obj$indices,unstan/mean(unstan))
  }
  #Add standardised index.
  coeffs = influ_coeffs(obj$model, obj$focus)
  ses = influ_ses(obj$model, obj$focus)
  base = mean(coeffs)
  obj$indices$stan = exp(coeffs-base)
  obj$indices$stanLower = exp(coeffs-base-2*ses)
  obj$indices$stanUpper = exp(coeffs-base+2*ses)
  
  #Create models with terms successively added
  obj$summary = NULL
  
  #TODO calculate influence statistics in this loop too
  for(termCount in 0:length(obj$terms)){
    if(termCount>0){
      term = obj$terms[termCount]
      #Update both formula and model
      newFormula = formula(obj$model)
      newFormula = update.formula(newFormula,formula(paste("~",paste(paste(obj$terms[1:termCount],collapse='+')))))
      model = update(obj$model,newFormula)
      #Get index for this model
      index = influ_effects(model, obj$focus)
      #Add column to .$indices
      obj$indices = cbind(obj$indices,index)
      #Give index column the right hand side of formula as name
      names(obj$indices)[ncol(obj$indices)] = if(termCount==1) term else paste('+',term)
    } else {
      term = 'intercept'
      model = update(obj$model,.~1)
    }
    
    type = class(model)[1]
    logLike =  switch(type,
                      survreg = model$loglik[2],
                      logLik(model)
    )
    fitted = switch(type,
                    survreg = predict(model,type='response'),
                    fitted(model)
    )
    
    #Sums-of-squared based R2 based on log(observed) and log(fitted)
    if(termCount==0) r2 = 0 else r2 = cor(observed,fitted)^2 # CHECK: original version if log(), but I think it is a mistake
    
    #Deviance pseudo-R2
    r2Dev = (model$null.deviance-model$deviance)/model$null.deviance
    if(length(r2Dev)==0) r2Dev = NA
    
    #Negelkerke pseudo-R2
    if(termCount==0) logLikeInterceptOnly = logLike
    n = length(observed)
    r2Negel = (1-exp((logLikeInterceptOnly-logLike)*(2/n)))/(1-exp(logLikeInterceptOnly*(2/n)))
    
    obj$summary = rbind(obj$summary,data.frame(
      term = term,
      k = length(coef(model)),
      logLike = logLike,
      aic = extractAIC(model)[2],
      r2 = r2,
      r2Dev = r2Dev,
      r2Negel = r2Negel
    ))
  }
  #R2 values presented as differences
  obj$summary = within(obj$summary,{
    k = c(1,diff(k))
    r2 = c(NA,diff(r2))
    r2Dev = c(NA,diff(r2Dev))
    r2Negel = c(NA,diff(r2Negel))
  })
  
  #Calculate predicted values
  preds = predict(obj$model,type='terms',se.fit=T)
  fit = as.data.frame(preds$fit)
  se.fit = as.data.frame(preds$se.fit)
  preds = cbind(fit,se.fit)
  names(preds) = c(paste('fit',names(fit),sep='.'),paste('se.fit',names(fit),sep='.'))
  obj$preds = cbind(obj$data, preds)
  #Calculate influences and statisitcs
  obj$influences = data.frame(level=levels(obj$model$model[,obj$focus]))
  overall = c(NA,NA) # NAs for null model and for focus term
  trend = c(NA,NA)
  for(term in obj$terms){
    if(term!=obj$focus){
      infl = aggregate(
        list(value = obj$preds[,paste('fit',term,sep='.')]),
        list(level = obj$preds[,obj$focus]),
        mean
      )
      overall = c(overall,with(infl,exp(mean(abs(value)))-1))
      trend = c(trend,with(infl,exp(cov(1:length(value),value)/var(1:length(value)))-1))
      names(infl) = c('level',term)
      obj$influences = merge(obj$influences,infl,all.x=T,by='level')
    }
  }
  #Insert statistics into summary table with NAs for NULL and focus terms
  obj$summary$overall = overall
  obj$summary$trend = trend
  
  return(obj)
}


influ_stanPlot <- function(obj){
  
  require(ggplot2)
  df_plot = obj$indices
  p1 = ggplot(data = df_plot, aes(x = as.integer(level), y = stan)) +
    geom_line() +
    geom_ribbon(aes(ymin = stanLower, ymax = stanUpper), linetype = 2, alpha = 0.2) +
    geom_point(aes(x = as.integer(level), y = unstan), color = 2) +
    ylab('Index') + xlab(obj$labels[[obj$focus]])
  return(p1)
  
}


influ_stepPlot <- function(obj, panels=T, setpar=T){

  startCol = 6
  cols = startCol:ncol(obj$indices)
  if(panels){
    if(setpar)par(mfrow=c(length(cols),1),mar=c(0,4,0,0),oma=c(4,1,1,1))
    levels = as.integer(obj$indices$level)
    xlim=range(levels)
    ylim=c(0,max(obj$indices[,cols]))
    for(col in cols){
      plot(NA,xaxt='n',las=1,ylab='Index',ylim=ylim,xlim=xlim)
      #abline(h=1,lty=2)
      for(prev in cols[1]:col){
        if(prev<col-1) lines(levels,obj$indices[,prev],col='grey',lwd=1.5)
        if(prev==col-1) lines(levels,obj$indices[,prev],lty=3,col='black',lwd=1.5)
        if(prev==col) points(levels,obj$indices[,prev],pch=16,cex=1.25,col='black',type='o')
      }
      legend('topleft',legend=names(obj$indices)[col],bty='n')
    }
    axis(side=1,at=xlim[1]:xlim[length(xlim)],labels=obj$indices$level)
  }
  else {
    levels = as.integer(obj$indices$level)
    xlim=range(levels)
    with(obj$indices,{
      plot(NA,ylim=c(0,max(obj$indices[,cols])),xlim=xlim,las=1,ylab='Index',xlab=obj$labels[[obj$focus]],xaxt='n')
      for(col in cols) lines(as.integer(level),obj$indices[,col],type='o',pch=col-startCol+1,cex=1.25)
      axis(side=1,at=xlim[1]:xlim[length(xlim)],labels=level)
      legend('top',legend=names(obj$indices)[cols],pch=cols-startCol+1,pt.cex=1.25,bty='n')
    })
  }
  
}


influ_influPlot = function(obj, panels=T, setpar=T){
  
  cols = 2:ncol(obj$influences)
  ylim=exp(range(obj$influences[,cols]))
  if(panels){
    if(setpar) par(mfrow=c(length(cols),1),mar=c(0,4,0,0),oma=c(4,1,1,1))
    levels = as.integer(obj$influences$level)
    xlim=range(levels)
    for(col in cols){
      plot(levels,exp(obj$influences[,col]),pch=16,cex=1.25,col='black',type='o',xaxt='n',las=1,ylab='Influence',ylim=ylim,xlim=xlim)
      abline(h=1,lty=2)
      legend('topleft',legend=names(obj$influences)[col],bty='n')
    }
    axis(side=1,at=xlim[1]:xlim[length(xlim)],labels=obj$influences$level)
  }
  else{
    levels = as.integer(obj$influences$level)
    xlim=range(levels)
    with(obj$influences,{
      plot(NA,ylim=ylim,xlim=xlim,las=1,ylab='Influence',xlab=obj$labels[[obj$focus]],xaxt='n')
      for(col in cols) lines(as.integer(level),exp(obj$influences[,col]),type='o',pch=col,cex=1.25)
      axis(side=1,at=xlim[1]:xlim[length(xlim)],labels=level)
      legend('top',legend=names(obj$influences)[cols],pch=cols,pt.cex=1.25,bty='n')
      abline(h=1,lty=2)
    })
  }
  
}


influ_stepAndInfluPlot <- function(obj){
  par(mfcol=c(ncol(obj$indices)-5,2),mar=c(0,5,0,0),oma=c(4,1,1,1))
  influ_stepPlot(obj = obj, setpar=F)
  #Create a blank plot in the top of the influences column for the row for the focus term.
  plot.new()
  influ_influPlot(obj = obj, setpar=F)
}


influ_cdiPlot <- function(obj, term, variable=NULL){
  par(oma=c(1,1,1,1),cex.lab=1.25,cex.axis=1.25)
  layout(matrix(c(1,1,0,2,2,3,2,2,3),3,3,byrow=TRUE))
  
  #Define levels of term on which coefficient and distribution plots will be based
  #This is done by search for each column name in the data as a whole word in the
  #each term. This allows for matching of terms like 'poly(log(var),3)' with 'var'
  if(is.null(variable)){
    for(name in names(obj$preds)){
      match = grep(paste('([^[:alnum:]_])+',name,'([^[:alnum:]_])+',sep=''),paste('(',term,')')) # ([^\\w])+ Means ignore any 'word' character (alphanumerics plus underscore)
      if(length(match)>0){
        variable = name
        break
      }
    }
  }
  if(is.null(variable)) stop('Unable to find a matching variable for term "',term,'". Please specify in the argument "variable".')
  levels = obj$preds[,variable]
  
  #Numeric terms are cut into factors
  if(is.numeric(levels)){
    breaks = pretty(levels,30)
    step = breaks[2]-breaks[1]
    labels = breaks+step/2
    breaks = c(breaks,breaks[length(breaks)]+step)
    levels = cut(levels,breaks,labels=labels,include.lowest=T)
  }
  
  #Reorder levels according to coefficients if necessary
  if(obj$orders[[term]]=='coef'){
    coeffs = aggregate(obj$preds[,paste('fit',term,sep='.')],list(levels),mean)
    names(coeffs) = c('term','coeff')
    coeffs = coeffs[order(coeffs$coeff),]
    levels = factor(levels,levels=coeffs$term,ordered=T)
  }
  
  #Coefficients
  coeffs = aggregate(obj$preds[,paste(c('fit','se.fit'),term,sep='.')],list(levels),mean)
  names(coeffs) = c('term','coeff','se')
  coeffs = within(coeffs,{
    lower = coeff-se
    upper = coeff+se
  })
  
  par(mar=c(0,5,3,0),las=1)
  with(coeffs,{
    xs = 1:max(as.integer(term))
    ylim = c(min(exp(lower)),max(exp(upper)))
    if(ylim[1]<0.5*min(exp(coeff))) ylim[1] = 0.5*min(exp(coeff))
    if(ylim[2]>2*max(exp(coeff))) ylim[2] = 2*max(exp(coeff))
    plot(as.integer(term),exp(coeff),xlim=range(xs),ylim=ylim,pch=2,cex=1.5,xaxt='n',ylab='',log='y')
    mtext('Coefficient',side=2,line=4,las=0,cex=0.8)
    abline(h=1,lty=2)
    abline(v=xs,lty=1,col='grey')
    segments(as.integer(term),exp(lower),as.integer(term),exp(upper))
    arrows(as.integer(term),exp(lower),as.integer(term),exp(upper),angle=90,code=3,length=0.05)
    axis(3,at=xs,labels=levels(term)[xs])
  })
  
  #Distribution
  distrs = aggregate(obj$preds[,1],list(levels,obj$preds[,obj$focus]),length)
  names(distrs) = c('term','focus','count')
  distrs = merge(distrs,aggregate(list(total=distrs$count),list(focus=distrs$focus),sum))
  distrs$prop = with(distrs,count/total)
  par(mar=c(5,5,0,0),las=1)
  xlab = obj$labels[[variable]]
  if(is.null(xlab)) xlab = variable
  ylab= obj$labels[[obj$focus]]
  if(is.null(ylab)) ylab = obj$focus
  with(distrs,{
    xs = 1:max(as.integer(term))
    ys = 1:max(as.integer(focus))
    plot(NA,xlim=range(xs),ylim=range(ys),xaxt='n',yaxt='n',xlab=xlab,ylab='')
    abline(v=xs,lty=1,col='grey')
    axis(1,at=xs,labels=levels(term)[xs])
    abline(h=ys,lty=1,col='grey')
    axis(2,at=ys,labels=levels(focus)[ys])
    mtext(ylab,side=2,line=4,las=0,cex=0.8)
    points(as.integer(term),as.integer(focus),cex=sqrt(prop)*12)
  })
  
  #Influence
  par(mar=c(5,0,0,3),las=1)
  ys = 1:nrow(obj$influences)
  with(obj$influences,{
    plot(NA,xlim=range(exp(get(term))),ylim=range(ys),yaxt='n',xlab='Influence')
    abline(h=ys,lty=1,col='grey')
    abline(v=1,lty=2)
    points(exp(get(term)),ys,type='o',pch=16,cex=1.8)
    axis(4,at=ys,labels=level)
  })
  
}


influ_cdiPlotAll <- function(obj, done = function(term){
  cat("cdiPlot for",term,". Press enter for next")
  scan()
}){
  for(term in obj$terms) {
    if(term!=obj$focus){
      influ_cdiPlot(obj = obj, term)
      done(term)
    }
  }
}
