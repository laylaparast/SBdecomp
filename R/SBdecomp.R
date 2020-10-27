pred.smooth <-function(zz,zi.one, bw=NULL,y1, weight = NULL) { 	
   	val.unique = unique(zz)
   	if(is.factor(zz) |  sum((val.unique %in% c(0,1)) == FALSE) == 0) {
  		val.unique = unique(zz)
  		num.unique = length(val.unique)
  		est = vector(length = length(zi.one))
  		for(mm in 1:num.unique) {
  			est[zi.one == val.unique[mm]] = mean(y1[zz == val.unique[mm]])
  		}
  		est[est >=0.9999] = 0.9999 
     	est[est <= 0.0001] = 0.0001
     	return(est)
  		
  	}
  if(!(is.factor(zz) |  sum((val.unique %in% c(0,1)) == FALSE) == 0)) {
  	if(is.null(bw)) { bw = bw.nrd(zz)}
  	if(is.null(weight)) {weight = rep(1,length(y1))}
    est = sum(weight*Kern.FUN(zz,zi.one,bw=bw)*y1)/sum(weight*Kern.FUN(zz,zi.one,bw=bw))
 	if(sum(is.na(est))>0){
       	c.mat = cbind(zi.one, est)
    	for(o in 1:length(est)) {
    		if(is.na(est[o])){
    			distance = abs(zi.one - zi.one[o])
    			c.temp = cbind(c.mat, distance)
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where predication is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			est[o] = new.est[1]   #in case there are multiple matches
    	}
  }}
     est[est >=0.9999] = 0.9999 
     est[est<= 0.0001] = 0.0001
    	  return(est) }
  }



VTM<-function(vc, dm){
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }
    
 
Kern.FUN <- function(zz,zi,bw) 
  { 
    out = (VTM(zz,length(zi))- zi)/bw
	dnorm(out)/bw
           
  }

sbdecomp = function(outcome, treatment, confounders, data = NULL, type = "inclusion", estimation = "parametric", Bonly = T, balance = T, n.trees = 20000, interaction.depth = 4, shrinkage = 0.005, verbose = FALSE,  stop.method = c("es.max"), cv.folds = 0, standard.error = F, boot.rep=500)
   {#confounders must be in a data.frame with names
	
	if(is.character(outcome) | is.character(treatment) | is.character(confounders)) {
		if(is.null(data)) {stop("Data not supplied.")}	
	}
	if(is.character(outcome)) {outcome.temp = data[,names(data) == outcome, drop = TRUE]; outcome = outcome.temp}
	if(is.character(treatment)) {treatment.temp = data[,names(data) == treatment, drop = TRUE]; treatment = treatment.temp}
	if(is.character(confounders)) {confounders.temp = data[,!is.na(match(names(data), confounders)), drop = TRUE]; confounders = confounders.temp}
	
	data.check=cbind(outcome, treatment, confounders)
	if(sum(is.na(data.check)) > 0) {
		sss = apply(is.na(data.check),2,sum)
		sss = sss[sss>0]
		sss = paste(names(sss), collapse=" ")
		message(paste("Missing values in variable(s):", sss))
		ttt = apply(is.na(data.check),1,sum)
		ttt = which(ttt > 0)
		outcome = outcome[-ttt]
		treatment = treatment[-ttt]
		confounders = confounders[-ttt,]
 		message(paste("Complete cases only used in analysis; ", length(ttt), " cases removed."))
	}
	names.covariates = names(confounders)
	ncovariates = length(names.covariates)
	is.f.vec = vector(length = ncovariates)
	for(kk in 1:ncovariates){
		is.f.vec[kk] = is.factor(confounders[,kk])
	}
	if(sum(is.f.vec) ==0 ) {message("No confounders identified as factors. If this is incorrect, please check that variables that should be factors are coded as factors in the supplied arguments.")}
	if(sum(is.f.vec) > 0) {message(paste("The following variable(s) were identified as factors:", names.covariates[is.f.vec]))}
	delta.naive = mean(outcome[treatment==1]) - mean(outcome[treatment==0])
	p.value.delta.naive = t.test(outcome[treatment==1], outcome[treatment==0])$p.value
	ci.delta.naive = as.vector(t.test(outcome[treatment==1], outcome[treatment==0])$conf.int)
	great.16 = which(nchar(names.covariates) > 16)
	if(length(great.16) > 0) {
		for(tt in great.16) {
			message(paste("Note: Variable name for ",names.covariates[tt]," shortened to ", substr(names.covariates[tt],1,16),sep="" ))
			names(confounders)[tt] = substr(names.covariates[tt],1,16)
			names.covariates[tt] = substr(names.covariates[tt],1,16)
		}
	}
	data.use = as.data.frame(cbind(outcome, treatment, confounders))
	
	if(!(estimation == "nonparametric" | estimation == "parametric"))    {stop("Estimation method must be either parametric or nonparametric.")}
   	if(!(type == "removal" | type == "inclusion"))    {stop("Type must be either removal or inclusion.")}

if(estimation == "nonparametric") {
   	ps.all=twang::ps(as.formula(paste("treatment ~", paste(names.covariates, collapse="+"))), data=data.use,
      n.trees = n.trees, interaction.depth = interaction.depth, verbose = verbose, estimand = "ATE",
      stop.method = stop.method, shrinkage = shrinkage, cv.folds=cv.folds)
      weights.np = twang::get.weights(ps.all,stop.method=stop.method)
      delta.ps =  weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - weighted.mean(outcome[treatment==0], weights.np[treatment== 0])
      	data.use$weights.np = weights.np
	  	design.ps <- svydesign(ids=~1, weights=~weights.np, data=data.use)
		glm1 <- svyglm(outcome ~ treatment, design=design.ps)
		p.value.delta.fully.adjusted = summary(glm1)$coef[2,4]
		ci.delta.fully.adjusted = as.vector(confint(glm1)[2,])
      	lambda.all = delta.naive-delta.ps
	  
	  balance.naive.mean = mean(abs(twang::bal.table(ps.all)$unw$std.eff.sz))
	  balance.naive.max = max(abs(twang::bal.table(ps.all)$unw$std.eff.sz))
	  balance.fully.adjusted.mean = mean(abs(twang::bal.table(ps.all)$es.max.ATE$std.eff.sz))
	  balance.fully.adjusted.max = max(abs(twang::bal.table(ps.all)$es.max.ATE$std.eff.sz))
	  
	  if(type == "removal")	{
	  	message("Note: the nonparametric estimation procedure is running; this procedure may be very time intensive (on the scale of several minutes to hours depending on the size of your dataset).")
		balance.ex.mean = vector(length = ncovariates)
		balance.ex.max = vector(length = ncovariates)
		lambda.each.ex.np = vector(length = ncovariates)
		B.each.ex.np = vector(length =  ncovariates)
		delta.each.ex.np = vector(length =  ncovariates)
	  	for(j in 1:ncovariates) {
			#single confounder removal
			names.covariates.remove = names.covariates[-j]
			if(length(names.covariates.remove) >1) {
			ps.all=twang::ps(as.formula(paste("treatment ~",  paste(names.covariates.remove, collapse="+"))),
   			data=data.use,n.trees = n.trees, interaction.depth = interaction.depth, verbose = verbose, estimand 			= "ATE", stop.method = stop.method, shrinkage = shrinkage, cv.folds=cv.folds)
			weights.np = twang::get.weights(ps.all,stop.method=stop.method)
			delta.ps.vec.np.temp = weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - weighted.mean(outcome[treatment==0], weights.np[treatment== 0])
			
			delta.each.ex.np[j] = delta.ps.vec.np.temp
			lambda.each.ex.np[j] = delta.ps-delta.ps.vec.np.temp	
			balance.ex.mean[j] = mean(abs(twang::bal.table(ps.all)$es.max.ATE$std.eff.sz))
			balance.ex.max[j] = max(abs(twang::bal.table(ps.all)$es.max.ATE$std.eff.sz))
			}
			if(length(names.covariates.remove) ==1) {
			single.var = data.use[,names(data.use) == names.covariates.remove]
			pred = sapply(single.var,pred.smooth,zz=single.var, y1=treatment)
			weights.np = treatment*(1/pred)+(1-treatment)*(1/(1-pred))
			delta.ps.vec.np.temp  = weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - weighted.mean(outcome[treatment==0], weights.np[treatment== 0])
			delta.each.ex.np[j] = delta.ps.vec.np.temp
			lambda.each.ex.np[j] = delta.ps-delta.ps.vec.np.temp
			dat.hold = data.frame(treatment = treatment, single.var=single.var)
			hold.weights = twang::dx.wts(x=weights.np, data = dat.hold, vars = "single.var", treat.var = "treatment", x.as.weights = TRUE, estimand = "ATE")
			balance.ex.mean[j] = mean(abs(twang::bal.table(hold.weights)[[2]]$std.eff.sz))
			balance.ex.max[j] = max(abs(twang::bal.table(hold.weights)[[2]]$std.eff.sz))
			}
	
		}
		B.sum = sum(abs(lambda.each.ex.np))
		B.each.ex.np = abs(lambda.each.ex.np)/B.sum
		lambda.each.ex.np = as.data.frame(t(lambda.each.ex.np))
		delta.each.ex.np = as.data.frame(t(delta.each.ex.np))
		names(lambda.each.ex.np) = names.covariates
		row.names(lambda.each.ex.np) = "lambda"
		B.each.ex.np = as.data.frame(t(B.each.ex.np))
		names(B.each.ex.np) = names.covariates
		row.names(B.each.ex.np) = "B"
		balance.ex.mean = as.data.frame(t(balance.ex.mean))
		names(balance.ex.mean) = names.covariates
		row.names(balance.ex.mean) = "balance.mean"
		balance.ex.max = as.data.frame(t(balance.ex.max))
		names(balance.ex.max) = names.covariates
		row.names(balance.ex.max) = "balance.max"
		names(delta.each.ex.np) = names.covariates
		row.names(delta.each.ex.np) = "delta"
		results = list("delta.naive" = delta.naive, "p.value.delta.naive" = p.value.delta.naive, "conf.int.delta.naive" = ci.delta.naive, "delta.fully.adjusted" = delta.ps, "p.value.delta.fully.adjusted" = p.value.delta.fully.adjusted, "conf.int.delta.fully.adjusted" = ci.delta.fully.adjusted, "B" = B.each.ex.np)
		if(!Bonly) {
			results = c(results,list("estimated.selection.bias" = lambda.all, "lambda" = lambda.each.ex.np, "delta.each" = delta.each.ex.np))
		}
		if(balance){
			results = c(results,list("balance.naive.mean" = balance.naive.mean, 
	  "balance.naive.max" = balance.naive.max, "balance.fully.adjusted.mean" = balance.fully.adjusted.mean,  "balance.fully.adjusted.max" = balance.fully.adjusted.max, "balance.mean" = balance.ex.mean, "balance.max" = balance.ex.max))
		}
		}
   	
   	if(type == "inclusion")	{
		balance.in.mean = vector(length =  ncovariates)
   		balance.in.max = vector(length =  ncovariates)
   		lambda.each.in.np = vector(length = ncovariates)
   		B.each.in.np = vector(length =  ncovariates)
		delta.each.in.np = vector(length = ncovariates)
		
	  	for(j in 1:ncovariates) {
			#single confounder removal
			names.covariates.one = names.covariates[j]
			single.var = data.use[,names(data.use) == names.covariates.one]
			pred = sapply(single.var,pred.smooth,zz=single.var, y1=treatment)
			weights.np = treatment*(1/pred)+(1-treatment)*(1/(1-pred))
			delta.ps.vec.np.temp  = weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - weighted.mean(outcome[treatment==0], weights.np[treatment== 0])
			delta.each.in.np[j] = delta.ps.vec.np.temp
			lambda.each.in.np[j] = delta.naive-delta.ps.vec.np.temp
			dat.hold = data.frame(treatment = treatment, single.var=single.var)
			hold.weights = twang::dx.wts(x=weights.np, data = dat.hold, vars = "single.var", treat.var = "treatment", x.as.weights = TRUE, estimand = "ATE")
			balance.in.mean[j] = mean(abs(twang::bal.table(hold.weights)[[2]]$std.eff.sz))
			balance.in.max[j] = max(abs(twang::bal.table(hold.weights)[[2]]$std.eff.sz))
		}
		B.sum = sum(abs(lambda.each.in.np))
		B.each.in.np = abs(lambda.each.in.np)/B.sum 
		lambda.each.in.np = as.data.frame(t(lambda.each.in.np))
		delta.each.in.np = as.data.frame(t(delta.each.in.np))
		names(lambda.each.in.np) = names.covariates
		row.names(lambda.each.in.np) = "lambda"
		B.each.in.np = as.data.frame(t(B.each.in.np))
		names(B.each.in.np) = names.covariates
		row.names(B.each.in.np) = "B"
		balance.in.mean = as.data.frame(t(balance.in.mean))
		names(balance.in.mean) = names.covariates
		row.names(balance.in.mean) = "balance.mean"
		balance.in.max = as.data.frame(t(balance.in.max))
		names(balance.in.max) = names.covariates
		row.names(balance.in.max) = "balance.max"
		names(delta.each.in.np) = names.covariates
		row.names(delta.each.in.np) = "delta"
		results = list("delta.naive" = delta.naive, "p.value.delta.naive" = p.value.delta.naive, "conf.int.delta.naive" = ci.delta.naive, "delta.fully.adjusted" = delta.ps, "p.value.delta.fully.adjusted" = p.value.delta.fully.adjusted, "conf.int.delta.fully.adjusted" = ci.delta.fully.adjusted, "B" = B.each.in.np)
		if(!Bonly) {
			results = c(results,list( "estimated.selection.bias" = lambda.all, "lambda" = lambda.each.in.np, delta.each = delta.each.in.np))
		}
		if(balance){
			results = c(results,list("balance.naive.mean" = balance.naive.mean, 
	  "balance.naive.max" = balance.naive.max, "balance.fully.adjusted.mean" = balance.fully.adjusted.mean,  "balance.fully.adjusted.max" = balance.fully.adjusted.max, "balance.mean" = balance.in.mean, "balance.max" = balance.in.max))
		}
	}
	}
	if(estimation == "parametric") {
	      
      p.model = glm(as.formula(paste("treatment ~", paste(names.covariates, collapse="+"))),
      data=data.use, family = binomial)
      pred = predict(p.model, type = "response")
      weights.np = treatment*(1/pred)+(1-treatment)*(1/(1-pred))
      delta.ps =  weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - weighted.mean(outcome[treatment==0], weights.np[treatment== 0])
      data.use$weights.np = weights.np
	  	design.ps <- svydesign(ids=~1, weights=~weights.np, data=data.use)
		glm1 <- svyglm(outcome ~ treatment, design=design.ps)
		p.value.delta.fully.adjusted = summary(glm1)$coef[2,4]
		ci.delta.fully.adjusted = as.vector(confint(glm1)[2,])

      lambda.all = delta.naive-delta.ps
	  
	  dat.hold = as.data.frame(cbind(treatment, data.use[,  names(data.use) %in% names.covariates]))
	  hold.weights = twang::dx.wts(x=weights.np, data = dat.hold, vars = names.covariates, treat.var = "treatment", x.as.weights = TRUE, estimand = "ATE")
	  balance.naive.mean = mean(abs(twang::bal.table(hold.weights)$unw$std.eff.sz))
	  balance.naive.max = max(abs(twang::bal.table(hold.weights)$unw$std.eff.sz))
	  balance.fully.adjusted.mean = mean(abs(twang::bal.table(hold.weights)[[2]]$std.eff.sz))
	  balance.fully.adjusted.max = max(abs(twang::bal.table(hold.weights)[[2]]$std.eff.sz))
	  
	  if(type == "removal")	{
		balance.ex.mean = vector(length = ncovariates)
		balance.ex.max = vector(length = ncovariates)
		lambda.each.ex.np = vector(length = ncovariates)
		B.each.ex.np = vector(length =  ncovariates)
		delta.each.ex.np = vector(length =  ncovariates)
	  	for(j in 1:ncovariates) {
			#single confounder removal
			names.covariates.remove = names.covariates[-j]
			if(length(names.covariates.remove) == 1) {names.covariate.remove.text = names.covariates.remove}
			if(length(names.covariates.remove) > 1) {names.covariate.remove.text = paste(names.covariates.remove, collapse="+")}
			p.model = glm(as.formula(paste("treatment~", names.covariate.remove.text)),
   			data=data.use, family = binomial)
			pred = predict(p.model, type = "response")
			weights.np = data.use$treatment*(1/pred)+(1-treatment)*(1/(1-pred))
			delta.ps.vec.np.temp = weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - weighted.mean(outcome[treatment==0], weights.np[treatment== 0])

			delta.each.ex.np[j] = delta.ps.vec.np.temp
			lambda.each.ex.np[j] = delta.ps-delta.ps.vec.np.temp
			if(length(names.covariates.remove) > 1) { dat.hold = as.data.frame(cbind(treatment, data.use[,  names(data.use) %in% names.covariates.remove])); hold.weights = twang::dx.wts(x=weights.np, data = dat.hold, vars = names.covariates.remove, treat.var = "treatment", x.as.weights = TRUE, estimand = "ATE")	}
			if(length(names.covariates.remove) == 1) { single.var = data.use[,names(data.use) == names.covariates.remove]; dat.hold = data.frame(treatment = treatment, single.var=single.var); hold.weights = twang::dx.wts(x=weights.np, data = dat.hold, vars = "single.var", treat.var = "treatment", x.as.weights = TRUE, estimand = "ATE")	}
			balance.ex.mean[j] =mean(abs(twang::bal.table(hold.weights)[[2]]$std.eff.sz))
			balance.ex.max[j] = max(abs(twang::bal.table(hold.weights)[[2]]$std.eff.sz)) 
	
		}
		B.sum = sum(abs(lambda.each.ex.np))
		B.each.ex.np = abs(lambda.each.ex.np)/B.sum
		lambda.each.ex.np = as.data.frame(t(lambda.each.ex.np))
		delta.each.ex.np = as.data.frame(t(delta.each.ex.np))
		names(lambda.each.ex.np) = names.covariates
		row.names(lambda.each.ex.np) = "lambda"
		B.each.ex.np = as.data.frame(t(B.each.ex.np))
		names(B.each.ex.np) = names.covariates
		row.names(B.each.ex.np) = "B"
		balance.ex.mean = as.data.frame(t(balance.ex.mean))
		names(balance.ex.mean) = names.covariates
		row.names(balance.ex.mean) = "balance.mean"
		balance.ex.max = as.data.frame(t(balance.ex.max))
		names(balance.ex.max) = names.covariates
		row.names(balance.ex.max) = "balance.max"
		names(delta.each.ex.np) = names.covariates
		row.names(delta.each.ex.np) = "delta"
		results = list("delta.naive" = delta.naive, "p.value.delta.naive" = p.value.delta.naive, "conf.int.delta.naive" = ci.delta.naive, "delta.fully.adjusted" = delta.ps, "p.value.delta.fully.adjusted" = p.value.delta.fully.adjusted, "conf.int.delta.fully.adjusted" = ci.delta.fully.adjusted, "B" = B.each.ex.np)
		if(!Bonly) {
			results = c(results,list("estimated.selection.bias" = lambda.all, "lambda" = lambda.each.ex.np, "delta.each" = delta.each.ex.np))
		}
		if(balance){
			results = c(results,list("balance.naive.mean" = balance.naive.mean, 
	  "balance.naive.max" = balance.naive.max, "balance.fully.adjusted.mean" = balance.fully.adjusted.mean,  "balance.fully.adjusted.max" = balance.fully.adjusted.max, "balance.mean" = balance.ex.mean, "balance.max" = balance.ex.max))
		}
		}
   	
   	if(type == "inclusion")	{
		balance.in.mean = vector(length =  ncovariates)
   		balance.in.max = vector(length =  ncovariates)
   		lambda.each.in.np = vector(length = ncovariates)
   		B.each.in.np = vector(length =  ncovariates)
		delta.each.in.np = vector(length =  ncovariates)
	  	for(j in 1:ncovariates) {
			#single confounder removal
			names.covariates.one = names.covariates[j]
			single.var = data.use[,names(data.use) == names.covariates.one]
			p.model = glm(as.formula(paste("treatment~", names.covariates.one)), data=data.use, family = binomial)
			pred = predict(p.model, type = "response")
			weights.np = treatment*(1/pred)+(1-treatment)*(1/(1-pred))
			delta.ps.vec.np.temp = weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - weighted.mean(outcome[treatment==0], weights.np[treatment== 0])
			delta.each.in.np[j] = delta.ps.vec.np.temp
			lambda.each.in.np[j] = delta.naive-delta.ps.vec.np.temp
			dat.hold = data.frame(treatment = treatment, single.var=single.var)
			hold.weights = twang::dx.wts(x=weights.np, data = dat.hold, vars = "single.var", treat.var = "treatment", x.as.weights = TRUE, estimand = "ATE")
			balance.in.mean[j] = mean(abs(twang::bal.table(hold.weights)[[2]]$std.eff.sz))
			balance.in.max[j] = max(abs(twang::bal.table(hold.weights)[[2]]$std.eff.sz))
		}
		B.sum = sum(abs(lambda.each.in.np))
		B.each.in.np = abs(lambda.each.in.np)/B.sum 
		lambda.each.in.np = as.data.frame(t(lambda.each.in.np))
		delta.each.in.np = as.data.frame(t(delta.each.in.np))
		names(lambda.each.in.np) = names.covariates
		row.names(lambda.each.in.np) = "lambda"
		B.each.in.np = as.data.frame(t(B.each.in.np))
		names(B.each.in.np) = names.covariates
		row.names(B.each.in.np) = "B"
		balance.in.mean = as.data.frame(t(balance.in.mean))
		names(balance.in.mean) = names.covariates
		row.names(balance.in.mean) = "balance.mean"
		balance.in.max = as.data.frame(t(balance.in.max))
		names(balance.in.max) = names.covariates
		row.names(balance.in.max) = "balance.max"
		names(delta.each.in.np) = names.covariates
		row.names(delta.each.in.np) = "delta"
		
		results = list("delta.naive" = delta.naive, "p.value.delta.naive" = p.value.delta.naive, "conf.int.delta.naive" = ci.delta.naive, "delta.fully.adjusted" = delta.ps, "p.value.delta.fully.adjusted" = p.value.delta.fully.adjusted, "conf.int.delta.fully.adjusted" = ci.delta.fully.adjusted, "B" = B.each.in.np)
		if(!Bonly) {
			results = c(results,list( "estimated.selection.bias" = lambda.all, "lambda" = lambda.each.in.np, "delta.each" = delta.each.in.np))
		}
		if(balance){
			results = c(results,list("balance.naive.mean" = balance.naive.mean, 
	  "balance.naive.max" = balance.naive.max, "balance.fully.adjusted.mean" = balance.fully.adjusted.mean,  "balance.fully.adjusted.max" = balance.fully.adjusted.max, "balance.mean" = balance.in.mean, "balance.max" = balance.in.max))
		}
	}
	}
	
	if(standard.error) {
		lambda.each.boot = matrix(nrow = boot.rep, ncol = ncovariates)
		B.each.boot = matrix(nrow = boot.rep, ncol = ncovariates)
		original.outcome = outcome
		original.treatment = treatment
		original.confounders = confounders
		message("Note: bootstrap procedure for standard error estimates is running; this procedure may be very time intensive (on the scale of several minutes to a full day depending on the size of your dataset).")
		progress.ind = round(boot.rep/10*c(1,2,3,4,5,6,7,8,9))
		for(mmm in 1:boot.rep) {
			index.want = sample(1:length(outcome), length(outcome), replace = T)
			outcome = original.outcome[index.want]
			treatment = original.treatment[index.want]
			confounders = original.confounders[index.want,]
			data.use = as.data.frame(cbind(outcome, treatment, confounders))
			delta.naive.vec = mean(outcome[treatment==1]) - mean(outcome[treatment==0])
			if(estimation == "nonparametric") {
				ps.all=twang::ps(as.formula(paste("treatment ~", paste(names.covariates, collapse="+"))), data=data.use,
      			n.trees = n.trees, interaction.depth = interaction.depth, verbose = verbose, estimand = "ATE",
      			stop.method = stop.method, shrinkage = shrinkage, cv.folds=cv.folds)
      			weights.np = twang::get.weights(ps.all,stop.method=stop.method)
      			delta.ps =  weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - 													weighted.mean(outcome[treatment==0], weights.np[treatment== 0])
      			lambda.all = delta.naive.vec-delta.ps
	  	  
	  			if(type == "removal")	{
	  				for(j in 1:ncovariates) {
						#single confounder removal
						names.covariates.remove = names.covariates[-j]
						if(length(names.covariates.remove) > 1) {
						ps.all=twang::ps(as.formula(paste("treatment ~",  paste(names.covariates.remove, collapse="+"))),
   						data=data.use,n.trees = n.trees, interaction.depth = interaction.depth, verbose = verbose, estimand 							= "ATE", stop.method = stop.method, shrinkage = shrinkage, cv.folds=cv.folds)
						weights.np = twang::get.weights(ps.all,stop.method=stop.method)
						delta.ps.vec.np.temp = weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - 										weighted.mean(outcome[treatment==0], weights.np[treatment== 0])
			
						lambda.each.boot[mmm,j] = delta.ps-delta.ps.vec.np.temp	
						}
						if(length(names.covariates.remove) == 1) {
						single.var = data.use[,names(data.use) == names.covariates.remove]
						pred = sapply(single.var,pred.smooth,zz=single.var, y1=treatment)
						weights.np = treatment*(1/pred)+(1-treatment)*(1/(1-pred))
						delta.ps.vec.np.temp  = weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - 										weighted.mean(outcome[treatment==0], weights.np[treatment== 0])

						lambda.each.boot[mmm,j] = delta.ps-delta.ps.vec.np.temp
						}
						}
					B.sum = sum(abs(lambda.each.boot[mmm,]))
					B.each.boot[mmm,] = abs(lambda.each.boot[mmm,])/B.sum
				}
   	
   				if(type == "inclusion")	{
			  		for(j in 1:ncovariates) {
					#single confounder removal
					names.covariates.one = names.covariates[j]
					single.var = data.use[,names(data.use) == names.covariates.one]
					pred = sapply(single.var,pred.smooth,zz=single.var, y1=treatment)
					weights.np = treatment*(1/pred)+(1-treatment)*(1/(1-pred))
					delta.ps.vec.np.temp  = weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - 										weighted.mean(outcome[treatment==0], weights.np[treatment== 0])

					lambda.each.boot[mmm,j] = delta.naive.vec-delta.ps.vec.np.temp
					}
					B.sum = sum(abs(lambda.each.boot[mmm,]))
					B.each.boot[mmm,] = abs(lambda.each.boot[mmm,])/B.sum	
				}
			} #end of nonparametric	
			if(estimation == "parametric") {
				p.model = glm(as.formula(paste("treatment ~", paste(names.covariates, collapse="+"))),
      			data=data.use, family = binomial)
      			pred = predict(p.model, type = "response")
      			weights.np = treatment*(1/pred)+(1-treatment)*(1/(1-pred))

      			delta.ps =  weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - 													weighted.mean(outcome[treatment==0], weights.np[treatment== 0])
      			lambda.all = delta.naive.vec-delta.ps
	  	  
	  			if(type == "removal")	{
	  				for(j in 1:ncovariates) {
						#single confounder removal
						names.covariates.remove = names.covariates[-j]
						if(length(names.covariates.remove) == 1) {names.covariate.remove.text = names.covariates.remove}
						if(length(names.covariates.remove) > 1) {names.covariate.remove.text = paste(names.covariates.remove, collapse="+")}
						p.model = glm(as.formula(paste("treatment~", names.covariate.remove.text)),
   						data=data.use, family = binomial)
						pred = predict(p.model, type = "response")
						weights.np = data.use$treatment*(1/pred)+(1-treatment)*(1/(1-pred))
						delta.ps.vec.np.temp = weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - 										weighted.mean(outcome[treatment==0], weights.np[treatment== 0])
			
						lambda.each.boot[mmm,j] = delta.ps-delta.ps.vec.np.temp	
						}
					B.sum = sum(abs(lambda.each.boot[mmm,]))
					B.each.boot[mmm,] = abs(lambda.each.boot[mmm,])/B.sum
				}
   	
   				if(type == "inclusion")	{
			  		for(j in 1:ncovariates) {
					#single confounder removal
					names.covariates.one = names.covariates[j]
					single.var = data.use[,names(data.use) == names.covariates.one]
					p.model = glm(as.formula(paste("treatment~", names.covariates.one)), data=data.use, family = binomial)
					pred = predict(p.model, type = "response")
					weights.np = treatment*(1/pred)+(1-treatment)*(1/(1-pred))
					delta.ps.vec.np.temp  = weighted.mean(outcome[treatment==1], weights.np[treatment == 1]) - 										weighted.mean(outcome[treatment==0], weights.np[treatment== 0])

					lambda.each.boot[mmm,j] = delta.naive.vec-delta.ps.vec.np.temp
					}
					B.sum = sum(abs(lambda.each.boot[mmm,]))
					B.each.boot[mmm,] = abs(lambda.each.boot[mmm,])/B.sum	
				}
			} #end of parametric	
		if(mmm == progress.ind[1]) {message("Bootstrapping is 10% complete.")}
		if(mmm == progress.ind[2]) {message("Bootstrapping is 20% complete.")}
		if(mmm == progress.ind[3]) {message("Bootstrapping is 30% complete.")}
		if(mmm == progress.ind[4]) {message("Bootstrapping is 40% complete.")}
		if(mmm == progress.ind[5]) {message("Bootstrapping is 50% complete.")}
		if(mmm == progress.ind[6]) {message("Bootstrapping is 60% complete.")}
		if(mmm == progress.ind[7]) {message("Bootstrapping is 70% complete.")}
		if(mmm == progress.ind[8]) {message("Bootstrapping is 80% complete.")}
		if(mmm == progress.ind[9]) {message("Bootstrapping is 90% complete.")}
		} #end of boot rep 
		B.standard.error = apply(B.each.boot,2, sd)
		lambda.standard.error = apply(lambda.each.boot, 2, sd)
		names(B.standard.error) = names.covariates
		names(lambda.standard.error) = names.covariates
		if(Bonly) {results = c(results,list("B.standard.error" = B.standard.error))}
	  	if(!Bonly) {results = c(results,list("B.standard.error" = B.standard.error, 
	  "lambda.standard.error" = lambda.standard.error))}
		}   #end standard error
	
	return(results)
	}
	
bar.sbdecomp = function(output.list, main = ""){
	max.y=max(output.list$B)
	barplot(as.numeric(output.list$B), names.arg = names(output.list$B),
ylab = "Proportion Explained", las=2, cex.names = 0.7, ylim = c(0,max.y+0.05), main = main)
}

