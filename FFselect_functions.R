
svdwrapper <- function(x) {
    gotit <- F
    try({svdx = svd(x); gotit = TRUE}, silent = TRUE)
    if(gotit) return(svdx)
    try({svdtx = svd(t(x)); gotit = TRUE}, silent = TRUE)
    if(!gotit) stop("svd(x) and svd(t(x)) both failed")
    warning("svd(x) failed but svd(t(x)) worked")
    temp <- svdtx$u
    svdtx$u <- svdtx$v
    svdtx$v <- temp
    svdtx
}

matmult<- function(mat) mat - rep(1, nrow(mat)) %*% t(colMeans(mat))
relmat_centered<-function(data){
    a<-matmult(data)
    b<-tcrossprod(a)
    c<-as.matrix(b/ncol(data))
    rownames(c)<-rownames(data)
    colnames(c)<-rownames(data)
    K_normalization_factor<-1/(sum(diag(c))/nrow(c))
    c<-c*K_normalization_factor
    return(c)
}

MME<-function(X,Z,GI,RI,y) {
    XX<-t(X) %*% RI %*% X
    XZ<-t(X) %*% RI %*% Z
    ZZ<-(t(Z) %*% RI %*% Z) + GI
    Xy<-t(X) %*% RI %*% y
    Zy<-t(Z) %*% RI %*% y
    LHS<-rbind(cbind(XX,XZ), cbind(t(XZ),ZZ))
    RHS<-rbind(Xy,Zy)
    C<-ginv(LHS)
    SOL<-C %*% RHS
    SSR<-t(SOL) %*% RHS
    SOLNS<-cbind(SOL,sqrt(diag(C)))
    return(list(LHS=LHS,RHS=RHS,C=C,SSR=SSR,SOLNS=SOLNS))
}


FFselect<-function(y=NULL, M=NULL, covariate.mat=NULL, env=NULL, allSnp.weight=0.05, max.features=50, adjacent.snps=3, LD.window=500000, ncores=1){
    test.env <- new.env()
    if(is.null(covariate.mat)){covariate.mat<-matrix(1, length(y), 1)}
    y<-as.vector(scale(y))
    test.env$covariate.mat<-covariate.mat
    
    K.snps<-relmat_centered(M[,1:ncol(M)])
    dcmp<-svd(K.snps)
    G.allSnps<-dcmp$u[,1:min(1000,ncol(dcmp$u))]%*%diag(sqrt(dcmp$d[1:min(1000,ncol(dcmp$u))]))
    normalization_factor<-1/(sqrt(sum(G.allSnps^2)/nrow(G.allSnps)))
    G.allSnps<-G.allSnps*normalization_factor
    test.env$G.allSnps<-G.allSnps
    rm(K.snps, dcmp)
    invisible(gc())
    
    if(!is.null(env)){
        env<-as.factor(as.numeric(env))
        G0<-model.matrix(~env-1)
        test.env$G0<-G0
        
        
        #adjust y for environmental effects (model includes polygenic effect to prevent overadjustment in scenario of environment confounded with genetic similarity)
        
        print("Estimating shared environment effects")
        max_env_weight<-function(x){
            fit <- lrgpr(formula=if(ncol(covariate.mat)>1){as.formula(y~covariate.mat[,-1])}else{as.formula(y~1)}, svdwrapper(cbind(G0*sqrt(x), G.allSnps*sqrt(1-x))), scale=FALSE)
            #print(paste("env.weight=", round(x, digits=3), " LL=", round(fit$logLik,digits=3), sep=""))
            return(c(fit$logLik))
        }
        
        env.weight<-optimize(max_env_weight, lower=0.001, upper=1, tol = 0.001, maximum=TRUE)$maximum
        test.env$env.weight<-env.weight
        
        fit <- lrgpr(formula=if(ncol(covariate.mat)>1){as.formula(y~covariate.mat[,-1])}else{as.formula(y~1)}, svdwrapper(cbind(G0*sqrt(env.weight), G.allSnps*sqrt(1-env.weight))), scale=FALSE)
        GI<-diag(c(rep(fit$sigSq_e/(fit$sigSq_a*env.weight),ncol(G0)),rep( fit$sigSq_e/(fit$sigSq_a*(1-env.weight)),ncol(G.allSnps))))
        Z<-cbind(G0, G.allSnps)
        RI<-diag(length(y))
        
        res<-MME(X=covariate.mat, Z=Z, GI=GI, RI=RI, y=y)
        env.intercept<-G0%*%res$SOLNS[(ncol(covariate.mat)+1):(ncol(G0)+(ncol(covariate.mat))),1]
        y2<-as.vector(scale(y-env.intercept))
        test.env$y2<-y2
    }else{
        G0<-matrix(1,length(y),1)
        env.weight<-0
        y2<-y
        test.env$G0<-G0
        test.env$env.weight<-env.weight
        test.env$y2<-y2
    }
    
    #obtain subset of top ranked SNPs for feature selection
    pg<-lrgpr(formula=if(ncol(covariate.mat)>1){as.formula(y2~covariate.mat[,-1])}else{as.formula(y2~1)}, svdwrapper(G.allSnps), scale=FALSE, diagnostic=TRUE)
    #ResultA<-glmApply(formula=if(ncol(covariate.mat)>1){as.formula(pg$residuals~covariate.mat[,-1] + SNP)}else{as.formula(pg$residuals~SNP)}, features=M, nthreads=ncores, progress=FALSE)
    ResultA<-glmApply(pg$residuals~SNP, features=M, nthreads=ncores, progress=FALSE)
    Result<-ResultA$pValues
    chisq <- qchisq(Result,1,lower.tail=FALSE)
    chisq[is.na(chisq)]<-0
    lambda <- median(chisq)/qchisq(0.5,1)
    newchisq <- chisq/lambda
    Result <- pchisq(newchisq, df=1,lower.tail=FALSE)
    names(Result)<-colnames(M)
    top_snps<-names(Result[Result<0.05])
    
    discoveries<-vector()
    counter<-1
    bic<-vector()
    g<-pg$sigSq_a
    dcmp<-svd(G.allSnps)
    Result <- lrgprApply(formula=if(ncol(covariate.mat)>1){as.formula(y2~covariate.mat[,-1] + SNP)}else{as.formula(y2~SNP)}, features=M[,!colnames(M)%in%discoveries & colnames(M)%in%top_snps], decomp=dcmp, reEstimateDelta=FALSE, scale=FALSE, nthreads=ncores, progress=FALSE)
    bonferonni.correction<--log10(0.05/ncol(M))
    
    print("Selecting feature SNPs")
    
    print(paste("step=0", " LL=",round(pg$logLik, digits=3), " BIC=",round(pg$BIC, digits=3), " geneticVar.goal=", round(pg$sigSq_a*allSnp.weight, digits=3), " geneticVar=",round(pg$sigSq_a, digits=3), " residualVar=", round(pg$sigSq_e, digits=3), sep=""))
    
    while((counter<=min(max.features, length(top_snps)) & g>(pg$sigSq_a*allSnp.weight)) | (-log10(min(Result))>bonferonni.correction) ){
        discoveries<-c(discoveries, names(Result)[which.min(Result)])
        
        snp.index<-match(discoveries, colnames(M))
        if(adjacent.snps==0){
            discoveriesPlusAdjacent<-discoveries
        }else{
            xx<-t(sapply(snp.index, function(x) seq((x-adjacent.snps),(x+adjacent.snps))))
            xx[xx<1]<-1
            xx[xx>length(colnames(M))]<-length(colnames(M))
            check.same.chr<-t(apply(xx,1, function(x) ldply(strsplit(colnames(M)[x], split="_"))[[2]]))
            snp.index.PlusAdjacent<-unique(xx[cbind(check.same.chr[,1]==check.same.chr[,1],check.same.chr[,2]==check.same.chr[,1],check.same.chr[,3]==check.same.chr[,1])])
            discoveriesPlusAdjacent<-colnames(M)[snp.index.PlusAdjacent]
        }
        xx<-as.matrix(if(ncol(covariate.mat)>1){cbind(covariate.mat[,-1], M[,colnames(M)%in%discoveriesPlusAdjacent])}else{M[,colnames(M)%in%discoveriesPlusAdjacent]})
        xx<-xx[,!is.na(lm(y2~xx)$coefficients[-1])]
        
        pvals<-lrgpr(y2~xx, decomp=dcmp,  scale=FALSE, diagnostic=TRUE)
        x<-xx%*%as.matrix(pvals$coefficients[-1])
        Result <- lrgprApply(y2~x+SNP, features=M[,!colnames(M)%in%discoveriesPlusAdjacent & colnames(M)%in%top_snps], decomp=dcmp, reEstimateDelta=FALSE, scale=FALSE, nthreads=ncores, progress=FALSE, delta=pvals$delta)
        
        g<-pvals$sigSq_a
        bic[counter]<-pvals$BIC
        print(paste("step=",counter, " LL=",round(pvals$logLik, digits=3), " BIC=",round(pvals$BIC, digits=3), " geneticVar.goal=", round(pg$sigSq_a*allSnp.weight, digits=3), " geneticVar=",round(pvals$sigSq_a, digits=3), " residualVar=", round(pvals$sigSq_e, digits=3), sep=""))
        counter<-counter+1
    }
    
    bic.min<-which.min(bic)
    if(bic.min!=length(discoveries)){
        discoveries<-discoveries[1:bic.min]
        xx<-t(sapply(snp.index, function(x) seq((x-adjacent.snps),(x+adjacent.snps))))
        xx[xx<1]<-1
        xx[xx>length(colnames(M))]<-length(colnames(M))
        check.same.chr<-t(apply(xx,1, function(x) ldply(strsplit(colnames(M)[x], split="_"))[[2]]))
        snp.index.PlusAdjacent<-unique(xx[cbind(check.same.chr[,1]==check.same.chr[,1],check.same.chr[,2]==check.same.chr[,1],check.same.chr[,3]==check.same.chr[,1])])
        discoveriesPlusAdjacent<-colnames(M)[snp.index.PlusAdjacent]
    }
    
    G1<-scale(M[,discoveriesPlusAdjacent], scale=FALSE)/sqrt(length(discoveriesPlusAdjacent))
    normalization_factor<-1/(sqrt(sum(G1^2)/nrow(G1)))
    G1<-G1*normalization_factor
    test.env$G1<-G1
    
    print("Optimizing weight1")
    max_G1_weight<-function(x){
        fit <- lrgpr(formula=if(ncol(covariate.mat)>1){as.formula(y2~covariate.mat[,-1])}else{as.formula(y2~1)}, svdwrapper(cbind(G1*sqrt(x), G.allSnps*sqrt(1-x))), scale=FALSE)
        print(paste("weight1=", round(x, digits=3), " LL=", round(fit$logLik,digits=3), sep=""))
        return(c(fit$logLik))
    }
    
    selectSnp.weight<-optimize(max_G1_weight, lower=0.001, upper=1, tol = 0.001, maximum=TRUE)$maximum
    test.env$selectSnp.weight<-selectSnp.weight
    
    if(!is.null(env)){
        print("Optimizing weight2")
        max_env_weight2<-function(x){
            fit <- lrgpr(formula=if(ncol(covariate.mat)>1){as.formula(y~covariate.mat[,-1])}else{as.formula(y~1)}, svdwrapper(cbind(G0*sqrt(x), (cbind(G1*sqrt(selectSnp.weight), G.allSnps*sqrt(1-selectSnp.weight)))*sqrt(1-x))), scale=FALSE)
            print(paste("weight2=", round(x, digits=3), " LL=", round(fit$logLik,digits=3), sep=""))
            return(c(fit$logLik))
        }
        
        env.weight2<-optimize(max_env_weight2, lower=max(0.001, env.weight-0.1), upper=min(1, env.weight+0.1), tol = 0.001, maximum=TRUE)$maximum
        test.env$env.weight2<-env.weight2
    }else{
        env.weight2<-0
        test.env$env.weight2<-env.weight2
    }
    
    print("Beginning GWAS")
    
    datLD<-matrix(NA, ncol(M), length(discoveries))
    rownames(datLD)<-colnames(M)
    colnames(datLD)<-discoveries
    row_chr<-ldply(strsplit(as.character(rownames(datLD)),split= "_"))[[2]]
    row_bp<-as.numeric(ldply(strsplit(as.character(rownames(datLD)),split= "_"))[[3]])
    col_chr<-ldply(strsplit(discoveries,split= "_"))[[2]]
    col_bp<-as.numeric(ldply(strsplit(discoveries,split= "_"))[[3]])
    for (i in 1:ncol(datLD)){
        datLD[,i]<-!(row_chr==col_chr[i] & row_bp>(col_bp[i]-LD.window)  & row_bp<(col_bp[i]+LD.window))
    }
    yy<-datLD
    keep_snps<-apply(yy, 1, function (xy) paste(colnames(datLD)[xy], collapse="*"))
    rm(yy, datLD)
    unique_keep_snps<-unique(keep_snps)
    keep_snp_index<-list()
    for (i in 1:length(unique_keep_snps)){keep_snp_index[[i]]<-which(keep_snps%in%unique_keep_snps[i])}
    
    Result<-array(1, ncol(M))
    names(Result)<-colnames(M)
    
    for (jj in 1:length(unique_keep_snps)){
        print(paste("GWAS subset" ,jj, "of", length(unique_keep_snps), "      Number of test snps=", length(keep_snp_index[[jj]])))
        snp.index<-match(strsplit(unique_keep_snps[jj], "*", fixed=TRUE)[[1]], colnames(M))
        if(adjacent.snps==0){
            snp.index.PlusAdjacent<-snp.index
        }else{
            xx<-t(sapply(snp.index, function(x) seq((x-adjacent.snps),(x+adjacent.snps))))
            xx[xx<1]<-1
            xx[xx>length(colnames(M))]<-length(colnames(M))
            check.same.chr<-t(apply(xx,1, function(x) ldply(strsplit(colnames(M)[x], split="_"))[[2]]))
            snp.index.PlusAdjacent<-unique(xx[cbind(check.same.chr[,1]==check.same.chr[,1],check.same.chr[,2]==check.same.chr[,1],check.same.chr[,3]==check.same.chr[,1])])
        }
        M2<-scale(M[,snp.index.PlusAdjacent], scale=FALSE)
        G1<-M2/sqrt(ncol(M2))
        normalization_factor<-1/(sqrt(sum(G1^2)/nrow(G1)))
        G1<-G1*normalization_factor
        
        fit <- lrgprApply(formula=if(ncol(covariate.mat)>1){as.formula(y~covariate.mat[,-1] + SNP)}else{as.formula(y~SNP)}, features=M[,keep_snp_index[[jj]]], svdwrapper(cbind(G0*sqrt(env.weight2), (cbind(G1*sqrt(selectSnp.weight), G.allSnps*sqrt(1-selectSnp.weight)))*sqrt(1-env.weight2))),  reEstimateDelta=FALSE, scale=FALSE, nthreads=ncores, progress=FALSE)
        
        Result[match(names(fit), names(Result))]<-fit
    }
    return(list(pvals=Result, featureSNP=discoveries, RM1weight=env.weight2, RM2weight=selectSnp.weight*(1-env.weight2) , RM3weight= (1-selectSnp.weight)*(1-env.weight2)))
}
