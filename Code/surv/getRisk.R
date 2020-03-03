getRisk <- function(i, surv.obj, exp.obj) {
  rownames(exp.obj) <- make.names(rownames(exp.obj))
   ##remove left-out sample:
     Xlearn <- exp.obj[, -i]
     Ylearn <- surv.obj[-i]
     ##step 1: feature selection, take top 200 cox scores:
       model.step1 = plusMinus(X=t(Xlearn), y=Ylearn, tuningpar="nfeatures", lambda=200)
       keep.genes <- names(model.step1@coefficients)[ abs(model.step1@coefficients) > 0 ]
       Xlearn.200genes <- Xlearn[keep.genes, ]
       exp.200genes <- exp.obj[keep.genes, ]
       ##step 2: principal components regression with PC1-5:
         sig.pc <- prcomp(t(Xlearn.200genes),scale.=TRUE)
         sig.coxph <- coxph(Ylearn~sig.pc$x[,1:5])
         ##final step: extract coefficients:
           sig.coef <- drop((sig.pc$rotation[,1:5] /sig.pc$scale)%*% sig.coxph$coefficients)
           model <- new("ModelLinear", coefficients=sig.coef, modeltype="plusminus")
           ##Make predictions:
             ret = predict(model, newdata=t(exp.200genes), type="lp")
             return( ret@lp[i] )
             }