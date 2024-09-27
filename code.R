##Codes for running lasso-cox regression
# 清除环境中的所有对象
rm(list = ls())

# 加载所需的库
library(glmnet)
library(readxl)
library(openxlsx)
# 读取数据
lasso <- read_xlsx("train.xlsx")
fitx <- as.matrix(lasso[, 4:1144])
fity <- as.matrix(lasso[, 2:3])

# 训练 Lasso-Cox 回归模型
fit <- glmnet(fitx, fity, family = "cox")
plot(fit, xvar = "lambda", label = TRUE)

# 使用交叉验证选择最佳的 lambda 值
set.seed(1111)
cvfit <- cv.glmnet(fitx, fity, family = "cox", type.measure = "deviance", nfolds = 10)
plot(cvfit)

# 提取最佳 lambda 对应的非零特征的系数
selected_coef <- coef(cvfit, s = cvfit$lambda.min)
non_zero_indices <- which(selected_coef != 0)

# 提取最终选择的特征名称
selected_features <- colnames(fitx)[non_zero_indices]

# 将选定特征的数据保存到 Excel 文件中
write.xlsx(selected_data, "selected_features.xlsx", row.Names = FALSE)


##Codes for cut-point of PSHCC determination
cutpoint <- read_xlsx("rad_features_train.xlsx")
library(survminer)
cut_point <- surv_cutpoint(data = cutpoint, time = "time", event = "status", variables = "sig",
                           minprop = 0.1, progressbar = TRUE)
#plot(cut_point, palette = c("#CD3333", "#1874CD"), legend = "")
plot(cut_point, palette = c("#EA8379", "#7DAEE0"), legend = "")
cut_point
#---------------------------------------------------------------------------------------------------------------------------------------------------

##Codes for generating survival curves
training_data <- read_xlsx("rad_features_train.xlsx")
validation_data <- read_xlsx("rad_features_test.xlsx")
library(survminer)
library(survival)
Signature_OS_training <- survfit(Surv(OS,OSstatus)~Sub,data=training_data)
ggsurvplot(Signature_OS_training,training_data,size=1,
           palette = "nejm",pval = TRUE,risk.table = TRUE,risk.table.col = "strata",
           break.x.by=10,break.y.by=0.2,xlim=c(0,72),axes.offset=T)


Signature_OS_validation <- survfit(Surv(OS,OSstatus)~Sub,data=validation_data)
ggsurvplot(Signature_OS_validation,validation_data,size=1,
           palette = "nejm",pval = TRUE,risk.table = TRUE,risk.table.col = "strata",
           break.x.by=10,break.y.by=0.2,xlim=c(0,72),axes.offset=T)
#-------------------------------------------------------------------------------------------------------------------------------------------------

##Codes for development and validation of nomogram
training_data <- read_xlsx("train.xlsx",sheet="nomo")
validation_data <- read_xlsx("test.xlsx",sheet="nomo")
library(rms)
library(survival)
dd_training <- datadist(training_data)
options(datadist="dd_training")

#Codes for constructing and plotting nomogram
f_OS <- cph(Surv(OS,OSstatus) ~ Sig_p+Sig_r, data = training_data, x = T, y = T, surv = T, time.inc = 60)
surv_OS <- Survival(f_OS)

nom_OS <- nomogram(f_OS, fun = list(function(x) surv_OS(24, x),
                                    function(x) surv_OS(36, x),
                                    function(x) surv_OS(60, x)),
                   fun.at = c(seq(.1,.9,by = .1),.95),
                   funlabel = c("2-year OS","3-year  OS","5-year  OS"),lp = 0)
plot(nom_OS)

png("nomogram_plot.png", width = 1200, height = 600) # 调整这里的 width 和 height 来改变图形尺寸

# 创建和绘制 nomogram
nom_OS <- nomogram(f_OS, fun = list(function(x) surv_OS(24, x),
                                    function(x) surv_OS(36, x),
                                    function(x) surv_OS(60, x)),
                   fun.at = c(seq(.1, .9, by = .1), .95),
                   funlabel = c("2-year OS", "3-year OS", "5-year OS"), lp = 0)
plot(nom_OS)

# 关闭图形设备，保存图像
dev.off()

#Codes for calculating C-index in training cohort
library("boot")
rcorrcens(Surv(OS,OSstatus) ~ predict(f_OS), data = training_data)
rcorrcens(Surv(DFS,DFSstatus) ~ predict(f_DFS), data = training_data)
#Codes for calculating C-index in validation cohort
f_OS_val <- cph(Surv(OS,OSstatus)~predict(f_OS, newdata = validation_data), x = T, y = T, surv = T,time.inc = 60, data = validation_data)
rcorrcens(Surv(OS,OSstatus) ~ predict(f_OS, newdata = validation_data), data = validation_data)

#Codes for calibration curve in training cohort
f_OS_24 <- cph(Surv(OS,OSstatus) ~ Sig_p+Sig_r, data = training_data, x = T, y = T, surv = T, time.inc = 24)

cal_OS_training_24 <-calibrate(f_OS_24, cmethod = "KM", method = "boot", u = 24, m = 22, B = 1000)
cal_OS_training_36 <-calibrate(f_OS_36, cmethod = "KM", method = "boot", u = 36, m = 22, B = 1000)
cal_OS_training <-calibrate(f_OS, cmethod = "KM", method = "boot", u = 60, m = 22, B = 1000)
plot(cal_OS_training,xlim = c(0,1),ylim = c(0,1), col =  c("#EA8379"), errbar.col =  c("#EA8379"))
plot(cal_OS_training_36,col =  c("#B395BD"), errbar.col =  c("#B395BD") , add = TRUE)
plot(cal_OS_training_24,col =  c("#7DAEE0"), errbar.col =  c("#7DAEE0") , add = TRUE)
legend("bottomright", legend = c("5-year OS ", "3-year OS", "2-year  OS"), col = c("#EA8379", "#B395BD", "#7DAEE0"), pch = 16)

coxfit1 <- cph(Surv(OS,OSstatus) ~ Sig_p+Sig_r,
               data = training_data, x=T,y=T,surv = T,
               time.inc = 60 # 1 年
)
coxfit2 <- cph(Surv(OS,OSstatus) ~ Sig_p+Sig_r,
               data = training_data, x=T,y=T,surv = T,
               time.inc = 24 # 1 年
)
coxfit2 <- cph(Surv(OS,OSstatus) ~ Sig_p+Sig_r,
               data = training_data, x=T,y=T,surv = T,
               time.inc = 36 # 1 年
)
# m=50表示每次计算50个样本，一般取4-6个点，u=365和上面的time.inc对应
cal1 <- calibrate(coxfit1, cmethod="KM", method="boot",u=60,m=22,B=1000)
cal2 <- calibrate(coxfit2, cmethod="KM", method="boot",u=24,m=22,B=1000) 
cal3 <- calibrate(coxfit2, cmethod="KM", method="boot",u=36,m=22,B=1000) 
plot(cal1,
     #lwd = 2, # 误差线粗细
     lty = 0, # 误差线类型，可选0-6
     errbar.col = c("#2166AC"), # 误差线颜色
     xlim = c(0.4,1),ylim= c(0.4,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) # 字体大小
lines(cal1[,c('mean.predicted',"KM")], 
      type = 'b', # 连线的类型，可以是"p","b","o"
      lwd = 3, # 连线的粗细
      pch = 16, # 点的形状，可以是0-20
      col = "tomato") # 连线的颜色
box(lwd = 1) # 边框粗细
abline(0,1,lty = 3, # 对角线为虚线
       lwd = 2, # 对角线的粗细
       col = "grey70" # 对角线的颜色
) 

plot(cal1,lwd = 1,lty = 0,errbar.col = c("#7DAEE0"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#7DAEE0"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#7DAEE0"), pch = 16)

plot(cal2,lwd = 1,lty = 0,errbar.col = c("#EA8379"),
     xlim = c(0,1),ylim= c(0,1),col = c("#EA8379"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#EA8379"), pch = 16)
# 假设cal2是已经通过某种方法（如calibrate函数）得到的校准数据
# 绘制校准曲线，首先使用plot函数
plot(cal3,lwd = 1,lty = 0,errbar.col = c("#B395BD"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B395BD"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#B395BD"), pch = 16)


abline(0,1, lwd = 1, lty = 3, col = c("grey70"))

legend("bottomright", #图例的位置
       legend = c("5-year OS","2-year OS","3-year OS"), #图例文字
       col =c("#7DAEE0","#EA8379","#B395BD"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.0,#图例字体大小
       bty = "n")#不显示图例边框

#Codes for calibration curve in validation cohort
f_OS_val_36 <- cph(Surv(OS,OSstatus) ~ Sig_p+Sig_r, data = validation_data, x = T, y = T, surv = T, time.inc = 36)
f_OS_val_24 <- cph(Surv(OS,OSstatus) ~ Sig_p+Sig_r, data = validation_data, x = T, y = T, surv = T, time.inc = 24)

cal_OS_validation_24 <-calibrate(f_OS_val_24 , cmethod = "KM", method = "boot", u = 24, m = 10, B = 1000)
cal_OS_validation_36 <-calibrate(f_OS_val_36, cmethod = "KM", method = "boot", u = 36, m = 10, B = 1000)
cal_OS_validation <-calibrate(f_OS_val, cmethod = "KM", method = "boot", u = 60, m = 10, B = 1000)
plot(cal_OS_validation,xlim = c(0,1),ylim = c(0,1), col =  c("#EA8379"), errbar.col =  c("#EA8379"))
plot(cal_OS_validation_36,col =  c("#B395BD"), errbar.col =  c("#B395BD") , add = TRUE)
plot(cal_OS_validation_24,col =  c("#7DAEE0"), errbar.col =  c("#7DAEE0") , add = TRUE)
legend("bottomright", legend = c("5-year OS ", "3-year OS", "2-year OS"), col = c("#EA8379", "#B395BD", "#7DAEE0"), pch = 16)


# 1年
coxfit1 <- cph(Surv(OS,OSstatus) ~ Sig_p+Sig_r,
               data = validation_data, x=T,y=T,surv = T,
               time.inc = 60 # 1 年
)
coxfit2 <- cph(Surv(OS,OSstatus) ~ Sig_p+Sig_r,
               data = validation_data, x=T,y=T,surv = T,
               time.inc = 24 # 1 年
)
coxfit2 <- cph(Surv(OS,OSstatus) ~ Sig_p+Sig_r,
               data = validation_data, x=T,y=T,surv = T,
               time.inc = 36 # 1 年
)
# m=50表示每次计算50个样本，一般取4-6个点，u=365和上面的time.inc对应
cal1 <- calibrate(coxfit1, cmethod="KM", method="boot",u=60,m=10,B=1000)
cal2 <- calibrate(coxfit2, cmethod="KM", method="boot",u=24,m=10,B=1000) 
cal3 <- calibrate(coxfit2, cmethod="KM", method="boot",u=36,m=10,B=1000) 
plot(cal1,
     #lwd = 2, # 误差线粗细
     lty = 0, # 误差线类型，可选0-6
     errbar.col = c("#2166AC"), # 误差线颜色
     xlim = c(0.4,1),ylim= c(0.4,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) # 字体大小
lines(cal1[,c('mean.predicted',"KM")], 
      type = 'b', # 连线的类型，可以是"p","b","o"
      lwd = 3, # 连线的粗细
      pch = 16, # 点的形状，可以是0-20
      col = "tomato") # 连线的颜色
box(lwd = 1) # 边框粗细
abline(0,1,lty = 3, # 对角线为虚线
       lwd = 2, # 对角线的粗细
       col = "grey70" # 对角线的颜色
) 

plot(cal1,lwd = 1,lty = 0,errbar.col = c("#7DAEE0"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#7DAEE0"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#7DAEE0"), pch = 16)

plot(cal2,lwd = 1,lty = 0,errbar.col = c("#EA8379"),
     xlim = c(0,1),ylim= c(0,1),col = c("#EA8379"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 3, col = c("#EA8379"), pch = 16)
# 假设cal2是已经通过某种方法（如calibrate函数）得到的校准数据
# 绘制校准曲线，首先使用plot函数
plot(cal2, lwd = 1, lty = 0, errbar.col = "#EA8379", xlim = c(0,1), ylim = c(0,1), col = "#EA8379", add = TRUE)

# 获取mean.predicted和KM的值，然后在y值上加上0.05来进行平移
# 使用lines函数添加平移后的曲线
shift_amount <- 0.005  # 定义平移量
lines(cal2[, 'mean.predicted'], cal2[, 'KM'] + shift_amount, type = 'b', lwd = 3, col = "#EA8379", pch = 16)

plot(cal3,lwd = 1,lty = 0,errbar.col = c("#B395BD"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B395BD"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 2, col = c("#B395BD"), pch = 16)

abline(0,1, lwd = 1, lty = 3, col = c("grey70"))

legend("bottomright", #图例的位置
       legend = c("5-year OS","2-year OS","3-year OS"), #图例文字
       col =c("#7DAEE0","#EA8379","#B395BD"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.0,#图例字体大小
       bty = "n")#不显示图例边框

library(rms)

# 为2年、3年、5年生存率创建Cox模型
f_OS_val_24 <- cph(Surv(OS, OSstatus) ~ Sig_p + Sig_r, data = validation_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 24)
f_OS_val_36 <- cph(Surv(OS, OSstatus) ~ Sig_p + Sig_r, data = validation_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)
f_OS_val_60 <- cph(Surv(OS, OSstatus) ~ Sig_p + Sig_r, data = validation_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 60)

library(rms)

# 假设您已经加载并准备好了数据和模型
# 首先，您需要确保cph和calibrate的正确应用

# 创建Cox模型
f_OS_val_24 <- cph(Surv(OS, OSstatus) ~ Sig_p + Sig_r, data = validation_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 24)
f_OS_val_36 <- cph(Surv(OS, OSstatus) ~ Sig_p + Sig_r, data = validation_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)
f_OS_val_60 <- cph(Surv(OS, OSstatus) ~ Sig_p + Sig_r, data = validation_data, x = TRUE, y = TRUE, surv = TRUE, time.inc = 60)

# 校准模型
cal_OS_validation_24 <- calibrate(f_OS_val_24, cmethod = "KM", method = "boot", u = 24, m = 10, B = 1000)
cal_OS_validation_36 <- calibrate(f_OS_val_36, cmethod = "KM", method = "boot", u = 36, m = 10, B = 1000)
cal_OS_validation_60 <- calibrate(f_OS_val_60, cmethod = "KM", method = "boot", u = 60, m = 10, B = 1000)

# 绘图
plot(cal_OS_validation_60, xlim = c(0, 1), ylim = c(0, 1), col = "#EA8379", errbar.col = "#EA8379", lwd = 2, cex = 1.2, xlab = "Predicted Probability", ylab = "Observed Probability")
# 添加偏移以区分曲线
plot(cal_OS_validation_36, col = "#B395BD", errbar.col = "#B395BD", lwd = 2, cex = 1.2, lty = 2, add = TRUE, xlim = c(0.01, 1.01))
plot(cal_OS_validation_24, col = "#7DAEE0", errbar.col = "#7DAEE0", lwd = 2, cex = 1.2, lty = 3, add = TRUE, xlim = c(-0.01, 0.99))

# 添加图例
legend("bottomright", legend = c("5-year OS", "3-year OS", "2-year OS"), col = c("#EA8379", "#B395BD", "#7DAEE0"), pch = c(19, 17, 15), lty = c(1, 2, 3), lwd = 2)


#Codes for time-independent ROC curve and AUROC comparison in training cohort
library(riskRegression)
library(survival)
Srv_OS <- Surv(training_data$OS, training_data$OSstatus)
coxmod_nomogram_OS <- coxph(Srv_OS ~ Sig_p+Sig_r, data=training_data,x=TRUE)
coxmod_Sig_p_OS <- coxph(Srv_OS ~ Sig_p, data=training_data,x=TRUE)
coxmod_Sig_r_OS <- coxph(Srv_OS ~ Sig_r, data=training_data,x=TRUE)

ROC_training_OS <- Score(list("nomo"=coxmod_nomogram_OS, "path_sig"=coxmod_Sig_p_OS,
                              "rad_sig"=coxmod_Sig_r_OS),
                         formula=Hist(OS,OSstatus)~1,data=training_data,
                         times=60,plots="ROC",summary="risk",contrasts=TRUE)
ROC_training_OS
plotROC(ROC_training_OS, col = c("#EA8379", "#7DAEE0", "#B395BD"), lwd = 1.5, cex.legend = 0.8)

#Codes for time-independent ROC curve and AUROC comparison in validation cohort
ROC_validation_OS <- Score(list("nomo"=coxmod_nomogram_OS, "pat_sig"=coxmod_Sig_p_OS,
                                "rad_sig"=coxmod_Sig_r_OS),
                           formula=Hist(OS,OSstatus)~1,data=validation_data,
                           times=60,plots="ROC",summary="risk",contrasts=TRUE)

ROC_validation_OS
plotROC(ROC_validation_OS, col = c("#EA8379", "#7DAEE0", "#B395BD"), lwd = 1.5)

#Codes for C-index comparison in training cohort
library(compareC)
f_OS_Sig_p <- cph(Surv(OS,OSstatus) ~ Sig_p, data = training_data, x = T, y = T, surv = T, time.inc = 60)
rcorrcens(Surv(OS,OSstatus) ~ predict(f_OS_Sig_p), data = training_data)

f_OS_Sig_r <- cph(Surv(OS,OSstatus) ~ Sig_r, data = training_data, x = T, y = T, surv = T, time.inc = 60)
rcorrcens(Surv(OS,OSstatus) ~ predict(f_OS_Sig_r), data = training_data)

compareC(timeX = training_data$OS, statusX = training_data$OSstatus,
         scoreY = predict(f_OS), scoreZ = predict(f_OS_Sig_p))
compareC(timeX = training_data$OS, statusX = training_data$OSstatus,
         scoreY = predict(f_OS), scoreZ = predict(f_OS_Sig_r))
compareC(timeX = training_data$OS, statusX = training_data$OSstatus,
         scoreY = predict(f_OS_Sig_p), scoreZ = predict(f_OS_Sig_r))

#Codes for C-index comparison in validation cohort
f_OS_Sig_p_val <- cph(Surv(OS,OSstatus) ~ Sig_p, data = validation_data, x = T, y = T, surv = T, time.inc = 60)
rcorrcens(Surv(OS,OSstatus) ~ predict(f_OS_Sig_p_val), data = validation_data)

f_OS_Sig_r_val <- cph(Surv(OS,OSstatus) ~ Sig_r, data = validation_data, x = T, y = T, surv = T, time.inc = 60)
rcorrcens(Surv(OS,OSstatus) ~ predict(f_OS_Sig_r_val), data = validation_data)

compareC(timeX = validation_data$OS, statusX = validation_data$OSstatus,
         scoreY = predict(f_OS_val), scoreZ = predict(f_OS_Sig_p_val))
compareC(timeX = validation_data$OS, statusX = validation_data$OSstatus,
         scoreY = predict(f_OS_val), scoreZ = predict(f_OS_Sig_r_val))
compareC(timeX = validation_data$OS, statusX = validation_data$OSstatus,
         scoreY = predict(f_OS_Sig_p_val), scoreZ = predict(f_OS_Sig_r_val))


#Codes for decision curve analysis in training cohort
source("stdca.R")
library(survival)
Srv_OS <- Surv(training_data$OS, training_data$OSstatus)
coxmod_radiopathomics <- coxph(Srv_OS ~ Sig_p+Sig_r, data=training_data,x=TRUE)
coxmod_pathomics <- coxph(Srv_OS ~ Sig_p, data=training_data,x=TRUE)
coxmod_radiomics <- coxph(Srv_OS ~ Sig_r, data=training_data,x=TRUE)

training_data$radiopathomics <- c(1- (summary(survfit(coxmod_radiopathomics,newdata=training_data), times=60)$surv))
training_data$pathomics <- c(1- (summary(survfit(coxmod_pathomics,newdata=training_data), times=60)$surv))
training_data$radiomics <- c(1- (summary(survfit(coxmod_radiomics,newdata=training_data), times=60)$surv))

# 创建包含所有预测因子的数据框
dca_data <- data.frame(OS = as.numeric(training_data$OS),
                       OSstatus = as.numeric(training_data$OSstatus),
                       radiopathomics = as.numeric(training_data$radiopathomics),
                       pathomics = as.numeric(training_data$pathomics),
                       radiomics = as.numeric(training_data$radiomics))
names(dca_data)

# 使用 stdca() 函数进行决策曲线分析，同时包含三个预测因子
dca_result <- stdca(data = dca_data,
                    outcome = "OSstatus",
                    ttoutcome = "OS",
                    timepoint = 60,
                    predictors = c("radiopathomics","pathomics","radiomics"),
                    cmprsk = FALSE,
                    smooth = TRUE,
                    xstop = 0.9)
dcr_training_Sig_p_r_OS <- stdca(data=dca_data, outcome="OSstatus", ttoutcome="OS", timepoint=60, 
                             predictors="radiopathomics", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
dcr_training_Sig_p_OS <- stdca(data=dca_data, outcome="OSstatus", ttoutcome="OS", timepoint=60, 
                                  predictors="pathomics", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
dcr_training_Sig_r_OS <- stdca(data=dca_data, outcome="OSstatus", ttoutcome="OS", timepoint=60, 
                             predictors="radiomics", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
plot(dcr_training_Sig_p_r_OS$net.benefit$none,xlim=c(0, 90), ylim=c(-0.1, 0.4), 
     xlab="Threshold probability (%)", ylab="Net benefit",type="l", lwd=1,col="#fccccb") 
lines(dcr_training_Sig_p_r_OS$net.benefit$all, type="l", lwd=1,col="#bdb5e1") 
lines(dcr_training_Sig_p_r_OS$net.benefit$radiopathomics, type="l",lwd=1,col="#7ac7e2")
lines(dcr_training_Sig_p_OS$net.benefit$pathomics, type="l",lwd=1,col="#f7df87")
lines(dcr_training_Sig_r_OS$net.benefit$radiomics, type="l",lwd=1,col="#e3716e")
legend("topright", legend = c("None", "All", "nomo", "path_sig", "rad_sig"), 
       col = c("#fccccb", "#bdb5e1", "#7ac7e2", "#f7df87", "#e3716e"), lty = 1, lwd = 1, cex = 0.7)


#Codes for decision curve analysis in validation cohort
validation_data$radiopathomics <- c(1- (summary(survfit(coxmod_radiopathomics,newdata=validation_data), times=60)$surv))
validation_data$pathomics <- c(1- (summary(survfit(coxmod_Sig_p_OS,newdata=validation_data), times=60)$surv))
validation_data$radiomics <- c(1- (summary(survfit(coxmod_Sig_r_OS,newdata=validation_data), times=60)$surv))
# 创建包含所有预测因子的数据框
dca_data_val <- data.frame(OS = as.numeric(validation_data$OS),
                           OSstatus = as.numeric(validation_data$OSstatus),
                           radiopathomics = as.numeric(validation_data$radiopathomics),
                           pathomics = as.numeric(validation_data$pathomics),
                           radiomics = as.numeric(validation_data$radiomics))
colnames(dca_data_val)

dcr_validation_Sig_p_r_OS <- stdca(data=dca_data_val, outcome="OSstatus", ttoutcome="OS", timepoint=60, 
                                  predictors="radiopathomics", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
dcr_validation_Sig_p_OS <- stdca(data=dca_data_val, outcome="OSstatus", ttoutcome="OS", timepoint=60, 
                                 predictors="pathomics", cmprsk=FALSE, smooth=TRUE, xstop=0.9)
dcr_validation_Sig_r_OS <- stdca(data=dca_data_val, outcome="OSstatus", ttoutcome="OS", timepoint=60, 
                               predictors="radiomics", cmprsk=FALSE, smooth=TRUE, xstop=0.9)


plot(dcr_validation_Sig_p_r_OS$net.benefit$none,xlim=c(0, 90), ylim=c(-0.1, 0.3), 
     xlab="Threshold probability (%)", ylab="Net benefit",type="l", lwd=1,col="#fccccb") 
lines(dcr_validation_Sig_p_r_OS$net.benefit$all, type="l", lwd=1,col="#bdb5e1") 
lines(dcr_validation_Sig_p_r_OS$net.benefit$radiopathomics, type="l",lwd=1,col="#7ac7e2")
lines(dcr_validation_Sig_p_OS$net.benefit$pathomics, type="l",lwd=1,col="#f7df87")
lines(dcr_validation_Sig_r_OS$net.benefit$radiomics, type="l",lwd=1,col="#e3716e")
legend("topright", legend = c("None", "All", "nomo", "path_sig", "rad_sig"), 
       col = c("#fccccb", "#bdb5e1", "#7ac7e2", "#f7df87", "#e3716e"), lty = 1, lwd = 1, cex = 0.7)
