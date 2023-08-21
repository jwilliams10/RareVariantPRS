rm(list = ls())

sampleids_all <- read.table("/data/williamsjacr/UKB_WES_lipids/Data/sampleids.txt", quote="\"", comment.char="")

set.seed(1330)

i <- 1:nrow(sampleids_all)

train_number <- round(max(i)*0.7) + 1
tune_val_number <- round(max(i)*0.15)

train <- sample(i, train_number)

i <- i[!(i %in% train)]

tune <- sample(i, tune_val_number)
validation <- i[!(i %in% tune)]


train <- sampleids_all[train,]
tune <- sampleids_all[tune,]
validation <- sampleids_all[validation,]

write.table(train,"/data/williamsjacr/UKB_WES_lipids/Data/train.txt",row.names = FALSE)
write.table(tune,"/data/williamsjacr/UKB_WES_lipids/Data/tune.txt",row.names = FALSE)
write.table(validation,"/data/williamsjacr/UKB_WES_lipids/Data/validation.txt",row.names = FALSE)
