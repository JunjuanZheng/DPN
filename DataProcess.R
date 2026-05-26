# ==========================================
# STEP 1: Data Preparation for SSN (R script)
# ==========================================
library(tidyr)

# =================================================================
# GSE95849 ExpressionSet 数据预处理脚本
# =================================================================
library(Biobase)
library(dplyr)
library(stringr)
library(limma) # 用于探针去重表达量平均


###########################################################################################
# GSE95849 数据集清洗与标准化（完全对照你的 GSE24290 逻辑流水线）
###########################################################################################
library(GEOquery)
library(dplyr)

# 确保工作目录一致
pathTemp = '/mnt2/wanggd_group/zjj/Part/DiabeticNeuralgia/BlukRNA_GSEdata_DPN'
setwd(pathTemp)

dataset = 'GSE95849'
gplFile = "GPL22448.soft.gz"  # GSE95849 对应的 Phalanx 芯片平台注释文件

# -----------------------------------------------------------------
# 1. 载入/下载数据
# -----------------------------------------------------------------
# 方式 A：直接读取你上一轮中已经加载好的本地 'GSE95849_eSet.Rdata'
load(file.path(pathTemp,'Data',dataset,'GSE95849_eSet.Rdata'))   # 执行后工作区会生成你的 `gset` 列表
gds <- gset[[1]]

# 方式 B（备用）：如果你需要重新从 GEO 下载，请取消下两行的注释：
# gds <- getGEO(dataset, destdir = "./Data")
# gds <- gds[[1]]

# 2. 提取临床表型信息并保存为 CSV
phenno = pData(gds)
write.csv(phenno, paste0('./Data/', dataset, '_phenno.csv'))

# 3. 提取芯片探针表达矩阵
exprset <- exprs(gds) %>% as.data.frame()
cat("原始表达矩阵维度：", dim(exprset), "\n")
print(exprset[1:4, 1:4])

# -----------------------------------------------------------------
# 4. 解析 GPL22448 注释文件
# -----------------------------------------------------------------
# 提示：请确保已将 "GPL22448.soft.gz" 下载并存放在你的 ./Data/ 目录下
lines <- readLines(paste0('./Data/',dataset,'/',gplFile))

# 寻找平台数据表开始的标志
data_start <- which(grepl("!platform_table_begin", lines)) + 1
if(length(data_start) == 0) {
  data_start <- which(grepl("^ID\t", lines))
}

# 读取 GPL 注释表
GPL <- read.table(paste0('./Data/',dataset,'/', gplFile), skip = data_start - 1, 
                  sep = '\t', header = TRUE, fill = TRUE, 
                  quote = "", stringsAsFactors = FALSE)

# ?? 避坑检查步：建议在这里执行 View(GPL) 或 colnames(GPL) 观察一下列名
# 不同的 GPL 平台 Gene Symbol 所在的列可能不同。
# 如果你想用固定的列名提取（更安全）：
# gpl <- GPL[, c("ID", "Gene_Symbol")] 
# 如果依然按你之前的代码通过位置提取（默认前两列为探针ID与基因名）：
gpl <- GPL[, c(1, 4)]
colnames(gpl) = c('ID', 'Gene symbol')
gpl = gpl[-1, ]

# -----------------------------------------------------------------
# 5. 合并与探针去重 (Avereps)
# -----------------------------------------------------------------
exp = exprset
exp$ID <- rownames(exp) # 增加新的一列存放探针ID

# 将表达矩阵与注释表依据 ID 进行 merge
exp_symbol <- merge(exp, gpl, by="ID")
exp_symbol <- na.omit(exp_symbol) # 剔除未能匹配到 Symbol 的探针

# 额外防错过滤：过滤掉 Symbol 为空字符串或 NA 的无效行，避免干扰去重
exp_symbol <- exp_symbol[exp_symbol$`Gene symbol` != "" & !is.na(exp_symbol$`Gene symbol`), ]

# 检查有多少重复的 Gene symbol
cat("重复基因统计：\n")
print(table(duplicated(exp_symbol$`Gene symbol`)))

# 对相同的基因名，其表达值取平均值
exp_unique <- avereps(exp_symbol[, -c(1, ncol(exp_symbol))], ID = exp_symbol$`Gene symbol`)

# 验证去重结果（应该全为 FALSE）
print(table(duplicated(rownames(exp_unique))))

# -----------------------------------------------------------------
# 6. 转换并保存
# -----------------------------------------------------------------
GSE95849 = data.frame(exp_unique)
cat("清洗后的基因级表达矩阵维度：", dim(GSE95849), "\n")

save(GSE95849, file = paste0('Data/', dataset,'/', dataset, '_Expression_20260525.Rdata'))
cat("GSE95849 干净的矩阵已成功保存至 Data/GSE95849_Expression.Rdata ！\n")

# =================================================================
# 1. 环境设置与数据加载
# =================================================================
library(limma)
library(dplyr)

setwd(pathTemp)

# 读取你清洗好并保存的基因表达矩阵 (已转换为Gene Symbol且去重)
# 如果你已经有清洗好的 Rdata，直接 load：
load(paste0(pathTemp,'/Data/' ,dataset,'/', 'GSE95849_Expression_20260525.Rdata')) # 载入后变量名为 GSE95849
expr_mat <- as.matrix(GSE95849)
library(limma)
#  进行 Log2 对数转换 (加 1 是为了防止 log2(0) 报错)
# 这步能把 25000 和 3 这种巨大差异，压缩到大概 14 和 2 的正常分布区间
expr_mat_log <- log2(expr_mat + 1)
#  进行分位数标准化
# 这步能强制所有样本(18个列)的表达量中位数和整体分布保持一致，消除系统误差
expr_mat_norm <- normalizeBetweenArrays(expr_mat_log)
# 检查标准化后的数据
head(expr_mat_norm)
# 接下来，将标准化后的矩阵赋给你要往下走的变量
expr_mat <- expr_mat_norm

# 读取临床表型数据
meta_data <- read.csv(paste0(pathTemp,'/Data/' ,dataset,'/','GSE95849_phenno.csv'), header=TRUE, row.names=1, stringsAsFactors=FALSE)

# =================================================================
# 2. 导入 STRING 背景网络并过滤表达矩阵
# =================================================================
bg <- read.table(paste0(pathTemp,'/Data/StringDataBase/', "9606.protein.links.v12.0.symbols.txt"), header=F, stringsAsFactors=FALSE)
valid_genes <- intersect(row.names(expr_mat), unique(c(bg$V1, bg$V2)))
expr_mat_filtered <- expr_mat[valid_genes, ]

# =================================================================
# 3. 严格定义分组与进展严重程度评分 (Severity)
# =================================================================
# 严格匹配你的 table() 输出结果
control_string <- "healthy"
dm_string      <- "diabetes mellitus"
dpn_string     <- "Diabetic peripheral neuropathy"

# 提取样本 ID
ref_samples    <- rownames(meta_data)[meta_data$disease.state.ch1 == control_string]
dm_samples     <- rownames(meta_data)[meta_data$disease.state.ch1 == dm_string]
dpn_samples    <- rownames(meta_data)[meta_data$disease.state.ch1 == dpn_string]
sample_samples <- c(dm_samples, dpn_samples) # 12个待测样本

# 构建下游分析必需的临床表型进展矩阵
sub_meta <- meta_data[sample_samples, "disease.state.ch1", drop=FALSE]
colnames(sub_meta) <- "Group"
sub_meta$Severity <- ifelse(sub_meta$Group == dm_string, 1, 2) # DM=1, DPN=2
write.table(sub_meta, file=paste0(pathTemp,'/Data/' ,dataset,'/',"SSN_Input/GSE95849_disease_severity.txt"), quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)

# =================================================================
# 4. 矩阵分离、asinh平滑变换与文件导出
# =================================================================
#dir.create("SSN_Input", showWarnings = FALSE)

# 4.1 导出参考组 (6个健康样本)
refs <- expr_mat_filtered[, ref_samples]
refs_asinh <- asinh(refs)
refs_df <- cbind(data.frame(ID=rownames(refs_asinh)), refs_asinh)
write.table(refs_df, file=paste0(pathTemp,'/Data/' ,dataset,'/',"SSN_Input/reference.txt"), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

# 4.2 导出测试组 (6个DM + 6个DPN = 12个样本)
samples <- expr_mat_filtered[, sample_samples]
samples_asinh <- asinh(samples)
samples_df <- cbind(data.frame(ID=rownames(samples_asinh)), samples_asinh)
write.table(samples_df, file=paste0(pathTemp,'/Data/' ,dataset,'/',"SSN_Input/samples.txt"), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

# 4.3 导出过滤后的背景网络
bg_filtered <- bg[(bg$V1 %in% valid_genes) & (bg$V2 %in% valid_genes), ]
write.table(bg_filtered[, c("V1", "V2")], 
            file=paste0(pathTemp,'/Data/' ,dataset,'/',"SSN_Input/STRING_background.txt"), 
            row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

cat("【步骤一完成】Python 所需输入文件已成功保存在 SSN_Input/ 目录下！\n")



# ==========================================
# STEP 2: 运行 Python 脚本构建单样本网络
# ==========================================
# https://github.com/xp-liu/SSN 
# 1. 切换到工作目录
cd /mnt2/wanggd_group/zjj/Part/DiabeticNeuralgia
# 2. 创建输出目录
mkdir -p SSN_Output
# 3. 运行 Python 脚本 (请将 /path/to/ 替换为你本地实际存放 SSN 仓库的绝对路径)
# 3. 运行 Python 脚本 (严格使用等号连接，且加上 -pvalue 参数)
python /mnt2/wanggd_group/zjj/Part/DiabeticNeuralgia/BlukRNA_GSEdata_DPN/Scripts/SSN_AnalysisCode/SSN-master/construct_single_network.py \
  -ref=./Data/GSE95849/SSN_Input/reference.txt \
  -sample=./Data/GSE95849/SSN_Input/samples.txt \
  -background=./Data/GSE95849/SSN_Input/STRING_background.txt \
  -out=./Data/GSE95849/SSN_Output \
  -pvalue=0.05


