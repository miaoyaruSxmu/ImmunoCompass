library(caret)        # 用于数据划分
library(glmnet)       # 用于 Elastic Net
library(randomForest) # 用于随机森林
library(tensorflow)
library(xgboost)      # 用于 XGBoost
library(e1071)        # 用于 SVR
library(brms)         # 用于贝叶斯回归 (功能强大，基于Stan)
library(keras)        # 用于神经网络
library(tidyverse)    # 用于数据处理和可视化 (ggplot2)

root <- gsub("\\\\", "/", Sys.getenv("IBI_TOOL_ROOT", unset = "."))
result_dir <- file.path(root, "Result")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(root, "X_clean.rda"))
load(file.path(root, "y.rda"))
X= as.matrix(X_clean[, !apply(X_clean, 2, function(x) any(is.na(x) | is.infinite(x)))])
train_indices <- createDataPartition(as.factor(y), p = 0.7, list = FALSE)

X_train <- X[train_indices, ]
y_train <- y[train_indices]

X_test <- X[-train_indices, ]
y_test <- y[-train_indices]

# 将矩阵转换为数据框，因为一些模型更喜欢数据框
df_train <- as.data.frame(X_train)
df_train$y <- y_train
df_test <- as.data.frame(X_test)
df_test$y <- y_test
# --- 存储所有模型的预测结果 ---
predictions <- data.frame(
  True_Index = y_test
)

# --- 1. Elastic Net (使用交叉验证调优) ---
cat("\nTraining Elastic Net...\n")
# 使用交叉验证找到最佳的lambda
Elastic_model <- caret::train(
  x = X_train,
  y = y_train,
  method = "glmnet",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(
    alpha = 0.5,  # 0表示岭回归 (L2正则化)
    lambda = 10^seq(-3, 3, length = 100)
  )
)

predictions$ElasticNet <- predict(Elastic_model, X_test)

cat("\nTraining ridge...\n")
# 使用交叉验证找到最佳的lambda
ridge_model <- caret::train(
  x = X_train,
  y = y_train,
  method = "glmnet",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(
    alpha = 0,  # 0表示岭回归 (L2正则化)
    lambda = 10^seq(-3, 3, length = 100)
  )
)

predictions$ridge <- predict(ridge_model, X_test)

cat("\nTraining lasso...\n")
# 使用交叉验证找到最佳的lambda
lasso_model <- caret::train(
  x = X_train,
  y = y_train,
  method = "glmnet",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(
    alpha = 1,  # 0表示岭回归 (L2正则化)
    lambda = 10^seq(-3, 3, length = 100)
  )
)

predictions$lasso <- predict(lasso_model, X_test)

# --- 2. 随机森林回归 ---
cat("Training Random Forest...\n")
rf_model <- randomForest(y ~ ., data = df_train, ntree = 500)
predictions$RandomForest <- predict(rf_model, newdata = df_test)


# --- 3. XGBoost 回归 ---
cat("Training XGBoost...\n")
dtrain <- xgb.DMatrix(data = X_train, label = y_train)
dtest <- xgb.DMatrix(data = X_test, label = y_test)
xgb_model <- xgboost(data = dtrain, nrounds = 100, objective = "reg:squarederror", verbose = 0)
predictions$XGBoost <- predict(xgb_model, dtest)
save(xgb_model, file = file.path(root, "xgb_model.rda"))

# --- 4. 支持向量机回归 (SVR) ---
cat("Training SVR...\n")
svr_model <- best_svr#svm(y ~ ., data = df_train, type = "eps-regression", kernel = "radial")
predictions$SVR <- predict(svr_model, newdata = df_test)


# # --- 5. 贝叶斯回归 (使用 brms) ---
# # brms非常强大，会进行MCMC采样，速度较慢
# cat("Training Bayesian Regression (brms)... This may take a moment.\n")
# # formula: y ~ . 表示使用所有其他列作为预测变量
# # family=gaussian() 表示我们假设残差是正态分布的，适合回归
# # silent=2 抑制编译信息
# brm_model <- brm(formula = y ~ .,
#                  data = df_train,
#                  family = gaussian(),
#                  chains = 2, iter = 2000, warmup = 500,
#                  silent = 2, refresh = 0, # 减少输出信息
#                  cores = parallel::detectCores() - 1)
# # 预测会给出完整的后验分布，我们取均值作为点估计
# brm_pred <- predict(brm_model, newdata = df_test)
# predictions$Bayesian <- brm_pred[, "Estimate"]


#predictions=predictions[,-6]

# 1. 计算性能指标
# 将长格式数据转换为适合计算的格式
library(dplyr)
library(tidyverse)
results_long <- predictions %>%
  pivot_longer(cols = -True_Index, names_to = "Model", values_to = "Predicted_Index")

# 计算每个模型的RMSE, MAE和Spearman相关性
performance_summary <- results_long %>%
  group_by(Model) %>%
  summarise(
    RMSE = sqrt(mean((Predicted_Index - True_Index)^2)),
    MAE = mean(abs(Predicted_Index - True_Index)),
    Spearman_Corr = cor(Predicted_Index, True_Index, method = "spearman")
  ) %>%
  arrange(desc(Spearman_Corr)) # 按相关性排序

cat("\n--- Model Performance Summary ---\n")
print(performance_summary)


# 2. 可视化比较
# 箱线图：展示每个模型对不同真实标签的预测值分布
# 一个好的模型应该能清晰地分开-1, 0, 1三组的预测分布
plot <- ggplot(results_long, aes(x = factor(True_Index), y = Predicted_Index, fill = factor(True_Index))) +
  geom_boxplot() +
  facet_wrap(~Model, nrow = 2) +
  labs(
    title = "Comparison of Model Predictions by True Immune Index",
    x = "True Immune Index",
    y = "Predicted Immune Index Score"
  ) +
  scale_fill_discrete(name = "True Index") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)




# 假设您已经运行了之前的模型训练代码，并且工作环境中已经存在`predictions`数据框

# 安装 fmsb 包，用于绘制雷达图
# install.packages("fmsb")
library(fmsb)
library(tidyverse) # 确保 tidyverse 已加载

# 1. 计算所有性能指标
# 将长格式数据转换为适合计算的格式
results_long <- predictions %>%
  # 贝叶斯模型有时会因MCMC问题产生NA，安全起见先过滤掉
  na.omit() %>%
  pivot_longer(cols = -True_Index, names_to = "Model", values_to = "Predicted_Index")

# 计算RMSE, MAE, Spearman Correlation 和新增的 Accuracy
full_performance_summary <- results_long %>%
  group_by(Model) %>%
  summarise(
    RMSE = sqrt(mean((Predicted_Index - True_Index)^2)),
    MAE = mean(abs(Predicted_Index - True_Index)),
    Spearman_Corr = cor(Predicted_Index, True_Index, method = "spearman"),
    # 将预测值四舍五入到最近的整数 (-1, 0, 1) 来计算分类准确率
    Accuracy = mean(round(Predicted_Index) == True_Index)
  )

cat("\n--- Full Model Performance Summary ---\n")
print(full_performance_summary)


# 2. 准备用于雷达图的数据
# 将模型名称作为行名
performance_for_radar <- as.data.frame(full_performance_summary)
rownames(performance_for_radar) <- performance_for_radar$Model
performance_for_radar <- performance_for_radar %>% select(-Model)

# 数据缩放函数
# 对于“越高越好”的指标 (e.g., Accuracy)
scale_higher_is_better <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
# 对于“越低越好”的指标 (e.g., RMSE)
scale_lower_is_better <- function(x) {
  (max(x) - x) / (max(x) - min(x))
}

# 应用缩放
scaled_performance <- performance_for_radar %>%
  mutate(
    RMSE = scale_lower_is_better(RMSE),
    MAE = scale_lower_is_better(MAE),
    Spearman_Corr = scale_higher_is_better(Spearman_Corr),
    Accuracy = scale_higher_is_better(Accuracy)
  )

# fmsb包要求数据框的第1行是最大值，第2行是最小值
# 在我们的缩放后，最大值为1，最小值为0
data_for_radarchart <- rbind(rep(1, ncol(scaled_performance)),
                             rep(0, ncol(scaled_performance)),
                             scaled_performance)

cat("\n--- Scaled Data for Radar Chart ---\n")
print(data_for_radarchart)



# 3. 绘制雷达图
# 定义颜色（可以自定义）
colors_border <- c("#377EB8", "#FF7F00", "#4DAF4A", "#E41A1C", "#984EA3", "#F781BF")
colors_in <- scales::alpha(colors_border, 0.2) # 使用半透明填充色

# 设置图形参数，留出足够空间给图例
par(mar = c(1, 1, 1, 1), xpd = TRUE)

radarchart(
  data_for_radarchart,
  axistype = 1, # 坐标轴样式

  # 自定义多边形
  pcol = colors_border,
  pfcol = colors_in,
  plwd = 2, # 线宽
  plty = 1, # 线型

  # 自定义网格
  cglcol = "grey",
  cglty = 1,
  axislabcol = "grey",
  caxislabels = seq(0, 1, 0.25), # 坐标轴刻度
  cglwd = 0.8,

  # 自定义坐标轴标签
  vlcex = 0.9 # 标签字体大小
)

# 添加图例
legend(
  x = "topright", # 图例位置
  inset = c(-0.1, 0), # 微调图例位置，防止与图形重叠
  legend = rownames(scaled_performance),
  bty = "n",
  pch = 20,
  col = colors_border,
  text.col = "black",
  cex = 1,
  pt.cex = 1.5
)

# 添加标题
title(main = "Model Performance Comparison Radar Chart", cex.main = 1.5, line = 0)






# 假设您已经运行了之前的代码，并生成了 full_performance_summary 数据框
library(tidyverse) # 确保已加载

# 1. 将数据从宽格式转换为长格式
# 这会创建一个包含三列的表: Model, Metric, Value
tidy_performance <- full_performance_summary %>%
  pivot_longer(
    cols = -Model,
    names_to = "Metric",
    values_to = "Value"
  )

# 2. 创建一个标准化的“性能分数”，范围在0-1之间
# 这样，无论指标是“越高越好”还是“越低越好”，我们都可以用统一的颜色方案
# 分数越高（越接近1），代表性能越好
scaled_tidy_performance <- tidy_performance %>%
  group_by(Metric) %>% # 对每个指标独立进行缩放
  mutate(
    # 使用if_else判断指标的方向性
    PerformanceScore = if_else(
      Metric %in% c("RMSE", "MAE"),
      # 对于“越低越好”的指标
      (max(Value, na.rm=T) - Value) / (max(Value, na.rm=T) - min(Value, na.rm=T)),
      # 对于“越高越好”的指标
      (Value - min(Value, na.rm=T)) / (max(Value, na.rm=T) - min(Value, na.rm=T))
    )
  ) %>%
  # 处理当所有值都相同时（分母为0）产生的NaN
  mutate(PerformanceScore = if_else(is.nan(PerformanceScore), 0.5, PerformanceScore)) %>%
  ungroup() # 完成分组操作

# 3. (可选但强烈推荐) 按模型的平均性能对模型进行排序
# 这样表现最好的模型会显示在图的顶部，更易于阅读
model_order <- scaled_tidy_performance %>%
  group_by(Model) %>%
  summarise(AvgScore = mean(PerformanceScore, na.rm=T)) %>%
  arrange(desc(AvgScore)) %>%
  pull(Model)

# 将Model列转换为因子，并按我们计算出的顺序排列
scaled_tidy_performance$Model <- factor(scaled_tidy_performance$Model, levels = rev(model_order))


cat("\n--- Tidy and Scaled Data for Heatmap ---\n")
head(scaled_tidy_performance)



# 4. 使用 ggplot2 绘制热图
heatmap_plot <- ggplot(
    scaled_tidy_performance,
    # 定义x轴、y轴和填充颜色
    aes(x = Metric, y = Model, fill = PerformanceScore)
  ) +
  # 使用 geom_tile 创建热图瓦片
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +

  # 在瓦片上添加原始数值标签，方便精确查看
  geom_text(aes(label = round(Value, 2)), color = "black", size = 4) +

  # 定义一个直观的颜色梯度（从红色到绿色代表从差到好）
  scale_fill_gradient2(
    low = "#F8696B",       # 差的颜色 (红)
    mid = "#FFEB84",       # 中等的颜色 (黄)
    high = "#63BE7B",      # 好的颜色 (绿)
    midpoint = 0.5,        # 颜色梯度的中点
    limits = c(0, 1),      # 颜色范围
    name = "Relative\nPerformance" # 图例标题
  ) +

  # 设置标题和标签
  labs(
    title = "Comprehensive Model Performance Comparison",
    subtitle = "Color indicates relative performance (Green is better), numbers are actual scores.",
    x = "Performance Metric",
    y = "Model"
  ) +

  # 使用一个简洁的主题
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # x轴标签倾斜
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_blank(), # 移除网格线
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# 显示图形
pdf(file.path(result_dir, "model_compare_heatmap.pdf"), 6, 6)
print(heatmap_plot)
dev.off()


# 假设您已经有了 `predictions` 数据框
# 将其转换为长格式以方便使用ggplot的facet_wrap
results_long <- predictions %>%
  na.omit() %>%
  pivot_longer(cols = -True_Index, names_to = "Model", values_to = "Predicted_Index")

# 绘制图形
predicted_vs_actual_plot <- ggplot(results_long, aes(x = True_Index, y = Predicted_Index)) +
  # 添加散点，并用“抖动”防止点完全重叠在-1, 0, 1三条线上
  geom_jitter(aes(color = factor(True_Index)), width = 0.2, height = 0.2, alpha = 0.6) +

  # 添加 y=x 的完美预测线（虚线）
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +

  # （可选）添加一条实际数据的线性回归线（平滑线），看整体趋势
  geom_smooth(method = "lm", aes(color="Trend"), se = FALSE) +

  # 按模型分面，每个模型一张小图
  facet_wrap(~ Model, nrow = 2) +

  # 美化
  labs(
    title = "Predicted vs. Actual Immune Index by Model",
    subtitle = "Points should be close to the dashed y=x line.",
    x = "True Immune Index",
    y = "Predicted Immune Index Score",
    color = "Group"
  ) +
  scale_color_manual(
      name = "Legend",
      values = c("-1" = "#377EB8", "0" = "#4DAF4A", "1" = "#E41A1C", "Trend" = "darkred"),
      labels = c("Immune Deficient", "Healthy", "Autoimmune", "Linear Trend")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")
pdf(file.path(result_dir, "model_compare_dotplot.pdf"), 8, 6)
print(predicted_vs_actual_plot)
dev.off()



# 在长格式数据中计算残差
results_with_residuals <- results_long %>%
  mutate(Residuals = True_Index - Predicted_Index)

# 绘制图形
residuals_plot <- ggplot(results_with_residuals, aes(x = Predicted_Index, y = Residuals)) +
  # 添加水平的 y=0 参考线
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +

  # 添加散点
  geom_point(aes(color = factor(True_Index)), alpha = 0.6) +

  # （可选）添加一条平滑曲线，看残差是否有趋势
  geom_smooth(se = FALSE, color = "darkblue") +

  # 按模型分面
  facet_wrap(~ Model, nrow = 2, scales="free_x") + # free_x允许各图x轴范围不同

  # 美化
  labs(
    title = "Residuals vs. Predicted Values by Model",
    subtitle = "Points should be randomly scattered around the red dashed line.",
    x = "Predicted Immune Index Score",
    y = "Residuals (True - Predicted)",
    color = "True Index Group"
  ) +
  theme_bw(base_size = 12)

print(residuals_plot)



# 我们在之前已经画过箱线图，这里提供一个结合了小提琴图和箱线图的更高级版本
violin_boxplot <- ggplot(results_long, aes(x = factor(True_Index), y = Predicted_Index, fill = factor(True_Index))) +
  # 先画小提琴图来展示分布密度
  geom_violin(trim = FALSE, alpha = 0.5) +

  # 在小提琴图内部叠加一个更简洁的箱线图
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +

  # 按模型分面
  facet_wrap(~ Model, nrow = 2) +

  # 美化
  labs(
    title = "Distribution of Predicted Scores for Each True Group",
    subtitle = "A good model shows clear separation between the three groups.",
    x = "True Immune Index Group",
    y = "Predicted Immune Index Score"
  ) +
  scale_fill_manual(
      name = "True Group",
      values = c("-1" = "#377EB8", "0" = "#4DAF4A", "1" = "#E41A1C"),
      labels = c("Immune Deficient", "Healthy", "Autoimmune")
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") # 通常图例和x轴信息重复，可省略

print(violin_boxplot)
# # --- 6. 神经网络回归 (使用 Keras) ---
# cat("Training Neural Network...\n")
# # 定义模型结构
# nn_model <- keras_model_sequential() %>%
#   layer_dense(units = 32, activation = "relu", input_shape = ncol(X_train)) %>%
#   layer_dropout(rate = 0.2) %>%
#   layer_dense(units = 16, activation = "relu") %>%
#   layer_dense(units = 1, activation = "linear") # 回归问题的输出层

# # 编译模型
# nn_model %>% compile(
#   loss = "mean_squared_error",
#   optimizer = optimizer_adam(),
#   metrics = list("mean_absolute_error")
# )

# # 训练模型
# history <- nn_model %>% fit(
#   X_train, y_train,
#   epochs = 50,
#   batch_size = 32,
#   validation_split = 0.2,
#   verbose = 0
# )

# # 预测
# predictions$NeuralNet <- predict(nn_model, X_test)[,1]


################################################################################
# Part 0: 环境设置与包加载
#
# 首次运行时，请确保安装了所有必要的包。
# install.packages(c("tidyverse", "caret", "glmnet", "randomForest", "xgboost",
#                    "e1071", "MASS", "pROC", "MLmetrics"))
################################################################################

library(tidyverse)    # 用于数据处理和可视化
library(caret)        # 用于数据划分和模型训练
library(glmnet)       # 用于 Elastic Net
library(randomForest) # 用于随机森林
library(xgboost)      # 用于 XGBoost
library(e1071)        # 用于 SVM
library(MASS)         # 用于序数回归 (polr)
library(MLmetrics)    # 用于计算F1分数等

# 清理环境
rm(list = ls())
set.seed(42) # 保证所有随机过程可复现


################################################################################
# Part 1: 数据模拟与准备
#
# 我们将创建一个包含三种有序状态的目标变量 y，并提供其三种格式：
# 1. y_numeric: -1, 0, 1 (用于回归)
# 2. y_factor: "Deficient", "Healthy", "Autoimmune" (用于分类)
# 3. y_ordered: 有序因子 (用于序数回归)
################################################################################

cat("--- Part 1: Simulating Data ---\n")

n_samples <- 400
n_features <- 50
X <- matrix(sample(c(-1, 1), n_samples * n_features, replace = TRUE), nrow = n_samples)
colnames(X) <- paste0("Pair_", 1:n_features)

# 创建一个潜在的“真实”关系
true_coeffs <- c(rep(0.6, 5), rep(-0.4, 3), rep(0, n_features - 8))
linear_combination <- X %*% true_coeffs + rnorm(n_samples, 0, 0.5)

# 根据潜在关系生成三种格式的目标变量
y_numeric <- case_when(
  linear_combination > 0.8  ~ 1,
  linear_combination < -0.8 ~ -1,
  TRUE                      ~ 0
)

y_factor <- factor(y_numeric,
                   levels = c(-1, 0, 1),
                   labels = c("Deficient", "Healthy", "Autoimmune"))

y_ordered <- ordered(y_factor, levels = c("Deficient", "Healthy", "Autoimmune"))

# 数据划分 (70% 训练, 30% 测试)
train_indices <- createDataPartition(y_factor, p = 0.7, list = FALSE)

# 准备所有格式的训练集和测试集
X_train <- X[train_indices, ]
X_test  <- X[-train_indices, ]

# 回归
y_train_numeric <- y_numeric[train_indices]
y_test_numeric  <- y_numeric[-train_indices]

# 分类
y_train_factor <- y_factor[train_indices]
y_test_factor  <- y_factor[-train_indices]

# 序数
y_train_ordered <- y_ordered[train_indices]
y_test_ordered  <- y_ordered[-train_indices]

# 将数据框格式也准备好，方便某些模型使用
df_train <- data.frame(X_train, y = y_train_factor)
df_test  <- data.frame(X_test, y = y_test_factor)


################################################################################
# Part 2: 模型训练、预测与评估
#
# 我们将对三大类方法进行建模，并将所有结果存储在一个列表中。
################################################################################

cat("--- Part 2: Training, Predicting, and Evaluating Models ---\n")

results_list <- list()

# --- 方法一：回归分析 (Regression) ---
cat("Training Regression Models...\n")
# 1.1 Elastic Net
cv_fit <- cv.glmnet(X_train, y_train_numeric, family = "gaussian", alpha = 0.5)
pred_reg_glmnet <- predict(cv_fit, newx = X_test, s = cv_fit$lambda.min)[,1]

# 1.2 随机森林回归
rf_reg_model <- randomForest(x = X_train, y = y_train_numeric, ntree = 500)
pred_reg_rf <- predict(rf_reg_model, newdata = X_test)

# 1.3 XGBoost 回归
dtrain_reg <- xgb.DMatrix(data = X_train, label = y_train_numeric)
dtest_reg  <- xgb.DMatrix(data = X_test, label = y_test_numeric)
xgb_reg_model <- xgboost(data = dtrain_reg, nrounds = 100, objective = "reg:squarederror", verbose = 0)
pred_reg_xgb <- predict(xgb_reg_model, dtest_reg)

# 评估回归模型
regression_predictions <- list(
  Reg_ElasticNet = pred_reg_glmnet,
  Reg_RandomForest = pred_reg_rf,
  Reg_XGBoost = pred_reg_xgb
)

for (model_name in names(regression_predictions)) {
  pred <- regression_predictions[[model_name]]

  # 将回归预测值四舍五入到最近的整数，以计算分类指标
  pred_class <- factor(round(pred), levels = c(-1, 0, 1), labels = c("Deficient", "Healthy", "Autoimmune"))

  results_list[[model_name]] <- data.frame(
    Model = model_name,
    Approach = "Regression",
    RMSE = sqrt(mean((pred - y_test_numeric)^2)),
    Spearman_Corr = cor(pred, y_test_numeric, method = "spearman"),
    Accuracy = Accuracy(pred_class, y_test_factor),
    Macro_F1 = F1_Score(y_true = y_test_factor, y_pred = pred_class, positive = NULL)
  )
}

# --- 方法二：多分类分析 (Multi-Class Classification) ---
cat("Training Classification Models...\n")
# 2.1 随机森林分类
rf_clf_model <- randomForest(x = X_train, y = y_train_factor, ntree = 500)
pred_clf_rf <- predict(rf_clf_model, newdata = X_test)

# 2.2 XGBoost 分类
dtrain_clf <- xgb.DMatrix(data = X_train, label = as.integer(y_train_factor) - 1) # XGBoost需要0,1,2...
dtest_clf  <- xgb.DMatrix(data = X_test, label = as.integer(y_test_factor) - 1)
xgb_clf_model <- xgboost(data = dtrain_clf, nrounds = 100, objective = "multi:softmax", num_class = 3, verbose = 0)
pred_clf_xgb_numeric <- predict(xgb_clf_model, dtest_clf)
pred_clf_xgb <- factor(pred_clf_xgb_numeric, levels = 0:2, labels = c("Deficient", "Healthy", "Autoimmune"))

# 评估分类模型
classification_predictions <- list(
  Clf_RandomForest = pred_clf_rf,
  Clf_XGBoost = pred_clf_xgb
)

for (model_name in names(classification_predictions)) {
  pred <- classification_predictions[[model_name]]

  # 将分类预测值转换为数值，以计算回归指标
  pred_numeric <- as.numeric(as.character(factor(pred, levels = c("Deficient", "Healthy", "Autoimmune"), labels = c(-1, 0, 1))))

  results_list[[model_name]] <- data.frame(
    Model = model_name,
    Approach = "Classification",
    RMSE = sqrt(mean((pred_numeric - y_test_numeric)^2)),
    Spearman_Corr = cor(pred_numeric, y_test_numeric, method = "spearman"),
    Accuracy = Accuracy(pred, y_test_factor),
    Macro_F1 = F1_Score(y_true = y_test_factor, y_pred = pred, positive = NULL)
  )
}

# --- 方法三：序数回归分析 (Ordinal Regression) ---
cat("Training Ordinal Regression Model...\n")
# 3.1 比例优势逻辑回归 (POLR)
# 安装和加载包 (如果还未加载)
# install.packages("ordinalNet")
# 安装和加载包 (如果还未加载)
# install.packages("ordinalNet")
library(randomForest)
rf_for_selection <- randomForest(x = X_train, y = y_train_factor)
importance <- importance(rf_for_selection)
top_features <- names(sort(importance[,1], decreasing = TRUE)[1:20])

# 使用筛选后的特征重新训练polr
df_train_ord_subset <- data.frame(X_train[, top_features], y = y_train_ordered)
polr_model_subset <- polr(y ~ ., data = df_train_ord_subset, Hess = TRUE)
pred_ord_polr <- predict(polr_model_subset, newdata = data.frame(X_test), type = "class")

# 3. 后续评估
# 这里的代码完全不需要改变，只需将您的原始代码中的 'Ord_ElasticNet' 替换掉 'Ord_POLR'
# 为了完整性，我重新列出
pred_numeric_ord <- as.numeric(as.character(factor(pred_ord_polr, levels = c("Deficient", "Healthy", "Autoimmune"), labels = c(-1, 0, 1))))

# 将结果存入 results_list
results_list[["Ord_POLR"]] <- data.frame(
  Model = "Ord_POLR",
  Approach = "Ordinal",
  RMSE = sqrt(mean((pred_numeric_ord - y_test_numeric)^2)),
  Spearman_Corr = cor(pred_numeric_ord, y_test_numeric, method = "spearman"),
  Accuracy = Accuracy(pred_ord_polr, y_test_factor),
  Macro_F1 = F1_Score(y_true = y_test_factor, y_pred = pred_ord_polr, positive = NULL)
)

# 之后再运行 Part 3 的代码，就能将 Ord_ElasticNet 的结果包含在最终的热图中了

################################################################################
# Part 3: 综合性能比较与可视化
#
# 我们将所有结果合并，并用热图进行可视化。
################################################################################

cat("--- Part 3: Consolidating and Visualizing Results ---\n")

# 合并所有结果
final_summary <- bind_rows(results_list)

cat("\n--- FINAL PERFORMANCE SUMMARY ---\n")
print(final_summary %>% arrange(desc(Macro_F1), desc(Spearman_Corr)))


# --- 准备热图数据 ---
# 1. 转换为长格式
tidy_summary <- final_summary %>%
  pivot_longer(
    cols = c(RMSE, Spearman_Corr, Accuracy, Macro_F1),
    names_to = "Metric",
    values_to = "Value"
  )

# 2. 创建标准化的性能分数 (0-1, 1为最佳)
scaled_tidy_summary <- tidy_summary %>%
  group_by(Metric) %>%
  mutate(
    PerformanceScore = if_else(
      Metric == "RMSE",
      (max(Value, na.rm=T) - Value) / (max(Value, na.rm=T) - min(Value, na.rm=T)), # RMSE越低越好
      (Value - min(Value, na.rm=T)) / (max(Value, na.rm=T) - min(Value, na.rm=T))  # 其他越高越好
    )
  ) %>%
  mutate(PerformanceScore = if_else(is.nan(PerformanceScore), 0.5, PerformanceScore)) %>%
  ungroup()

# 3. 按模型的平均性能排序
model_order <- scaled_tidy_summary %>%
  group_by(Model) %>%
  summarise(AvgScore = mean(PerformanceScore, na.rm=T)) %>%
  arrange(desc(AvgScore)) %>%
  pull(Model)

scaled_tidy_summary$Model <- factor(scaled_tidy_summary$Model, levels = rev(model_order))


# --- 绘制最终的热图 ---
final_heatmap <- ggplot(
    scaled_tidy_summary,
    aes(x = Metric, y = Model, fill = PerformanceScore)
  ) +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  geom_text(aes(label = round(Value, 3)), color = "black", size = 3.5) +
  scale_fill_gradient2(
    low = "#F8696B", mid = "#FFEB84", high = "#63BE7B",
    midpoint = 0.5, limits = c(0, 1), name = "Relative\nPerformance"
  ) +
  facet_grid(~ Approach, scales = "free_x", space = "free_x") + # 按建模方法分面
  labs(
    title = "Ultimate Showdown: Comparing Regression, Classification, and Ordinal Models",
    subtitle = "Color indicates relative performance (Green is best). Numbers are precise scores.",
    x = "Performance Metric", y = "Model"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size=12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title.position = "plot",
    plot.title = element_text(face="bold", size=16),
    legend.position = "right"
  )

print(final_heatmap)
