# MTMR
Multi-reader, Multi-test Analysis Tool

为便于实际应用, 我给出了mtmr_app.R的分析工具, 支持根据已有数据作假设检验和给出最优样本分配方案两个主要功能。还提供了功效计算和相关系数估计两个附带功能。在完成假设检验后,如果希望得到指定的AUC差值的功效, 可以使用函数mtmr_power()来计算功效; 如果在样本量预估时不确定相关系数, 可以先使用相关领域的已有数据来得到相关系数作为输入, 然后再完成样本量的预估。

对于输入数据, 要按照以下格式。设有m 个恶性病人, n 个良性病人, r 个医生, h 种
医疗影像(一般设为2)。那么需要两个打分矩阵X(𝑚+𝑛)×𝑟ℎ 和Y(𝑚+𝑛)×𝑟ℎ 作为输入。X 对
应恶性病人的得分矩阵, 第i 行即为第i 个病人的所有得分, 其中第(𝑙 − 1)𝑟 + 1 至第𝑙𝑟
列依次为r 名医生使用第l 种方法的打分。Y 对应恶性病人的得分矩阵, 含义同上。
下面对这几个功能对应的函数给出详细说明。

• mtmr_test(): 完成假设检验的函数, 输入参数为
x_rate: 恶性病人的得分矩阵(必须给出)
y_rate: 良性病人的得分矩阵(必须给出)
level: 置信度, 默认值为0.95
times: 如果选用bootstrap, 随机模拟的次数(默认为20000, 建议不要再大, 否则
运行较慢)
h: 医疗影像的种类数(一般设置为2)
fixed_reader: bool 值(默认为TRUE)。TRUE 表示认为医生是固定的, 否则认为医
生是随机选取的。
输出结果为一个list:
第一项为检验结论(那个影像方法更好), 第二项为AUC 差值的区间估计, 第三项
为第一种影像方法辅助下的AUC 的区间估计, 第四项为第二种影像方法辅助下
的AUC 的区间估计。

• mtmr_power(): 完成功效估计的函数, 输入参数为
r: 医生数量(必须给出)
m: 恶性病例数量(必须给出)
n: 良性病例数量(必须给出)
AUC: 不同影像方法辅助下的AUC 均值(必须给出)
delta: 认为两种影像有区别时的AUC 差值的阈值(必须给出)
