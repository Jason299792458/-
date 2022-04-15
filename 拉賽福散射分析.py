#%% 必要模組
import numpy as np
import matplotlib.pyplot as plt
#%% 實驗原始數據
N = [515, 163, 40, 20, 20, 20, 20, 15]
time = [60, 60, 60, 150, 443, 718, 1472, 2735]
theta = [2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20]
#%% 分析所需參數設定
V = []
X = []

for i in range(0, len(N)): #計算N/t和csc^4
    x = 1/np.sin(theta[i]*(np.pi/180))**4
    X.append(x)
    v = N[i]/time[i]
    V.append(v)
    
def LSM(ListX, ListY, n, c): 
#最小平方法
#其中(n為計算資料個數，c為從list後面減少計算資料個數)
    meanX = 0 #x的平均值
    meanY = 0 #y的平均值
    b1 = 0 #相關係數分子
    b2 = 0 #相關係數分母
    res = 0 #殘差平方和
    tot = 0 #總平方和
    for i in range(c, n): #平均值計算
        meanX += (1/(n - c))*ListX[i]
        meanY += (1/(n - c))*ListY[i]
    for i in range(c, n): #相關係數計算
        b1 += (ListX[i] - meanX)*(ListY[i] - meanY)
        b2 += (ListX[i] - meanX)**2
    b = b1/b2
    a = meanY - b*meanX
    for i in range(c, n):#R^2計算
        res += (ListY[i] - (ListX[i]*b + a))**2
        tot += (ListY[i] - meanY)**2
    R = 1 - res/tot
    return [a, b, -a/b, R, meanX, meanY]

o1 = len(N) #計算資料個數
o2 = 1 #list後面減少計算資料個數(通常用於刪除誤差項)
#%% 資料繪圖
print('R^2 =', round(LSM(X, V, o1, o2)[3], 3))

Xlsm = [-1000, 100000] #設定回歸線x的頭尾
Ylsm = [LSM(X, V, o1, o2)[1]*Xlsm[0] + LSM(X, V, o1, o2)[0], #設定回歸線y的頭尾
        LSM(X, V, o1, o2)[1]*Xlsm[1] + LSM(X, V, o1, o2)[0]]

plt.plot(Xlsm, Ylsm, '--', c = 'r')
plt.scatter(X, V, s = 25, c = 'black')
plt.xlabel("1/sin^4",)
plt.ylabel("N/t")
plt.xlim(-1000, 19000)
plt.ylim(-0.1, 3.5)
plt.show()

plt.scatter(theta, V)
plt.xlabel("angle(deg)",)
plt.ylabel("N/t")
plt.show()
    

