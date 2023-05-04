## coding=utf-8
import random
import sys
import math
import numpy.matlib as np
import GPCDE_conf as conf


class Tsort:
    def __init__(self):
        self.x = []
        self.xn = 0
        self.x_ind = []

    def set_value(self, _x, _xn):
        self.x = [0.0] * _xn
        self.x_ind = [0] * _xn
        self.xn = _xn
        for i in xrange(0, _xn):
            self.x[i] = _x[i]
            self.x_ind[i] = i

    def swap(self, l, h):
        t = self.x[l]
        self.x[l] = self.x[h]
        self.x[h] = t
        tt = self.x_ind[l]
        self.x_ind[l] = self.x_ind[h]
        self.x_ind[h] = tt

    def qsort(self, l, h):
        if h < l + 2:
            return
        e = h
        p = l
        while l < h:
            l += 1
            while l < e and self.x[l] <= self.x[p]:
                l += 1
            h -= 1
            while h > p and self.x[h] >= self.x[p]:
                h -= 1
            if l < h:
                self.swap(l, h)
        self.swap(h, p)
        self.qsort(p, h)
        self.qsort(l, e)


# shape[0] = rows行数
# shape[1] = Columns列数
class FCM:
    def __init__(self):
        self.expo = 2
        self.eps = 0.05
        self.maxIter = 10000
        self.alpha = 2

        self.L1 = conf.L1
        self.L2 = conf.L2
        self.stoped = False
        self.vn = 0   # 簇的个数
        self.U = None # 隶属度
        self.V = None # 簇的中心
        self.X = None
        self.P = []
        self.pin = 0
        self.minJ = 1e9

    def setParameter(self, L1, L2):
        self.L1 = L1
        self.L2 = L2

    def setClusterNum(self):
        if self.X.shape[0] < self.L1:
            self.vn = 1
            self.pin = self.X.shape[0]
        elif self.X.shape[0] == self.L1:
            self.vn = (self.X.shape[0] - self.L1) / self.L2 + 1
            self.pin = self.L1
        else:
            self.vn = (self.X.shape[0] - self.L1) / self.L2 + 2
            self.pin = self.L1

    def initFCM(self, _X):
        self.X = np.matrix(_X)
        self.stoped = False
        self.setClusterNum()
        # self.setParameter(80, 20)
        self.U = np.zeros((self.X.shape[0], self.vn))
        self.V = np.zeros((self.vn, _X.shape[1]))

        for i in range(self.U.shape[0]):
            for j in range(self.U.shape[1]):
                t = random.random()
                self.U[i, j] = t

        for i in range(self.U.shape[0]):
            Sum = 0
            for j in range(self.U.shape[1]):
                Sum += self.U[i, j]
            for j in range(self.U.shape[1]):
                if abs(Sum) < 1e-9:
                    self.U[i, j] = 10000
                    continue
                self.U[i, j] = self.U[i, j] / Sum

    def stepFCM(self):
        D = np.zeros((self.X.shape[0], self.V.shape[0]))
        flag = True
        for i in range(self.U.shape[0]):
            for j in range(self.U.shape[1]):
                self.U[i, j] = pow(self.U[i, j], self.alpha)
        self.getV()
        self.matrixDistance(D)
        Sum = 0.0
        for i in range(self.X.shape[0]):
            for j in range(self.V.shape[0]):
                Sum += D[i, j] * D[i, j] * self.U[i, j]
                D[i, j] = pow(D[i, j], -2 / (self.expo - 1))

        for i in range(self.X.shape[0]):
            s = 0.0
            for j in range(self.V.shape[0]):
                s += D[i, j]
            for j in range(self.V.shape[0]):
                temp = self.U[i, j]
                self.U[i, j] = D[i, j] / s
                if abs(temp - self.U[i, j]) > self.eps:
                    flag = False
        self.stoped = flag
        return Sum

    def getV(self):
        Sum = [0.0] * self.U.shape[1]
        for j in range(self.U.shape[1]):
            Sum[j] = 0.0
            for i in range(self.U.shape[0]):
                Sum[j] += self.U[i, j]

        for j in range(self.V.shape[0]):
            for k in range(self.X.shape[1]):
                self.V[j, k] = 0
                for i in range(self.X.shape[0]):
                    self.V[j, k] += self.U[i, j] * self.X[i, k]
                self.V[j, k] /= Sum[j]

    def matrixDistance(self, D):
        if self.V.shape[1] != self.X.shape[1]:
            print("These Matrixs do not have the same columns")
        for i in range(self.X.shape[0]):
            for j in range(self.V.shape[0]):
                Sum = 0.0
                for k in range(self.X.shape[1]):
                    Sum += pow(self.V[j, k] - self.X[i, k], 2)
                    D[i, j] = Sum ** 0.5

    def mFCM(self, _X):
        self.initFCM(_X)
        oldJ = 0.0
        for i in range(self.maxIter):
            J = self.stepFCM()
            if J < self.minJ:
                self.minJ = J
            if self.stoped:
                break
            if abs(oldJ - J) < 1e-6:
                break
            oldJ = J
        self.getClusterSets()

    def getClusterSets(self):
        xNum = self.X.shape[0]
        tmp = [0.0] * xNum
        self.P = []
        for i in range(self.vn):
            t = Tsort()
            for j in range(xNum):
                tmp[j] = self.U[j, i]
            t.set_value(tmp, xNum)
            t.qsort(0, xNum)
            # tp = [0] * self.L1
            if xNum > self.L1:
                tL = self.L1
            else:
                tL = xNum
            tp = [0] * tL
            for j in range(tL):
                idx = xNum - j - 1
                idx = t.x_ind[idx]
                tp[j] = idx
            self.P.append(tp)

    def FindBestSet(self, var):
        Min = 1e9
        mini = 0
        for i in range(self.vn):
            s = 0.0
            for j in range(self.X.shape[1]):
                s += (self.V[i, j] - var[j]) * (self.V[i, j] - var[j])
            if s < Min:
                Min = s
                mini = i
        return mini

    def getWeightSum(self, var, w):
        dist = [0.0] * self.vn
        distsum = 0.0
        for i in range(self.vn):
            dist[i] = 0.0
            for j in range(self.X.shape[1]):
                dist[i] += (self.V[i, j] - var[j]) * (self.V[i, j] - var[j])
            dist[i] **= 0.5
            distsum += 1 / dist[i]

        w = [0.0] * self.vn
        for i in range(self.vn):
            w[i] = (1 / dist[i]) / distsum


if __name__ == '__main__':
    b = np.matrix([[1,10], [2,8], [3,9], [5,5], [10,1], [8,3], [0.5,1]])
    fcm = FCM()
    fcm.mFCM(b)
    #fcm.FindBestSet([5, 0, 1])
    #print fcm.FindBestSet([2, 100])
    #print a

    print "得到的聚类集合：", fcm.P
    #print "个体：", fcm.X
    print "隶属度：", fcm.U
    print "聚类中心：", fcm.V

    print "得到的聚类集合：", fcm.pin
    print len(fcm.P[0])
