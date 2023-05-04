####################################################
# Note: This is Expected Improvement of Constraint Violation
# Last modified: 2018-02-05
# Developers: Jiao Ruwang    ruwangjiao@gmail.com
####################################################

# -*- coding: utf-8 -*-
import os
import sys
import random
import copy
import math
import numpy.matlib as np
import GPCDE_tools
import GPCDE_conf as conf
import poptools
from scipy import stats, integrate
import FCM
#import output

WORKING_DIR  = os.getcwd()
PROBLEM_DIR  = WORKING_DIR + r"/PROBLEM"
RESULT_DIR   = WORKING_DIR + r"/RESULT"
LOCAL_PATH   = [WORKING_DIR, PROBLEM_DIR, RESULT_DIR]
sys.path.extend(LOCAL_PATH)


def init(_genecount, _evaluator, _upper, _lower):
    ### 功能： 用拉丁超立方采样初始化种群 ###
    #output._init_output()
    global _pop, pop_size
    pop_size = 11*_genecount - 1
    _pop = GPCDE_tools.get_init_pop_lhs(pop_size, _genecount, _evaluator, _upper, _lower)
    _pop = GPCDE_tools.cal_violation(_pop)
    return _pop

fcm = FCM.FCM()
def fcm_pop(genecount):
    ### 功能： 将种群_pop进行模糊聚类,聚类后的种群为cluster_pop ###
    global cluster_pop
    X = np.zeros((pop_size, genecount))    # X为模糊聚类所需参数矩阵
    for i in xrange(0, pop_size):
        for j in xrange(0, genecount):
            X[i, j] = _pop[i]['genes'][j]
    fcm.mFCM(X)  # 模糊聚类函数接口
    clusterNum = fcm.vn  # 聚类个数
    cluster_pop = []
    print "clusterNum:", clusterNum
    if clusterNum == 1:
        cluster_pop.append(_pop)
    else:
        for i in xrange(0, clusterNum):
            tmp_pop = [{} for ii in xrange(fcm.L1)]
            for j in xrange(0, fcm.pin):
                tmp_pop[j] = _pop[fcm.P[i][j]]
            cluster_pop.append(tmp_pop)

def evaluator_likelihood(ind):
    ### 功能： 求最大似然函数 ###
    index = ind['extrainfo']["flag"]   # kk=0表示目标，否则表示第kk个约束
    for i in xrange(0, fcm.pin):
        for j in xrange(0, i + 1):
            total_R[_numCluster][index][i, j] = GPCDE_tools.cal_R(i, j, ind, cluster_pop[_numCluster])
    for i in xrange(0, fcm.pin):
        for j in xrange(i, fcm.pin):
            total_R[_numCluster][index][i, j] = total_R[_numCluster][index][j, i]
    mu = GPCDE_tools.cal_mu(total_R[_numCluster][index], total_y[_numCluster][index], _one)
    sigma2 = GPCDE_tools.cal_sigma(total_y[_numCluster][index], mu, total_R[_numCluster][index], _one, fcm.pin)

    #### maximizing the likelihood function ###
    val = sigma2 ** (fcm.pin / 2.0) * math.sqrt( math.fabs( np.linalg.det(total_R[_numCluster][index]) ) )
    # val = (fcm.pin / 2.0) * math.log10(sigma2) + math.log10(math.sqrt(math.fabs(np.linalg.det(total_R[_numCluster][index]))))

    r = {'objectives': [val], "constraints": []}
    if ind.get("id") is not None:
        r['id'] = ind['id']
    r['valid'] = True
    return r

def evaluator_EI(ind, obj_best, vio_best):
    ### 功能：用约束EI进行评估， 分为有可行解和无可行解两种情况 ###

    bestSet = fcm.FindBestSet(ind['genes']) # 找到最好的聚类集合
    
    mu, sigma2, uncertainty, predict_y = [], [], [], []

    ### 计算mu ###
    for i in xrange(constraints_num+1):
        temp1 = GPCDE_tools.cal_mu(total_R[bestSet][i], total_y[bestSet][i], _one)
        mu.append(temp1)

    ### 计算sigma的平方 ###
    for i in xrange(constraints_num+1):
        temp2 = GPCDE_tools.cal_sigma(total_y[bestSet][i], mu[i], total_R[bestSet][i], _one, fcm.pin)
        sigma2.append(temp2)
    
    ### 计算协相关矩阵 r=[Corr(x*,x1), Corr(x*,x2),..., Corr(x*,xn)] ###
    for j in xrange(constraints_num + 1):
        for i in xrange(0, fcm.pin):
            total_r[bestSet][j][i, 0] = GPCDE_tools.cal_r(i, ind, total_theta_p[bestSet][j], cluster_pop[bestSet])

    ### 计算方差 ###
    for i in xrange(constraints_num+1):
        temp3  = GPCDE_tools.cal_s2(total_R[bestSet][i], sigma2[i], total_r[bestSet][i], _one)
        uncertainty.append(temp3)
        
    ### 计算预测目标均值 ###
    for i in xrange(constraints_num+1):
        temp4 = GPCDE_tools.cal_predict_f(total_R[bestSet][i], mu[i], total_y[bestSet][i], total_r[bestSet][i], _one)
        predict_y.append(temp4)
        
    ### 无可行解情况下EI的被积函数 ###
    def fun_integral(tempx):
        temp = 1.0
        for i in xrange(1, constraints_num+1):
            temp *= stats.norm.cdf( (tempx-predict_y[i])/uncertainty[i] )
        return temp
    
    ### 存在可行解的情况 ###
    if vio_best == 0.0:
        temp = 1.0
        tempnormal = (obj_best - predict_y[0])/uncertainty[0]
        EIvalue = (obj_best - predict_y[0]) * stats.norm.cdf(tempnormal) + uncertainty[0] * stats.norm.pdf(tempnormal)
        for i in xrange(1, constraints_num+1):
            temp *= stats.norm.cdf(-predict_y[i]/uncertainty[i])
        EIvalue *= temp
    ### 不存在可行解的情况 ###
    else:
        temp = vio_best
        for i in xrange(1, constraints_num+1):
            temp *= stats.norm.cdf(-predict_y[i]/uncertainty[i])
        tmpvalue, err = integrate.quad( fun_integral, 0, vio_best)
        EIvalue = tmpvalue - temp

    r = {'objectives': [EIvalue], 'violation': 0, "constraints":[]}
    if ind.get("id") is not None:
        r['id'] = ind['id']
    r['valid'] = True
    return r

def model_building(evaluator_likelihood):
    ### 功能： 对目标和每个约束分别建模 ###
    best_theta_p = []
    # flag = 0, 代表目标； 否则flag = j, 代表第j个约束
    for mm in xrange(0, constraints_num + 1):
        _pop_theta_p = GPCDE_tools.init_theta_p(conf.popsize, _genecount, evaluator_likelihood, flag = mm)
        _pop_theta_p = GPCDE_tools.DE_cal_theta_p(conf.generation, _pop_theta_p, conf.popsize, evaluator_likelihood, _genecount*2)
        pop_temp = sorted(_pop_theta_p, key = lambda s: s['objective'])
        best_theta_p.append(pop_temp[0])
        for j in xrange(0, fcm.pin):
            for k in xrange(0, fcm.pin):
                total_R[_numCluster][mm][j, k] = GPCDE_tools.cal_R(j, k, best_theta_p[mm], cluster_pop[_numCluster])
        a = 1e-6 + total_R[_numCluster][mm]
        total_R[_numCluster][mm] = np.dot(a, 1.0/(1e-6+1.0))  # 正则化，防止出现奇异矩阵
    total_theta_p.append(best_theta_p)
  

    
def run(problem_initialize, evaluator, outputfreq = 1, current_g=0, MaxK = 1, K = 1, condition = lambda x : False):
    global total_R, _one, total_y, total_r, total_theta_p, _genecount, _evaluator, _upper, _lower, constraints_num
    global lObj, lVio, pop_size, _numCluster
    _evaluator, _genecount, _upper, _lower, constraints_num = evaluator, problem_initialize[0], problem_initialize[1], problem_initialize[2], problem_initialize[4]
    _pop = init(_genecount, _evaluator, _upper, _lower)
    lObj, lVio = [], []
    initModelSize = len(_pop)
    for g in xrange(conf.FEs):
        print "Current generation is:", g
        fcm_pop(_genecount)
        clusterNum = fcm.vn   # 聚类的个数
        obj_best, vio_best = GPCDE_tools.find_best_obj_vio(_pop) # 找出已评估个体中最好的目标值和违约值
        lObj.append('%.8f' % obj_best)
        lVio.append('%.8f' % vio_best)
        print "Current best solution (obj & vio):",  obj_best, vio_best
        _one = np.ones((fcm.pin, 1)) # fcm.pin为每个簇中所包含的个体数
        total_R, total_y, total_r, total_theta_p = [], [], [], []
        for ij in xrange(clusterNum):
            _numCluster = ij # 表示第几个簇
            R, y, r = GPCDE_tools.init_model(cluster_pop, constraints_num, fcm.pin, ij)
            total_R.append(R)
            total_y.append(y)
            total_r.append(r)
            
            ### 对目标和每个约束条件分别建模 ###
            model_building(evaluator_likelihood)

        ### 随机初始化父代，用DE产生子代，评估函数为EI ###
        pop_temp = poptools.get_init_pop(conf.popsize, _genecount)
        pop_temp = GPCDE_tools.cal_pheno(pop_temp, _upper, _lower, _genecount, len(pop_temp))
        GPCDE_tools.evaluate_EI(pop_temp, evaluator_EI, obj_best, vio_best)
        pop_child = GPCDE_tools.DE_generate_offspring(conf.generation, pop_temp, conf.popsize, evaluator_EI, _genecount, obj_best, vio_best, _upper, _lower)

        pop_child_tmp = sorted(pop_child, key = lambda s: s['objective'])
        best_ind = pop_child_tmp[-1] # 选取EI值最大的一个个体
        best_ind['id'] = pop_size
        rst = evaluator(best_ind) # 用真实评估函数评估所选出来的EI值最大的一个个体
        if not rst['valid']:
            pass
        else:
            if rst.get('extrainfo'):
                best_ind['extrainfo'].update(rst['extrainfo'])
            best_ind['objective']  = rst['objectives'][0]
            best_ind['constraint'] = rst['constraints']
            best_ind['valid'] = rst['valid']
            best_ind["pheno"] = []
            for iii in xrange(_genecount):
                t = best_ind["genes"][iii] * (_upper[iii] - _lower[iii]) + _lower[iii]
                best_ind["pheno"].append(t)
            print "best ind pheno:", best_ind["pheno"]
            [best_ind] = GPCDE_tools.cal_violation([best_ind])
        _pop.append(best_ind) # 将真实评估过的个体加入到_pop中
        pop_size += 1
        print 'Predict best solution (obj & vio):', best_ind["objective"], best_ind["violation"]
        print ""
        #output.output(_pop, i, generation_predict, outputfreq)
    return leave()

def leave():
    _pop.sort(cmp = poptools.compare)
    return _pop[0]['objective'], lObj, lVio

def get_average(res):
    ### 功能： 对结果求平均值 ###
    c = sum(res)
    ave = float(c)/len(res)
    return ave

def get_variance(res, ave):
    ### 功能： 对结果求标准差 ###
    sumvar = 0.0
    for i in range(len(res)):
        sumvar = sumvar + pow(float(res[i])-ave,2)
    var = pow(sumvar/len(res), 0.5)
    return var


if __name__ == '__main__':
    import g01, g02, g04, g06, g07, g08, g09, g10, g12, g14, g16, g18, g19, g24
    module =[g06]
    for m in module:
        print "this is",m.__name__,"problem"
        problem_initialize = m.problem_initialize()
        print "D is ",problem_initialize[0]
        t = 25
        res = []
        initFile = open(RESULT_DIR+"/"+str(m.__name__)+".txt",'w')
        print 70*"="
        initFile.write("this is GPCDE:")
        initFile.write('\n')
        initFile.close()
        while t > 0:
            avr=(run(m.problem_initialize(), m.evaluate))
            res.append(avr[0])
            res1 = avr[1]
            res2 = avr[2]
            initFile = open(RESULT_DIR+"/"+str(m.__name__)+".txt",'a')
            initFile.write('run is: '+str(t))
            initFile.write('\n')
            initFile.write(str(avr[0]))
            initFile.write('\n')
            initFile.write("All objective: "+str(res1))
            initFile.write('\n')
            initFile.write("All violation: "+str(res2))
            initFile.write('\n')
            t-=1
            initFile.close()
        initFile = open(RESULT_DIR+"/"+str(m.__name__)+".txt",'a')
        initFile.write('\n')
        average = get_average(res)
        initFile.write('average is: '+ str(average))
        initFile.write('\n')
        vari = get_variance(res, average)
        initFile.write("varience is: " +str(vari))
        print 70*"="
        initFile.close()
