# -*- coding: utf-8 -*-
import random
import sys
import copy
import math
import numpy.matlib as np
import poptools
import GPCDE_conf as conf

def get_init_pop(size, genecount, evaluator):
    ### 功能： 得到初始种群 ###
    pop = []
    for i in range(size):
        ind = {}
        ind['extrainfo'] = {}
        ind['vioclassify'] = {}
        ind['genes'] = [random.random() for j in range(genecount)]
        ind['extrainfo']['generation'] = 0
        pop.append(ind)
    #poptools.evaluate_pop_automatic(pop, evaluator, _fill_result)
    poptools.evaluate_pop(pop, evaluator)
    return pop

def get_init_pop_lhs(size, genecount, evaluator, upper, lower):
    ### 功能： 根据拉丁超立方采样得到初始种群 ###
    pop = []
    for i in range(size):
        ind = {}
        ind['extrainfo'] = {}
        ind['vioclassify'] = {}
        ind['genes'] = [random.random() for j in range(genecount)]
        ind['extrainfo']['generation'] = 0
        pop.append(ind)

    select_list = []
    segmentSize = 1.0 / size
    for i in range(genecount):
        tmp = []
        for j in range(size):
            segmentMin = j * segmentSize
            tmp.append(segmentMin + (random.random() * segmentSize))
        select_list.append(tmp)

    for i in range(size):
        for j in range(genecount):
            n = random.randint(0, len(select_list[j])-1)
            pop[i]['genes'][j] = select_list[j][n]
            del select_list[j][n]
    pop = cal_pheno(pop, upper, lower, genecount, len(pop))
    poptools.evaluate_pop_automatic(pop, evaluator, _fill_result)
    return pop

def cal_pheno(pop, upper, lower, n, size):
    ### 功能： 计算pheno值 ###
    for k in xrange(size):
        pop[k]["pheno"] = []
        for i in xrange(n):
            t = pop[k]["genes"][i] * (upper[i] - lower[i]) + lower[i]
            pop[k]["pheno"].append(t)
    return pop

def cal_violation(pop):
    ### 功能： 计算违约值 ###
    m = len(pop[0]["constraint"])
    for i in xrange(len(pop)):
        vio = 0.0
        for j in xrange(m):
            if pop[i]["constraint"][j] > vio:
                vio = pop[i]["constraint"][j]
        pop[i]["violation"] = vio
    return pop

def find_best_obj_vio(_pop):
    ### 功能： 找到当前最好个体的目标值和违约值 ###
    pop = copy.deepcopy(_pop)
    pop.sort(cmp = poptools.compare)
    obj_min = pop[0]["objective"]
    vio_min = pop[0]["violation"]
    return obj_min, vio_min

def init_model(cluster_pop, constraints_num, pin, index):
    ### 功能： 初始化模型（返回的依次为：协相关函数、目标值、协相关矩阵） ###
    #_one = np.ones((_pop_model_size, 1))
    _C = [np.zeros((pin, pin)) for kk in xrange(constraints_num+1)]
    _y = [np.zeros((pin, 1)) for kk in xrange(constraints_num+1)]
    _r = [np.zeros((pin, 1)) for kk in xrange(constraints_num+1)]
    for j in xrange(0, pin):
        _y[0][j, 0] = cluster_pop[index][j]['objective']
        for kk in xrange(constraints_num):
            _y[kk+1][j, 0] = cluster_pop[index][j]["constraint"][kk]
    return _C, _y, _r

def cal_R(i, j, ind, pop):
    ### 功能： 求两个个体之间的相关函数 R(i,j) ###
    ### R(i,j)=correlation function Corr(xi,xj), i,j=1,2,...,n. ###
    genes = ind['genes']
    genecount = len(ind["genes"])/2
    genes_theta = genes[:genecount]
    genes_pi = genes[-genecount:]
    theta = [2 * k for k in genes_theta]
    pi = [k + 1 for k in genes_pi]
    tmp = 0.0
    for m in xrange(0, genecount):
        tmp += theta[m] * math.pow(math.fabs((pop[i]['genes'][m] -pop[j]['genes'][m])), pi[m])
    tmpr = math.exp(-tmp)
    return tmpr

def cal_r(i, ind, best_theta_pi_ind, pop):
    ### 功能： 计算r矩阵 ###
    ### r=[Corr(x*,x1), Corr(x*,x2),..., Corr(x*,xn)] ###
    pheno = ind["genes"]
    genecount = len(ind["genes"])
    b_genes = best_theta_pi_ind['genes']
    b_genes_theta = b_genes[:genecount]
    b_genes_pi = b_genes[-genecount:]
    b_theta = [2 * k for k in b_genes_theta]
    b_pi = [k + 1 for k in b_genes_pi]
    tmp = 0.0
    for m in xrange(0, genecount):
        tmp += b_theta[m] * math.pow(math.fabs((pheno[m] - pop[i]['genes'][m])), b_pi[m])
    r = math.exp(-tmp)
    return r

def cal_mu(C, f, _one):
    ### 功能： 计算mu ###
    tmp1 = np.dot(np.dot(_one.T, C.I), f)
    tmp2 = np.dot(np.dot(_one.T, C.I), _one)
    mu  = tmp1/ tmp2
    return mu

def cal_sigma(f, mu, C, _one, _model_size):
    ### 功能： 计算sigma的平方 ###
    tmp3 = np.dot( np.dot( (f - np.dot(_one, mu)).T, C.I), (f - np.dot(_one, mu)) )
    sigma2 = tmp3[0, 0] / _model_size
    return sigma2

def cal_s2(C, sigma2, _r, _one):
    ### 功能： 计算方差 ###
    tmp4 = np.dot(np.dot(_r.T, C.I), _r)
    tmp5 = np.dot(np.dot(_one.T, C.I), _r)
    tmp8 = np.dot(np.dot(_one.T, C.I), _one)
    tmp6 = (1.0 - tmp5[0,0])**2
    tmp7 = 1.0 - tmp4[0,0] + tmp6 / tmp8[0,0]
    s2 = sigma2 * tmp7
    if s2 <= 0.0:
        s2 = 1.0e-15
    uncertainty = math.sqrt(s2)
    return uncertainty

def cal_predict_f(C, mu, f, _r, _one):
    ### 功能： 计算预测目标均值 ###
    tmp9 = np.dot(np.dot(_r.T, C.I), (f - np.dot(_one, mu)))
    predict_y = mu + tmp9[0, 0]
    return predict_y

def _fill_result(inds, rsts):
    for rst in rsts:
        i = rst['id']
        ind = inds[i]
        if not rst['valid']:
            #ind['violation'] = sys.float_info.max
            ind['objective'] = sys.float_info.max
        else:
            #ind['violation'] = rst['violation']
            ind["constraint"] = rst["constraints"]
            ind['objective'] = rst['objectives'][0]
        if rst.get('extrainfo'):
            ind['extrainfo'].update(rst['extrainfo'])

def evaluate_EI(pop_child, evaluator_EI, obj_min, vio_min, fill_result = _fill_result):
    ### 功能： 根据EI求适应值 ###
    results = []
    for i in range(len(pop_child)):
        pop_child[i]['id'] = i
        results.append(evaluator_EI(pop_child[i], obj_min, vio_min))
    fill_result(pop_child, results)
            
def init_theta_p(size, genecount, evaluator, flag):
    ### 功能：初始化theta和p ###
    global _pop_theta_p
    _pop_theta_p = poptools.get_init_pop(size, genecount * 2)
    for p in _pop_theta_p:
        p['extrainfo']["flag"] = flag
    poptools.evaluate_pop(_pop_theta_p, evaluator, fill_result = _fill_result)
    return _pop_theta_p

F, CR   = conf.F, conf.CR
def DE_cal_theta_p(generation, pop, size, evaluator, genecount):
    ### 功能： 根据差分演化算法求解theta和p #######
    _tmp = copy.deepcopy(pop)
    for g in xrange(1, generation + 1):
        for i in xrange(size):
            n = random.sample(range(size), 3)
            j = random.randint(0, genecount - 1)
            for k in xrange(genecount):
                ra = random.random()
                if ra > CR or k == genecount - 1:
                    x = pop[n[0]]['genes'][j]
                    y = pop[n[1]]['genes'][j]
                    z = pop[n[2]]['genes'][j]
                    r = x + (y - z) * F
                    if r > 1 or r < 0:
                        r = random.random()
                    _tmp[i]['genes'][j] = r
                else:
                    _tmp[i]['genes'][j] = pop[i]['genes'][j]
                j = (j + 1) % genecount
        #_tmp = cal_pheno(_tmp, _upper, _lower, genecount/2, len(_tmp))
        poptools.evaluate_pop(_tmp, evaluator, fill_result = _fill_result)
        for i in range(size):
            if _tmp[i]["objective"] < pop[i]["objective"]:
                pop[i], _tmp[i] = _tmp[i], pop[i]
    return pop

def DE_generate_offspring(generation, pop, size, evaluator, genecount, obj_min, vio_min, _upper, _lower):
    ### 功能： 根据差分演化算法获得子代种群，目标函数为EI（极大化） ####
    _tmp = copy.deepcopy(pop)
    for g in xrange(1, generation + 1):
        for i in xrange(size):
            n = random.sample(range(size), 3)
            j = random.randint(0, genecount - 1)
            for k in xrange(genecount):
                ra = random.random()
                if ra > CR or k == genecount - 1:
                    x = pop[n[0]]['genes'][j]
                    y = pop[n[1]]['genes'][j]
                    z = pop[n[2]]['genes'][j]
                    r = x + (y - z) * F
                    if r > 1 or r < 0:
                        r = random.random()
                    _tmp[i]['genes'][j] = r
                else:
                    _tmp[i]['genes'][j] = pop[i]['genes'][j]
                j = (j + 1) % genecount
            _tmp[i]['extrainfo']['generation'] = g
        _tmp = cal_pheno(_tmp, _upper, _lower, genecount, len(_tmp))
        evaluate_EI(_tmp, evaluator, obj_min, vio_min)
        for i in range(size):
            if _tmp[i]["objective"] > pop[i]["objective"]:
                pop[i], _tmp[i] = _tmp[i], pop[i]
    return pop
