import random
import sys
import os
import shutil
#from deap import dtm

try:
    import global_conf
    import cst#zhoudong
except:
    pass

def get_init_pop(size, genecount):
    pop = []
    for i in range(size):
        ind = {}
        ind['extrainfo']={}
        ind['vioclassify']={}
        ind['genes'] = [random.random() for j in range(genecount)]
        ind['extrainfo']['generation']=0
        pop.append(ind)
    return pop

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
            

def evaluate_pop_automatic(pop, evaluator, fill_result=_fill_result):
    '''try:
        workerId = dtm.getWorkerId()
    except AttributeError:
        workerId = -1
        
    if workerId < 0:
        evaluate_pop(pop, evaluator, fill_result)
    else:
        evaluate_pop_imap_unordered(pop, evaluator, fill_result)
        '''
    evaluate_pop(pop, evaluator, fill_result)

def evaluate_pop(pop, evaluator, fill_result=_fill_result):
    results = []
    for i in range(len(pop)):
        pop[i]['id'] = i
        results.append(evaluator(pop[i]))
    fill_result(pop, results)

def evaluate_pop_imap_unordered(pop, evaluator, fill_result=_fill_result):
    iters = []
    for i in range(len(pop)):
        pop[i]['id'] = i

        # robust
        if global_conf.ROBUST_ENABLE:
            from robust import robust_problem
            iters += robust_problem.robust(pop[i])
        # normal
        else:
            iters.append(pop[i])
            
    imapNotOrderedObj = dtm.imap_unordered(evaluator, iters, len(iters))
    results = [i for i in imapNotOrderedObj]

    # robust result select
    if global_conf.ROBUST_ENABLE:
        from robust import robust_problem
        results = robust_problem.select(results)

    ##zhoudong appended start
    for i in range(len(results)):
        if "err" in results[i].keys():
            results[i]=cst.getresult(results[i])
    
    fill_result(pop, results)

def evaluate_pop_map(pop, evaluator, fill_result=_fill_result):
    for i in range(len(pop)):
        pop[i]['id'] = i
    results = dtm.map(evaluator, pop)
    fill_result(pop, results)

def compare(a, b):
    ret = None
    if abs(a['violation'] - b['violation']) < 1e-6:
        if abs(a['objective'] - b['objective']) < 1e-6:
            ret = 0
        elif a['objective'] < b['objective']:
            ret = -1
        else:
            ret = 1
    elif a['violation'] < b['violation']:
        ret = -1
    else:
        ret = 1
    return ret

def condition(pop):
    pop.sort(cmp=compare)
    if pop[0]['violation'] < 1e-6:
        return True

def copyFiles(sourceDir,  targetDir): 
      for file in os.listdir(sourceDir):
          sourceFile = os.path.join(sourceDir, file)
          targetFile = os.path.join(targetDir, file) 
          if os.path.isfile(sourceFile): 
              if not os.path.exists(targetDir):  
                  os.makedirs(targetDir)  
              if not os.path.exists(targetFile) or(os.path.exists(targetFile) and (os.path.getsize(targetFile) != os.path.getsize(sourceFile))):  
                  open(targetFile, "wb").write(open(sourceFile, "rb").read()) 
          if os.path.isdir(sourceFile): 
              First_Directory = False 
              copyFiles(sourceFile, targetFile)

def getindfoldername(filename):#2_1
    if os.path.basename(os.path.dirname(filename))=='result':
        return os.path.basename(filename)
    else:
        flname=os.path.dirname(filename)
        return getindfoldername(flname)
def getindfolderpath(filename):#result\2_1
    if os.path.basename(os.path.dirname(filename))=='result':
        return filename
    else:
        flname=os.path.dirname(filename)
        return getindfolderpath(flname)

              
if __name__ == '__main__':
    a = {}
    b = {}
    a['violation'] = 5
    b['violation'] = 5
    a['objective'] = 7
    b['objective'] = 5
    print compare(a, b)
