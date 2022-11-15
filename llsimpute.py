'''
작성자 : 오세진
이 스크립트는 missing gene imputation을 위한 LLS algorithm을 python으로 제작한 것임.
LLSimpute는 굉장히 간단한 방식으로 상당히 효과적으로 missing gene을 imputation함.
방식은 pearson correlation을 통해 training data에서 가장 유사한 유전자들 k개를 찾은 후
해당 유전자들에서 linear combination값을 통해 해당 유전자의 값을 유추하는 것임.
'''
import pandas as pd
import random
import numpy as np
def _pairwise_correlation(A, B):
    am = A - np.mean(A, axis=0, keepdims=True)
    bm = B - np.mean(B, axis=0, keepdims=True)
    return am.T @ bm /  (np.sqrt(
        np.sum(am**2, axis=0,
               keepdims=True)).T * np.sqrt(
        np.sum(bm**2, axis=0, keepdims=True)))


def _identify_neighbors(ref=None,missing_genes=None,maxK=100):
    if ref is None:
        print(
            '''
            이 기능은 missing_gene들에 대해 가장 유의미하게 선형적으로 연관된
            유전자들들 찾기위한 것임
            parameters
            ref = 레퍼런스로 사용될 유전자 발현값 매트릭스 index는 유전자 column은 샘플들
            missing_genes = 유추하고자 하는 유전자 list ['a','b','c']
            maxK = 각 missing gene별로 최대 몇 개의 correlative 유전자들을 후보군으로 잡을지를 나타냄.

            # return
            DataFrame으로 column에는 MG, probe가 존재함
            MG는 missing gene들을 의미함
            probe는 string으로 maxK까지의 유전자들이 있음.
            유전자의 정렬 순서는 pearson correlation값의 절댓값이 큰 순으로 정렬됨.
            여기서 probe는 missing gene은 제외됨.
            '''
        )
        return
    # missing gene과 probe matrix 분리
    mg_mat=ref.loc[missing_genes,:].T
    probes=list(set(ref.index)-set(missing_genes))
    probe_mat=ref.loc[probes,:].T
    # Correlation coefficient 계산
    print('Pearson correlation coefficient 산출 중')
    cor_coef=_pairwise_correlation(A=mg_mat.values, B=probe_mat.values)
    cor_coef=np.abs(cor_coef)
    # 각 probe 별로 가장 correlative한 probe gene 고르기
    result=pd.DataFrame([])
    result['MG']=missing_genes
    result['Candidate Probes']=''
    print(f'Missing gene별 cancidate probe (k = {maxK}) 정리 중')
    for row in range(result.shape[0]): #row=0
        idx=np.argsort(-cor_coef[row,:])[:maxK]
        result.loc[row,'Candidate Probes']=','.join([probes[i] for i in idx])
    print(f'1차 선별 완료')
    return result
        
# mgcp=identify_neighbors(ref=ref,missing_genes=missing_genes,maxK=100)
from sklearn import linear_model
from kneed import KneeLocator # elbow point를 찾기위함
from tqdm import tqdm
def _select_best_probes(ref=None,mgcp=None,r2_cut=0.5):
    if ref is None:
        print(
            '''
            이 기능은 identify_neighbors()를 통해 찾아진 missing gene에대해 가장 좋은
            예측값을 제공해줄 것으로 추정되는 probe들 중 error (RMSE) 값을 최소로 하는 모델을 찾기위한 것임.
            # parameters
            ref = 레퍼런스로 사용될 유전자 발현값 매트릭스 index는 유전자 column은 샘플들
            mgcp = identify_neighbors() 후 나오는 mg별 candidate probe들 테이블
            r2_cut = Rsquare cutoff로 bestK로 선정된 값의 score가 지정된 R^2 score 미만이면 bestK 대신 R^2 score를 넘는 것을 bestK로 선택함.
            R^2 score를 넘는 것이 없을 경우엔 그냥 bestK로 내보냄.
            # return
            DataFrame으로 column에는 MG, probe가 존재함
            MG는 missing gene들을 의미함
            probe는 선정된 최종 probe들이 들어가게됨. maxK까지 증가시키면서 R^2 score가 유의미한 증가를 멈추는 부분까지 진행함.
            찾는 방식은 KneeLocator로 진행함. 여기서 계산된 elbow-point k값에서 1을 더해서 k+1개의 probe를 선택하여 모델을 만듬.
            ''')
        return
    # 모델 생성
    clf=linear_model.LinearRegression()
    # 계산 시작
    mgcp_fin=mgcp.copy()
    print('최적 K값 찾는 중')
    for idx in tqdm(range(mgcp.shape[0])): #idx=3
        y=ref.loc[mgcp.loc[idx,'MG'],:].values
        probe=mgcp.loc[idx,'Candidate Probes'].split(',')
        scores=[]
        for k in range(1,len(probe)): #k=1
            #print(k)
            x=ref.loc[probe[:k],:].values
            a=clf.fit(x.T,y)
            score=a.score(x.T,y)
            scores.append(score)
        #adj_score=[(1-((1-scores[idx])*(x.shape[1]-1)/(x.shape[1]-2-idx))) for idx in range(len(scores))]
        if scores[0]!=1: # 데이터 프로세싱을 어떻게 하느냐에 따라 유전자 값이 2개간에 동일할 수도 있음. 그런 것을 계산없이 넘기기위한 스크립트
            bestK = KneeLocator(range(1,len(probe)),scores, S=1.0, curve="concave", direction="increasing").knee
            bestK = min([bestK+1,len(scores)]) # bestK+1과 candidate probe 전체 중 무엇을 쓸지 결정.
        else:
            bestK=1
        if scores[bestK-1]<r2_cut and any(np.array(scores)>=r2_cut):
            bestK=np.where(np.array(scores)>=r2_cut)[0][0]+1
        #print(bestK)
        mgcp_fin.loc[idx,'Candidate Probes']=','.join(probe[:bestK])
        mgcp_fin.loc[idx,'R-square']=scores[bestK-1]
    mgcp_fin=mgcp_fin.rename(columns={'Candidate Probes':'Final Probes'})
    return mgcp_fin

#select_best_probes(ref=ref,mgcp=mgcp)
# target=ref.iloc[:18700,:]
import warnings
def impute_missing_gene(ref=None,target=None,maxK=100,r2_cut=0.5):
    if ref is None:
        print(
            '''
            이 기능은 적정한 양의 probe수를 선택한 후에
            model을 만들어서 새로운 dataset에서 imputation을 진행함.
            # parameter
            ref : reference table로 index는 유전자 id
            target : missing gene이 생긴 target table로 역시 index는 유전자 id
            maxK : 얼만큼의 correlative한 probe에 대해 평가를 할지 확인.
            ''')
        return
    # missing gene 선별
    missing_genes=list(set(ref.index)-set(target.index))
    # missing-gene에대한 correlative probe (cp) 산출
    mgcp=_identify_neighbors(ref=ref,missing_genes=missing_genes,maxK=maxK)
    # missing-gene에대한 best cp 산출
    mgcp=_select_best_probes(ref=ref,mgcp=mgcp,r2_cut=0.5)
    # missing-gene을 target table에서 유추하기
    pred=pd.DataFrame([])
    warnings.simplefilter('ignore')
    clf=linear_model.LinearRegression()
    for idx,gene in enumerate(mgcp['MG']): #gene='100652748'
        probes=mgcp.loc[idx,'Final Probes'].split(',')
        train_x=ref.loc[probes,:].values
        train_y=ref.loc[gene,:].values
        a=clf.fit(train_x.T,train_y)
        target_x=target.loc[probes,:].values
        pred_y=clf.predict(target_x.T)
        pred[gene]=pred_y
    pred=pred.applymap(lambda x: x if x>0 else 0).T # 음수로 나온 것들은 0으로 변환
    pred.columns=target.columns
    target_fin=pd.concat([target,pred],axis=0)
    print('예측 완료. 결과물 반환',target_fin.shape)
    warnings.resetwarnings()
    return target_fin

# target=ref.iloc[:18700,:]
#impute_missing_gene(ref=ref,target=target,maxK=100,r2_cut=0.5)