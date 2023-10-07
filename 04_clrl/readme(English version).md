# Common Latent Representation Learning  #
## Data processing  ##
**Step1. Feature extraction**

In language identification, MFCC or SDC is common used.
For MFCC extraction, we provide the matlab-based vision, i.e., extract_mfcc.m
For SDC extraction, the kaldi-based baseline[asv-subtools](https://github.com/Snowdar/asv-subtools "asv-subtools") can be used. The extracted features are saved in “.mfcc” format, which is convenient for matlab. After extracting the SDC of each audio using kaldi, it can be read using the relevant code provided in get_M&Y/from_kaldi.


**Step2. get GMM mean supervector M and label Y**

- ubm：gmm_em.m
- GMM mean supervector M：spkEnroll.m
- label Y：onehot.m

**ps.**
For more details, please refer to get_M&Y/from_matlab/demo_getM&Y.m


## CLRL for language identification ##
please refer to i-vector/clrl_em.m
demo can be found in demo_PLDA.m or demo_CDS.m


**ps.**

In order to improve the comparability of CLRL, we also provide fa_em.m and ppca_em.m codes to facilitate comparison and verification





