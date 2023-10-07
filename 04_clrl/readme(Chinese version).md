# Common Latent Representation Learning  #
## 数据处理  ##
**Step1. 特征提取**

在语种识别中，我们经常选择提取音频的MFCC特征/SDC特征。在本代码中，提供了基于matlab的MFCC特征提取的代码*extract_mfcc.m*

对于SDC特征的提取，在此建议使用基于kaldi的基线系统[asv-subtools](https://github.com/Snowdar/asv-subtools "asv-subtools")，并将提取得到的特征保存为.mfcc格式，便于matlab使用。在使用kaldi提取得到每个音频的sdc特征后，可以使用get_M&Y/from_kaldi中提供的相关代码进行读取。

**Step2. 得到GMM均值超矢量M与类别标签Y**

- ubm：gmm_em.m
- GMM均值超矢量M：spkEnroll.m
- 类别标签Y：onehot.m

**ps.**
对于数据处理的具体执行流程，可以参照*get_M&Y/from_matlab*中所提供的相关示例代码：*demo_getM&Y.m*

## CLRL语种识别 ##
我们在ivector文件夹下提供了CLRL的具体实现代码*clrl_em.m*

其具体的实现流程可以参照同文件夹目录下的*demo_PLDA.m* 和 *demo_CDS.m*代码

**ps.**

为了提高CLRL的可对比性，在此我们还提供了fa_em.m 与 ppca_em.m代码，便于进行对比实现验证





