# PyASR

**Python을 이용한 조상 서열 복원**

PyASR은 [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) ("Phylogenetic Analysis by Maximum Likelihood")에 대한 현대적인 Python 인터페이스를 제공하며, 특히 조상 단백질/DNA 서열 복원에 맞춰져 있습니다.

**참고:** PyASR은 현재 단백질 복원만 지원합니다. 이 기능은 개발 중입니다.

## 기본 사용법

```python
import phylopandas as pd
import dendropy as d
import pyasr

# phylopandas를 사용하여 조상 세트를 읽습니다.
df_seqs = pd.read_fasta('test.fasta')

# dendropy를 사용하여 트리를 읽습니다.
tree = d.Tree.get(path='tree.newick', schema='newick')

# 트리에서 노드를 복원합니다.
tree, df_seqs, df_anc = pyasr.reconstruct(df_seqs, tree, working_dir='test', alpha=1.235)

# 조상 데이터프레임을 CSV 파일로 저장합니다.
df_anc.to_csv('ancestors.csv')
```

ToyTree 라이브러리 덕분에 JupyterLab 내부에서 트리와 함께 조상을 나란히 시각화할 수 있습니다.

<img src="docs/jlab-example.png" align="middle">

## 설치

이 패키지는 PyPi에 배포되어 있습니다. pip를 사용하여 설치할 수 있습니다:
```
pip install pyasr
```

개발 버전을 얻으려면:
```
git clone
cd
pip install -e .
```

## 의존성

실제 복원 계산은 [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html)을 사용하여 수행됩니다. 이를 위해서는 PAML이 설치되어 있고 `codeml`/`baseml` 실행 파일이 `$PATH` 환경 변수로 내보내져야 합니다. PAML 설치 지침은 PAML 웹사이트에서 찾을 수 있습니다.

PyASR이 작동하려면 다음 Python 의존성이 필요합니다.

- Pandas
- Biopython
- PhyloPandas
- DendroPy
