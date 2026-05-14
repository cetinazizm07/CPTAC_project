library(reticulate)

# 1. Yeni oluşturduğumuz sorunsuz Conda ortamını kullanacağımızı kesinleştiriyoruz
use_condaenv("cptac_stable_env", required = TRUE)

py_run_string("
import cptac
import pandas as pd

co = cptac.Coad()

# Proteomik için çalışan 'bcm' kaynağını kullanmaya devam ediyoruz
prot = co.get_proteomics(source='bcm') 

# Klinik veri için jeneratör hatasını bypass eden otomatik kaynak bulucu:
clin = None
olasi_kaynaklar = ['washu', 'bcm', 'umich', 'pdc', 'mssm', 'harmonized', 'broad']

for kaynak in olasi_kaynaklar:
    try:
        clin = co.get_clinical(source=kaynak)
        print(f'Başarılı! Klinik veri şu kaynaktan çekildi: {kaynak}')
        break
    except Exception:
        # Eğer bu kaynakta hata verirse (veya yoksa) sonrakine geç
        continue

# Sütunlardaki MultiIndex yapısını R'ın anlayabileceği düz bir metne çeviriyoruz
if isinstance(prot.columns, pd.MultiIndex):
    prot.columns = ['_'.join(filter(None, map(str, col))).strip() for col in prot.columns.values]
    
if clin is not None and isinstance(clin.columns, pd.MultiIndex):
    clin.columns = ['_'.join(filter(None, map(str, col))).strip() for col in clin.columns.values]
")

# Python'daki (py$) verileri R ortamındaki değişkenlere (data.frame) aktarıyoruz
prot_data <- py$prot
clin_data <- py$clin

cat('Proteomik Verisi Boyutları:', dim(prot_data), '\n')
cat('Klinik Veri Boyutları:', dim(clin_data), '\n')