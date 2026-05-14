library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
# 1. Klinik verideki mükerrer sütun isimlerini otomatik olarak benzersiz yapıyoruz
# (Örn: İkinci kez geçen 'history_source' ismi 'history_source...74' gibi benzersiz olur)
colnames(clin_data) <- make.unique(colnames(clin_data))
colnames(prot_data) <- make.unique(colnames(prot_data))

# 2. Artık sütun isimleri tamamen benzersiz olduğu için satır isimlerini güvenle sütuna çevirebiliriz
prot_df <- prot_data %>% rownames_to_column(var = "CPTAC_Barcode")
clin_df <- clin_data %>% rownames_to_column(var = "CPTAC_Barcode")

# 3. Tabloları ortak barkod üzerinden birleştiriyoruz
merged_data <- inner_join(clin_df, prot_df, by = "CPTAC_Barcode")

# 4. İlgilendiğimiz klinik özellikleri ve tüm protein sütunlarını filtreleyip çekiyoruz
final_metadata_protein <- merged_data %>%
  select(
    CPTAC_Barcode,
    age, sex, bmi,
    histologic_type, histologic_grade,
    tumor_stage_pathological,               # Patolojik Evre
    `Survival status (1, dead; 0, alive)`,   # Hayatta kalma durumu
    `Overall survival, days`,               # Hayatta kalma süresi (gün)
    everything()                            # Geri kalan tüm protein sütunları
  )

# Sonucu doğrulayalım
cat("Birleştirilmiş ve Temizlenmiş Veri Seti Boyutu:", dim(final_metadata_protein), "\n")

#EDA
# 1. Sayısal değişkenler (Yaş ve BMI) için genel özet
cat("--- Sayısal Değişkenlerin Özeti ---\n")
final_metadata_protein %>% 
  select(age, bmi) %>% 
  summary()

# 2. Kategorik değişkenler (Cinsiyet ve Tümör Evresi) için frekans tablosu
cat("\n--- Cinsiyet Dağılımı ---\n")
table(final_metadata_protein$sex, useNA = "ifany")

cat("\n--- Tümör Evresi (Stage) Dağılımı ---\n")
table(final_metadata_protein$tumor_stage_pathological, useNA = "ifany")

cat("\n--- Sağkalım Durumu Dağılımı ---\n")
table(final_metadata_protein$`Survival status (1, dead; 0, alive)`, useNA = "ifany")

#Calulating NA ratio of proteins at samples
na_percentages <- colMeans(is.na(final_metadata_protein[, 10:ncol(final_metadata_protein)])) * 100

# NA yüzdesi dağılımını bir histogram ile görelim
ggplot(data.frame(NA_Rate = na_percentages), aes(x = NA_Rate)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Proteinlerdeki Eksik Veri (NA) Oranı Dağılımı",
       x = "Eksik Veri Oranı (%)",
       y = "Protein Sayısı")

#Principal Component Analysis
# 1. Sadece sayısal proteinleri filtrelediğimizden emin oluyoruz
protein_only_numeric <- final_metadata_protein %>% 
  select(where(is.numeric)) %>% 
  select(-any_of(c("age", "bmi", "Overall survival, days", "Survival status (1, dead; 0, alive)")))

# 2. HİÇ eksik verisi olmayan proteinleri seçiyoruz
clean_proteins <- protein_only_numeric %>% 
  select(where(~all(!is.na(.))))

# --- KONTROL ADIMI ---
# Eğer temiz protein tablosunun sütun sayısı 0 ise, hiçbir protein tamamen eksiksiz değildir.
# Bu durumda PCA'in çalışabilmesi için basit bir ortalama ataması (Mean Imputation) yapmalıyız.
if (ncol(clean_proteins) == 0) {
  cat("Uyarı: Hiç eksik verisi olmayan protein bulunamadı. Eksik değerler ortalama ile dolduruluyor...\n")
  
  # Her sütundaki NA değerlerini o sütunun ortalamasıyla dolduruyoruz
  imputed_proteins <- protein_only_numeric
  for(i in 1:ncol(imputed_proteins)) {
    imputed_proteins[is.na(imputed_proteins[,i]), i] <- mean(imputed_proteins[,i], na.rm = TRUE)
  }
  
  # Varyansı sıfır olan (hiç değişmeyen) sabit sütunlar varsa PCA'i bozmaması için eliyoruz
  clean_proteins <- imputed_proteins %>% select(where(~var(.) > 0))
}

# 3. PCA Hesaplama
pca_result <- prcomp(clean_proteins, scale. = TRUE)

# 4. PCA Sonuçlarını Görselleştirme İçin Birleştirme
pca_df <- data.frame(pca_result$x[, 1:2]) %>% 
  bind_cols(final_metadata_protein %>% select(tumor_stage_pathological, sex))

# 5. Grafiği Çizdirme
ggplot(pca_df, aes(x = PC1, y = PC2, color = tumor_stage_pathological)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Proteomik Profil PCA Skor Grafiği",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"),
       color = "Tümör Evresi")

#Heatmap
library(dplyr)
library(pheatmap)

# 1. Protein verilerini seçiyoruz
protein_matrix <- final_metadata_protein %>% 
  select(where(is.numeric)) %>% 
  select(-any_of(c("age", "bmi", "Overall survival, days", "Survival status (1, dead; 0, alive)")))

# 2. Hiç NA içermeyen proteinleri ayıklıyoruz
clean_proteins <- protein_matrix %>% select(where(~all(!is.na(.))))

# 3. Veri tipini kesin olarak sayısal matrise zorlayıp nesneyi yeniden oluşturuyoruz
clean_proteins_numeric <- as.data.frame(lapply(clean_proteins, as.numeric))
rownames(clean_proteins_numeric) <- rownames(clean_proteins)

# 4. Şimdi varyans hesaplamasını güvenle yapabiliriz
protein_vars <- apply(clean_proteins_numeric, 2, var)

# 5. En yüksek varyansa sahip ilk 50 proteini seçiyoruz
top_50_proteins <- names(sort(protein_vars, decreasing = TRUE)[1:50])

# 6. Isı haritası matrisini oluşturuyoruz (Satır=Protein, Sütun=Hasta)
heatmap_matrix <- t(as.matrix(clean_proteins_numeric[, top_50_proteins]))

# 7. Klinik veri için anotasyon tablosu hazırlığı
annotation_col <- data.frame(Stage = final_metadata_protein$tumor_stage_pathological)
rownames(annotation_col) <- colnames(heatmap_matrix)

# 8. Isı Haritasını Çizdirme
pheatmap(heatmap_matrix, 
         scale = "row", 
         annotation_col = annotation_col,
         show_colnames = FALSE,
         main = "En Değişken 50 Proteinin Kümeleme Isı Haritası")