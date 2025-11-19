import pandas as pd

# load table
df = pd.read_csv("out/EOPC_assoc.pancreas_cancer.glm.logistic.hybrid", sep="\t")

# select for additive test
df = df[df['TEST'] == 'ADD']

# select for specific columns
df_cleaned = df[['#CHROM','POS','ID','REF','ALT','A1_FREQ','OR','L95','U95','P']]

# rename columns
df_cleaned = df_cleaned.rename(columns={'#CHROM':'chromosome','POS':'basePairCoordinate','ID':'variantID','REF':'referenceAllele','ALT':'alternateAlleles','A1_FREQ':' alleleFrequency','OD':'oddsRatio','L95':'ci_lower','U95':'ci_upper','P':'pvalue'})

# sort by p-values and take the top 100 variants
ranked_df = df_cleaned.sort_values(by='pvalue',ascending=True).head(100)

# write to csv
ranked_df.to_csv("out/top100candidateVariants.csv",index=False)