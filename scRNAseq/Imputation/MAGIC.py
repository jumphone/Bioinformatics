import magic
import pandas as pd
magic_operator = magic.MAGIC()
X = pd.read_csv("input.txt",sep='\t')

#X_magic = magic_operator.fit_transform(X, genes=['A2M'])


X_magic = magic_operator.fit_transform(X, genes='all_genes')

X_magic.to_csv('output.csv',sep='\t')
