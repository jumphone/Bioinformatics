import magic
import numpy as np
import pandas as pd

bmmsc_data = magic.io.load_csv('MATRIX.txt')
libsize = bmmsc_data.sum(axis=1)

bmmsc_data = magic.preprocessing.library_size_normalize(bmmsc_data)
bmmsc_data = np.sqrt(bmmsc_data)
bmmsc_data.head()


magic_op = magic.MAGIC(t=4,k=5)
bmmsc_magic = magic_op.fit_transform(bmmsc_data, genes='all_genes')
bmmsc_magic.head()


bmmsc_magic.to_csv('magic.csv',sep='\t')



