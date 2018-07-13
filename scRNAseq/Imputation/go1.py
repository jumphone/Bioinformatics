import magic
import numpy as np
import pandas as pd

bmmsc_data = magic.io.load_csv('MATRIX.txt')
libsize = bmmsc_data.sum(axis=1)

#plt.hist(libsize, bins=50)
#plt.axvline(1000, c='r')
#plt.show()



