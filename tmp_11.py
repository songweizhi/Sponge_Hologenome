# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.spatial.distance import squareform
# from scipy.cluster.hierarchy import dendrogram, linkage
#
#
# mat = np.array([[0.0, 2.0, 0.1], [2.0, 0.0, 2.0], [0.1, 2.0, 0.0]])
#
# print(mat)
#
#
# dists = squareform(mat)
# linkage_matrix = linkage(dists, "single")
# dendrogram(linkage_matrix, labels=["0", "1", "2"])
# plt.title("test")
# plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# data = pd.read_csv('Wholesale customers data.csv')
# data.head()

df_A = pd.read_csv('/Users/songweizhi/Desktop/TEST.MATRIX', sep='\t', header=1)
print(df_A)


df_txt = '/Users/songweizhi/Desktop/Mantel_similarity_matrix.txt'



