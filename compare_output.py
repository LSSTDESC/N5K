import numpy as np
import matplotlib.pyplot as plt



m_clgg = np.load("tests/matter_clgg.npz")
c_clgg = np.load("tests/test_clgg.npz")

m_ls = m_clgg['ls']
m_cls = m_clgg['cls']
c_ls = c_clgg['ls']
c_cls = c_clgg['cls']
plt.loglog(m_ls,m_cls[0]/5.64101940381949,label="MATTER")
plt.loglog(c_ls,c_cls[0],label="REF")
plt.legend()
plt.show()
