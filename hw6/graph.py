import matplotlib.pyplot as plt
import numpy as np

n = [100, 1000, 10000, 100000]

m = 0.5
s = 1.5
def gaussian():
    for i in range(len(n)):
        x = s * np.random.randn(n[i]) + m 
        plt.subplot(2,len(n), i+1)
        plt.title('samples : ' + str(n[i]))
        plt.hist(x, 100)


plt.figure(figsize=(15, 5))
gaussian()
plt.show()
