import matplotlib.pyplot as plt
import pandas as pd
import math

inputData = pd.read_csv('3DTD_Error.csv', sep=",")
for i in range(len(inputData) - 1):
    error = inputData.err[i]
    deltaX = inputData.deltaX[i]
    print((math.log(error) - math.log(inputData.err[i + 1])) / (math.log(1/(1+inputData.deltaX[i])) - math.log(1/(1+inputData.deltaX[i+1]))))
    # Your code here to process error and deltaX

plt.plot(1/(1+inputData.deltaX), inputData.err, 'b--o', label="error")
plt.plot(1/(1+inputData.deltaX), 1/(1+inputData.deltaX) ** 1, 'k--', label="deltaX")
plt.plot(1/(1+inputData.deltaX), 1/(1+inputData.deltaX)** 2, 'k-', label="deltaX^2")
plt.xlabel('delta X')
plt.ylabel('Error')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()


