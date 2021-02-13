import matplotlib.pyplot as plt
import numpy as np


# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.





# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    data = np.loadtxt('./../src/Hdiffdata.txt')
    x = data[:,0]
    y = data[:,1]
    #x = 1.0 / x
    plt.plot(x,y)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Number of Steps")
    plt.ylabel(r'$\Delta H$')
    plt.savefig("HdiffPlot.png")
    plt.show()



# See PyCharm help at https://www.jetbrains.com/help/pycharm/
