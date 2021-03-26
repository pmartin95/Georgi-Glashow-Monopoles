import matplotlib.pyplot as plt
import numpy as np


# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.





# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    data_dict = {}

    with open("HMCsim1.txt", 'r') as datafile:
        lines = datafile.readlines()
        for line in lines:
            l = line.split()
            if l[0] in data_dict.keys():
                data_dict[l[0]].append(l[3])
            else:
                data_dict[l[0]] = [l[3]]

    del data_dict["#g"]
    print(data_dict)
    #x = 1.0 / x
    # plt.plot(x,y)
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.xlabel("Number of Steps")
    # plt.ylabel(r'$\Delta H$')
    # plt.savefig("HdiffPlot.png")
    # plt.show()



# See PyCharm help at https://www.jetbrains.com/help/pycharm/
