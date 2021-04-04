import matplotlib.pyplot as plt
import numpy as np
from statistics import mean
def extractPlaquettes(inputfile):
    data_dict = {}
    with open(inputfile, 'r') as datafile:
        lines = datafile.readlines()
        for line in lines:
            l = line.split()
            if l[0] in data_dict.keys():
                data_dict[l[0]].append(l[3])
            else:
                data_dict[l[0]] = [l[3]]

    del data_dict["#g"]


    print(data_dict)

    x =[float(i) for i in data_dict.keys()]
    y = []
    for key in data_dict.keys():
        y.append( mean([ float(j) for j in data_dict[key]]  ))
    x = [4.0 / z**2 for z in x   ]
    z = [np.log(1.0/k) for k in x ]
    y = [ np.log(4*i) for i in y  ]
    plt.plot(x,y)
    plt.plot(x,z)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("1/g^2")
    plt.ylabel(r'$\chi (1,1)$')
    plt.savefig("HdiffPlot.png")
    plt.show()


def jackknife(values):
    knifeset = []
    for i in range(len(values)):
        knifeset.append( [ x for x != values[i]  ] )
    knifeset_averages = []
    for i in range(len(knifeset)):
        knifeset_averages.append(mean(knifeset[i]))

    jackknife_average = mean(knifeset_averages)
    knifeset_var2 = []
        for i in range(len(knifeset_averages))
        knifeset_var2.append( (jackknife_average -  knifeset_averages[i]  )**2    )

    jackknife_variance = ( len(knifeset_var2) - 1  ) * np.sum(knifeset_var2)
    jackknife_error = np.sqrt(jackknife_variance /len(jackknife_variance) )

    return jackknife_average, jackknife_error
