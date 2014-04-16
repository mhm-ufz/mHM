
# script for checking pet results
#
# written: Matthias Zink, Feb 2014
#
import numpy as np
import ufz
import matplotlib.pyplot as plt


    #  1 : unknown  unknown  interception var  F64       912   1       1   1
    #  2 : unknown  unknown  snowpack    var  F64       912   1       1   1
    #  3 : unknown  unknown  SWC_L01     var  F64       912   1       1   1
    #  4 : unknown  unknown  SWC_L02     var  F64       912   1       1   1
    #  5 : unknown  unknown  SWC_L03     var  F64       912   1       1   1
    #  6 : unknown  unknown  SM_L01      var  F64       912   1       1   1
    #  7 : unknown  unknown  SM_L02      var  F64       912   1       1   1
    #  8 : unknown  unknown  SM_L03      var  F64       912   1       1   1
    #  9 : unknown  unknown  SM_Lall     var  F64       912   1       1   1
    # 10 : unknown  unknown  sealedSTW   var  F64       912   1       1   1
    # 11 : unknown  unknown  unsatSTW    var  F64       912   1       1   1
    # 12 : unknown  unknown  satSTW      var  F64       912   1       1   1
    # 13 : unknown  unknown  aET         var  F64       912   1       1   1
    # 14 : unknown  unknown  Q           var  F64       912   1       1   1
    # 15 : unknown  unknown  QD          var  F64       912   1       1   1
    # 16 : unknown  unknown  QIf         var  F64       912   1       1   1
    # 17 : unknown  unknown  QIs         var  F64       912   1       1   1
    # 18 : unknown  unknown  QB          var  F64       912   1       1   1
    # 19 : unknown  unknown  recharge    var  F64       912   1       1   1


filename = '/home/matze/ufz/source_code/mhm/pet-module/diff_pet_trunk'

data     = ufz.sread(filename, skip=1, strarr=True)

date  = data[:,2]
param = data[:,4].astype('int')
abso   = data[:,11].astype('float')
rela   = data[:,12].astype('float')

# plot abso and rel for aET
maske=(param==-13)
plt.plot(abso[maske])
#plt.plot(rela[maske])
plt.show()
# for iPar in range(np.min(param),0):
#   maske = (param == iPar)
#   plt.plot(rela[maske])

