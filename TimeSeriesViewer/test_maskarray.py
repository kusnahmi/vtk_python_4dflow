import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt

a = np.arange(20).astype(float).reshape(2, 5, 2, 1)
b = ma.greater(a, 6)
c = ~b

print(a)
print(b)
print(c)
print(5 * c)

print(a * b)
print(a * c)
print(a * b * 2 + a * c / 2)
print(a * b * 2 + a * ~b / 2)

#
# >>> a = np.arange(4)
# >>> a
# array([0, 1, 2, 3])
# >>> ma.masked_where(a <= 2, a)
# masked_array(data=[--, --, --, 3],
#              mask=[ True,  True,  True, False],
#        fill_value=999999)

print('===========================================')

a = np.arange(2, 98)
b = a.reshape(3, 4, 2, 2, 2)
c = np.pad(b, ((1, 1), (1, 1), (1, 1), (0, 0), (0, 0)), 'constant', constant_values=99)
print(a)
print(b)
print(b.shape)
print(c)
print(c.shape)

print('===========================================')

a = np.arange(9).reshape(3, 3)
b = np.zeros([3, 3])
b[:, 2] = 1
c = a * b
print(a)
print(b)
print(c)
print(a[0:-1,:])
print(a[0:-2,:])
print(a[:,:-1])
print(a[:,:-2])
print(a[:,-1])

plt.figure()
plt.imshow(a)
plt.matshow(a)
plt.show()

print('===========================================')
a = np.arange(20).astype(float).reshape(2, 5, 2)
x = 1
m = np.mean(a, x)
# m_repeat = np.repeat(m, a.shape[x], axis=x).reshape(a.shape)
m_repeat = np.repeat(m, a.shape[x]).reshape([2,2,5]).transpose([0,2,1])
print(a)
print(a.shape)
print(m)
print(m.shape)
print(m_repeat)
print(m_repeat.shape)

print('===========================================')
a = np.arange(20).astype(float).reshape(2, 5, 2)
x = 1
m = np.mean(a, x, keepdims=True)
print(a)
print(a.shape)
print(m)
print(m.shape)

az=a-m
print(az)
# m_repeat = np.repeat(m, a.shape[x], axis=x).reshape(a.shape)
# m_repeat = np.repeat(m, a.shape[x]).reshape([2,2,5]).transpose([0,2,1])
# print(m_repeat)
# print(m_repeat.shape)
