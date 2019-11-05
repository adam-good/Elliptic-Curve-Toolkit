load('./ECC.sage')

k = 2
field.<a> = GF(2^k)
poly = charpoly(a,'x')
print(f'GF({2**k}) Polynomial: {poly}')

a1=1
a2=0
a3=0
a4=0
a6=1

curve = WeierStrassCurve(a1,a2,a3,a4,a6,field)
#print(list(curve.elems()))
P = WeierStrassPoint(a, a, curve)
Q = -P


print(f'{P*3}')
# x1 = a
# y1 = curve.calculate(x1)[0]
# P = WeierStrassPoint(x1,y1,curve)