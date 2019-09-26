load('./ECC.sage')

A = 32
B = 37
p = 11

curve = FiniteEllipticCurve(A,B,p)
P,Q = sample(curve.elems(), 2)

n = int(123)
s = curve.identity()
for i in range(n):
    s = s + Q
    # print(s)
print(f'{s}')
print(f'{Q*n}')
print(f'{n*Q}')
print(curve.identity())
print(curve.elems())