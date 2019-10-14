load('./ECC.sage')

A = 32
B = 37
p = 11

curve = FiniteEllipticCurve(A,B,p)
P,Q = sample(curve.elems(), 2)

elgamal = EC_ElGamalSystem(curve, P)
priv = 32
pub = elgamal.generate_public_key(priv)

plaintext = Q
ciphertext = elgamal.Encrypt(priv, pub, plaintext)
message = elgamal.Decrypt(ciphertext, priv)

print(f'Plaintext:\t {plaintext}')
print(f'Ciphertext:\t {ciphertext}')
print(f'Message:\t {message}')