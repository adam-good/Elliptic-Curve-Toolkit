'''
    Filename:       ECC.sage
    Author:         Adam Good
    Description:    Provides a basic framework for Elliptic Curve Operations with relation to Cryptography.
'''

def square_root_mod(c,p):
    """Solves the equation x^2 = x mod p for x
    
    Arguments:
        c {int} -- The value to take the square root of
        p {int} -- The Modulus
    
    Returns:
        (int, int) -- A tuple containing the positive and negative roots of the equation (None if they don't exist)
    """
    L = legendre_symbol(c,p)
    if L == -1:
        return None
    elif L == 0:
        return (0, 0)

    if p % 4 == 3:
        exponent = (p+1)//4
        root = power_mod(c,exponent,p)
        neg_root = -root % p
        return (root, neg_root)
    else:
        test_exponent = (p-1)//4
        if power_mod(c, test_exponent, p) == (1 % p):
            root = power_mod(c,(p+3)//8,p)
            neg_root = -root % p
            return (root, neg_root)
        elif power_mod(c, test_exponent, p) == (-1 % p):
            root = ( power_mod(2,(p-1)//4,p) * power_mod(c, (p+3)//8, p) ) % p
            neg_root = -root % p
            return (root, neg_root)

class EllipticPoint:
    def __init__(self, x, y, curve):
        """Initialize an EllipticPoint
        
        Arguments:
            x {int} -- The x value of the point on an Elliptic Curve
            y {int} -- The y value of the point on an Elliptic Curve (X^3 + AX^2 + B mod p)
            curve {EllipticCurve} -- The Elliptic Curve that the point exists on/in
        
        Raises:
            Exception: The point does not exists in/on the curve.
        """
        self.curve = curve
        if self.curve.field != None and (x,y) != (oo,oo): 
            (x, y) = (x % curve.field, y % curve.field)
        
        (self.x, self.y) = (x, y)
        self.values = (self.x, self.y)

        if not self.is_valid():
            raise Exception("Point not in curve!!")

    def is_on_curve(self, curve):
        """Checks if this point is on a given elliptic curve
        
        Arguments:
            curve {EllipticCurve} -- An elliptic curve to be checked
        
        Returns:
            bool -- True if point is on curve, false otherwise
        """
        if self.y in curve.calculate(self.x):
            return True
        else:
            return False

    def is_valid(self):
        """Checks if the point is a valid point on the given curve
        
        Returns:
            bool -- True if the point exists on the curve; False otherwise
        """
        if self.curve == None:
            if self.x == oo and self.y == oo:
                return True
            else:
                return False
        elif self.curve.calculate(self.x) == None:
            return False
        else:
            if self.y in self.curve.calculate(self.x):
                return True
            else:
                return False

    def get_order(self):
        current = self
        for q in range(self.curve.field + 1):
            if current == self.curve.identity():
                return q+1
            current = current + self

    def __eq__(self,point):
        """Overload the Equality Operator 
        
        Arguments:
            point {EllipticPoint} -- The Elliptic Point that this Point is being compared to
        
        Returns:
            bool -- True if points are equal; False otherwise
        """
        if self.curve != point.curve:
            return False
        x1,y1 = self.values
        x2,y2 = point.values
        if (x1 == x2) and (y1 == y2):
            return True
        else:
            return False

    def __add__(self, point):
        """An Overload of the Addition Operator (self + point)
        
        Arguments:
            point {EllipticPoint} -- The point that this point is being added to
        
        Raises:
            Exception: The points do not exist on the same curve
            Exception: The curve is singular
        
        Returns:
            EllipticPoint -- The solution from adding points on an Elliptic Curve
        """
        if self.curve != point.curve:
            raise Exception("Points not in same curve")

        if self.curve.is_singular:
            raise Exception("Cannot add points on a singular curve!")

        curve = self.curve

        if self == curve.identity():
            return point
        elif point == curve.identity():
            return self
        else:
            x1,y1 = self.values
            x2,y2 = point.values

            A = self.curve.A
            B = self.curve.B

            if self == -point:
                return curve.identity()
            elif self == point:
                # TODO: This will break on an EC over a Finite Field where the inverse of the numerators DNE
                slope = (3*x1^2 + A) / (2 * y1)
            elif self != point:
                # TODO: This will break on an EC over a Finite Field where the inverse of the numerators DNE
                slope = (y2 - y1) / (x2 - x1)

            x3 = (slope^2 - x1 - x2)
            y3 = (slope * (x1 - x3) - y1)

            return EllipticPoint(x3,y3,self.curve)

    def __sub__(self, point):
        """An overload of the subtraction operator (self - point)
        
        Arguments:
            point {EllipticPoint} -- The Elliptic Point that will be subtracted from this point
        
        Returns:
            EllipticPoint -- The solution of subtracting points on an Elliptic Curve
        """
        neg = -point
        return self+neg

    def __neg__(self):
        """An overload of the arithmetic negation operator (-self)
        
        Returns:
            EllipticPoint -- Returns the negation of self (reflection of self across the X axis)
        """
        return EllipticPoint(self.x, -self.y, self.curve)

    def __mul__(self, other):
        """An overload of the multiplication operator (self * other)
        
        Arguments:
            other {int} -- The value to be self is to be multiplied by.

        Raises:
            Exception: Cannot multiply non-int type by EllipticPoint
        
        Returns:
            EllipticPoint -- The result of multiplying an Elliptic Point by a Scalar
        """
        if isinstance(other, int):
            return EllipticPoint.double_and_add(self.curve, self, other)
        elif isinstance(other, sage.rings.integer.Integer):
            return self.__mul__(int(other))
        else:
            raise Exception(f'Cannot multiply type {type(other)} by {type(self)}')

    def __rmul__(self, other):
        """An overload of the multiplication operator [reversed order] (other * self)
        
        Arguments:
            other {int} -- The value that self will be multiplied by
        
        Returns:
            EllipticPoint -- The result of multiplying an EllipticPoint by a Scalar value
        """
        return self.__mul__(other)

    def __str__(self):
        """Overload of the string representation of an EllipticPoint
        
        Returns:
            str -- A string reptresentation of an Elliptic Point
        """
        if self == self.curve.identity():
            return "𝒪"
        return "EC(" + str(self.x) + ", " + str(self.y) + ")"

    def __repr__(self):
        """Overload of the string representation of an EllipticPoint
        
        Returns:
            str -- A string reptresentation of an Elliptic Point
        """
        if self == self.curve.identity():
            return "𝒪"
             #return u'\u1D4AA'
        return "EC(" + str(self.x) + ", " + str(self.y) + ")"

    def double_and_add(curve, P, n):
        """An implementation of the Double and Add algorithm used for efficient repeated Elliptic Point addition
        
        Arguments:
            curve {EllipticCurve} -- An Elliptic Curve that the Elliptic Point exists on
            P {EllipticPoint} -- An Elliptic Point to be added repeatadly
            n {int} -- The number of times to add P
        
        Returns:
            EllipticPoint -- The result of the repeated addition
        """
        Q = P
        R = curve.identity()
        while n > 0:
            if n % 2 == 1:
                R = Q+R
            Q = Q+Q
            n = n//2
        return R

    def ECDLP_collision_attack(curve, P, Q, list_length):
        raise Exception("This is super broke, don't use it")
        p = curve.field

        js = [randint(1,p) for i in range(list_length)]
        L1 = { (j*P).values:j for j in js }

        solution = None
        for _ in range(list_length):
            k = randint(1,p)
            l2 = (Q + k*P).values
            if l2 in L1.keys():
                j = L1[l2]
                solution = (j-k) % p
                break

        return solution

    def ECDLP_bruteforce_attack(curve, P, Q, tries):
        for n in range(1, tries):
            if Q == n*P:
                return n
            
class EllipticCurve:
    def __init__(self,A,B):
        """Initialize an Elliptic Curve defined by Y^2 = X^3 + AX + B
        
        Arguments:
            A {int} -- The coefficient of X^1 in the Elliptic Curve Equation
            B {int} -- The coefficient of X^0 in the Elliptic Curve Equation
        """
        self.A = A
        self.B = B

        self.discriminate = -16*(4*self.A^3 +26*self.B^2)
        # If the discriminate is 0, the curve is singular. Therefore, this boolean expression works
        self.is_singular = self.discriminate == 0

    def calculate(self, x):
        """Calculate the Y 
        
        Arguments:
            x {[type]} -- [description]
        
        Returns:
            [type] -- [description]
        """
        if x == oo:
            return oo,oo
        else:
            return sqrt(x^3 + self.A*x + self.B)

    def identity(self):
        return EllipticPoint(oo, oo, self)

    def __eq__(self, curve):
        if curve == None:
            return False
        
        if self.A == curve.A and self.B == curve.B:
            return True
        else:
            return False

class FiniteEllipticCurve(EllipticCurve):
    def __init__(self,A,B,p):
        super(FiniteEllipticCurve,self).__init__(A,B)
        self.A = A
        self.B = B
        self.field = p

    def calculate(self, x):
        if x == oo:
            return oo,oo
        else:
            return square_root_mod(x^3 + self.A*x + self.B, self.field)

    def get_bit(self, x, y):
        if y not in self.calculate(x):
            raise Exception("This y value does not correspond to the given x value")

        if 0 <= y and y < self.field//2:
            return 0
        elif self.field//2 <= y < self.field:
            return 1

    def elems(self):
        curve = [self.identity()]

        for x in range(self.field):
            LHS = self.calculate(x)
            if LHS != None:
                pos_root, neg_root = LHS
                if pos_root == neg_root:
                    point = EllipticPoint(x, pos_root, self)
                    curve.append( point )
                else:
                    pos_point = EllipticPoint(x, pos_root, self)
                    neg_point = EllipticPoint(x, neg_root, self)
                    curve.append( pos_point )
                    curve.append( neg_point )
        return curve

    def __eq__(self, curve):
        if super(FiniteEllipticCurve, self).__eq__(curve) == True and self.field == curve.field:
            return True
        else:
            return False

class EC_MVElGamalSystem():
    def __init__(self, curve, point):
        self.curve = curve
        self.point = point
        self.p = self.curve.field

        if self.curve.is_singular:
            raise Exception("Can't do ElGamal Encryption on Singular Curve")
        
        if not self.point.is_on_curve(self.curve):
            raise Exception("Point must exist on Curve for ElGamal")

    def generate_public_key(self, secret_key):
        return secret_key * self.point

    def encrypt(self, plaintext, public_key, k=None):
        m1, m2 = plaintext
        if k == None:
            k = randint(1, self.p)
        R = k * self.point

        S = k * public_key
        xs, ys = S.values

        c1 = xs*m1 % self.p
        c2 = ys*m2 % self.p

        return (R,c1,c2)

    def decrypt(self, secret_key, ciphertext):
        R, c1, c2 = ciphertext

        T = secret_key * R
        xt,yt = T.values

        m1 = (inverse_mod(xt, self.p) * c1) % self.p
        m2 = (inverse_mod(yt, self.p) * c2) % self.p

        return (m1,m2)

class EC_ElGamalSystem():
    def __init__(self, curve, point):
        self.curve = curve
        self.point = point

        if self.curve.is_singular:
            raise Exception("Can't do ElGamal Encryption on Singular Curve")
        
        if not self.point.is_on_curve(self.curve):
            raise Exception("Point must exist on Curve for ElGamal")

    def generate_public_key(self, private_key):
        return private_key * self.point
        
    def Encrypt(self, private_key, public_key, message):
        c1 = private_key * self.point
        c2 = message + private_key * public_key
        return (c1, c2)

    def Decrypt(self, ciphertext, private_key):
        c1,c2 = ciphertext
        return c2 - private_key * c1

class EC_DiffieHellmanSystem():
    def __init__(self, curve, point):
        self.curve = curve
        self.point = point

    def generate_public_key(self, private_key):
        return private_key * self.point

    def compute_secret(self, private_key, public_key):
        return private_key * public_key

class ECDSA():
    def __init__(self, curve, point):
        self.curve = curve
        self.point = point
        self.order = self.point.get_order()

        if not is_prime(self.order):
            raise Exception("The Order of the public point must be prime")

    def verification_key(self, signing_key):
        return signing_key * self.point

    def sign(self, signing_key, ephemeral_key, document):
        document = document % self.order
        ephemeral_key = ephemeral_key % self.order

        S1 = (ephemeral_key * self.point).values[0] % self.order
        S2 = ((document + signing_key * S1) * inverse_mod(ephemeral_key, self.order)) % self.order

        return (S1, S2)

    def verify(self, verification_key, document, S1, S2):
        v1 = (document * inverse_mod(S2, self.order)) % self.order
        v2 = (S1 * inverse_mod(S2, self.order)) % self.order

        x,y = (v1 * self.point + v2*verification_key).values
        x,y = x%self.order, y%self.order

        return x == S1

        
