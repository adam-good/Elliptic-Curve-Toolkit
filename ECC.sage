from random import sample

def square_root_mod(c,p):
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
        self.curve = curve
        if self.curve.field != None and (x,y) != (oo,oo): 
            (x, y) = (x % curve.field, y % curve.field)
        
        (self.x, self.y) = (x, y)
        self.values = (self.x, self.y)

        if not self.is_valid():
            raise Exception("Point not in curve!!")

    def is_valid(self):
        if self.curve == None:
            if self.x == oo and self.y == oo:
                return True
            else:
                return False
        else:
            if self.y in self.curve.calculate(self.x):
                return True
            else:
                return False

    def __eq__(self,point):
        if self.curve != point.curve:
            return False
        x1,y1 = self.values
        x2,y2 = point.values
        if (x1 == x2) and (y1 == y2):
            return True
        else:
            return False

    def __add__(self, point):
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
                slope = (3*x1^2 + A) / (2 * y1)
            elif self != point:
                slope = (y2 - y1) / (x2 - x1)

            x3 = (slope^2 - x1 - x2)
            y3 = (slope * (x1 - x3) - y1)

            return EllipticPoint(x3,y3,self.curve)

    def __sub__(self, point):
        neg = -point
        return self+neg

    def __neg__(self):
        return EllipticPoint(self.x, -self.y, self.curve)

    def __mul__(self, other):
        if isinstance(other, int):
            return EllipticPoint.double_and_add(self.curve, self, other)
        else:
            raise Exception(f'Cannot multiply type {type(other)} by {type(self)}')

    def __rmul__(self, other):
        return self.__mul__(other)

    def __str__(self):
        if self == self.curve.identity():
            return "𝒪"
             #return u'\u1D4AA'
        return "EC(" + str(self.x) + ", " + str(self.y) + ")"

    def __repr__(self):
        if self == self.curve.identity():
            return "𝒪"
             #return u'\u1D4AA'
        return "EC(" + str(self.x) + ", " + str(self.y) + ")"

    def double_and_add(curve, P, n):
        Q = P
        R = curve.identity()
        while n > 0:
            if n % 2 == 1:
                R = Q+R
            Q = Q+Q
            n = n//2
        return R

class EllipticCurve:
    def __init__(self,A,B):
        self.A = A
        self.B = B

        self.discriminate = -16*(4*self.A^3 +26*self.B^2)
        # If the discriminate is 0, the curve is singular. Therefore, this boolean expression works
        self.is_singular = self.discriminate == 0

    def calculate(self, x):
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


