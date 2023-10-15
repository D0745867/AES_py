class GF:
    @staticmethod
    def G4_mul(x, y):
        """
        GF(2^2)的乘法运算，正规基{W^2, W}
        """
        a = (x & 0x02) >> 1
        b = x & 0x01
        c = (y & 0x02) >> 1
        d = y & 0x01
        e = (a ^ b) & (c ^ d)
        return (((a & c) ^ e) << 1) | ((b & d) ^ e)

    @staticmethod
    def G4_mul_N(x):
        """
        GF(2^2)的乘N操作，N = W^2
        """
        a = (x & 0x02) >> 1
        b = x & 0x01
        p = b
        q = a ^ b
        return (p << 1) | q

    @staticmethod
    def G4_mul_N2(x):
        """
        GF(2^2)的乘N^2操作，N = W^2
        """
        a = (x & 0x02) >> 1
        b = x & 0x01
        return ((a ^ b) << 1) | a

    @staticmethod
    def G4_inv(x):
        """
        GF(2^2)的求逆操作，该操作和GF(2^2)的平方操作等价
        """
        a = (x & 0x02) >> 1
        b = x & 0x01
        return (b << 1) | a

    @classmethod
    def G16_mul(cls, x, y):
        """
        GF(2^4)的乘法操作，正规基{Z^4, Z}
        """
        a = (x & 0xc) >> 2
        b = x & 0x03
        c = (y & 0xc) >> 2
        d = y & 0x03
        e = cls.G4_mul(a ^ b, c ^ d)
        e = cls.G4_mul_N(e)
        p = cls.G4_mul(a, c) ^ e
        q = cls.G4_mul(b, d) ^ e
        return (p << 2) | q

    @classmethod
    def G16_sq_mul_u(cls, x):
        """
        GF(2^4)的平方后乘u操作, u = N^2Z, N = W^2
        """
        a = (x & 0xc) >> 2
        b = x & 0x03
        p = cls.G4_inv(a ^ b)  # G4平方和求逆等价
        q = cls.G4_mul_N2(cls.G4_inv(b))
        return (p << 2) | q

    @classmethod
    def G16_inv(cls, x):
        """
        GF(2^4)的求逆操作
        """
        a = (x & 0xc) >> 2
        b = x & 0x03
        c = cls.G4_mul_N(cls.G4_inv(a ^ b))
        d = cls.G4_mul(a, b)
        e = cls.G4_inv(c ^ d)
        p = cls.G4_mul(e, b)
        q = cls.G4_mul(e, a)
        return (p << 2) | q

    @classmethod
    def G256_inv(cls, x):
        """
        GF(2^8)的求逆操作
        """
        a = (x & 0xF0) >> 4
        b = x & 0x0F
        c = cls.G16_sq_mul_u(a ^ b)
        d = cls.G16_mul(a, b)
        e = cls.G16_inv(c ^ d)
        p = cls.G16_mul(e, b)
        q = cls.G16_mul(e, a)
        return (p << 4) | q

    @staticmethod
    def G256_new_basis(x, b):
        """
        x在新基b下的表示
        """
        y = 0
        for i in range(8):
            if x & (1 << (7 - i)):
                y ^= b[i]
        return y


class AesBox:
    def __init__(self, A, g2b, b2g):
        self.A = A
        self.g2b = g2b
        self.b2g = b2g

    def Aes(self, x):
        gf = GF()
        t = gf.G256_new_basis(x, self.g2b)
        t = gf.G256_inv(t)
        t = gf.G256_new_basis(t, self.b2g)
        t = gf.G256_new_basis(t, self.A)
        return t ^ 0x63

    def AesTable(self):
        aes = self
        sBox = []
        for i in range(256):
            sBox.append(aes.Aes(x=i))  # 生成sBox
        return sBox


# def AES_SBOX(x):
#     t = G256_new_basis(x, g2b)
#     t = G256_inv(t)
#     t = G256_new_basis(t, b2g)
#     t = G256_new_basis(t, A)  # 仿射变换乘
#     return t ^ 0x63

data_A = [0b10001111, 0b11000111, 0b11100011, 0b11110001, 0b11111000, 0b01111100, 0b00111110, 0b00011111]  # 仿射变换矩阵乘
data_g2b = [0b10011000, 0b11110011, 0b11110010, 0b01001000, 0b00001001, 0b10000001, 0b10101001, 0b11111111]
data_b2g = [0b01100100, 0b01111000, 0b01101110, 0b10001100, 0b01101000, 0b00101001, 0b11011110, 0b01100000]

sBox = []
AES = AesBox(data_A, data_g2b, data_b2g)
for i in range(256):
    sBox.append(AES.Aes(i))  # 生成sbox

# sbox = AES.AesTable()
for i, s in enumerate(sBox):
    print(f'%02x' % s, ', ', end='')
    if (i + 1) % 16 == 0:
        print()
