s_box = (
    0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
    0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
    0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
    0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
    0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
    0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
    0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
    0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
    0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
    0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
    0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C, 0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
    0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
    0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
    0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
    0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
    0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16,
)


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


import time

# def AES_SBOX(x):
#     t = G256_new_basis(x, g2b)
#     t = G256_inv(t)
#     t = G256_new_basis(t, b2g)
#     t = G256_new_basis(t, A)  # 仿射变换乘
#     return t ^ 0x63

data_A = [0b10001111, 0b11000111, 0b11100011, 0b11110001, 0b11111000, 0b01111100, 0b00111110, 0b00011111]  # 仿射变换矩阵乘
data_g2b = [0b10011000, 0b11110011, 0b11110010, 0b01001000, 0b00001001, 0b10000001, 0b10101001, 0b11111111]
data_b2g = [0b01100100, 0b01111000, 0b01101110, 0b10001100, 0b01101000, 0b00101001, 0b11011110, 0b01100000]

start1 = time.time()
sBox = []
AES = AesBox(data_A, data_g2b, data_b2g)
for i in range(256):
    sBox.append(AES.Aes(i))  # 生成sbox

# sbox = AES.AesTable()
for i, s in enumerate(sBox):
    print(f'%02x' % s, ', ', end='')
    if (i + 1) % 16 == 0:
        print()
end1 = time.time()
###################################################

start2 = time.time()
sBox = []
AES = AesBox(data_A, data_g2b, data_b2g)
for i in range(256):
    sBox.append(s_box[i])  # 生成sbox

# sbox = AES.AesTable()
for i, s in enumerate(sBox):
    print(f'%02x' % s, ', ', end='')
    if (i + 1) % 16 == 0:
        print()
end2 = time.time()
print("SBox執行時間：%f 秒" % (end2 - start2))
print("non-SBox執行時間：%f 秒" % (end1 - start1))
