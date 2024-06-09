from sympy import Matrix
import math

def Kramer(M_Coeff, sv_ch) :
    det_All = M_Coeff.det()
    sol = []
    for i in range(M_Coeff.shape[0]):
        M_Coeff_i = M_Coeff.copy()
        M_Coeff_i[:, i] = sv_ch
        sol.append(M_Coeff_i.det() / det_All)

    return sol

f = open('output.txt', 'w')
p_h2 = 125
p_o2 = 125
p_h2o = 125
p_oh = 125
p_h = 125
p_o = 125
Mt = 125

T = 2900
alpha = 0.6
P = 250
R0 = 1.987

I_h2=20332.5
I_o2=22491.7
I_h2o=-28603.8
I_oh=29985.6
I_h=65027.2
I_o=72578.7

S_h2=48.1685
S_o2=67.6481
S_h2o=68.1315
S_oh=61.0793
S_h=38.6927
S_o=49.9258

Cp_h2=8.81281
Cp_o2=9.48155
Cp_h2o=13.5928
Cp_oh=8.82819
Cp_h=4.96821
Cp_o=4.99041

f.write("Исходные данные по варианту №61 ст.группы 23СЖ01 Ахметова А.Р. " + "\n" +
        "T = " + str(T) + "\n" +
        "p = " + str(P) + "\n" +
        "alpha = " + str(alpha) + "\n" +
        "I_h2 = " + str(I_h2) + "\n" +
        "I_o2 = " + str(I_o2) + "\n" +
        "I_h2o = " + str(I_h2o) + "\n" +
        "I_oh = " + str(I_oh) + "\n" +
        "I_h = " + str(I_h) + "\n" +
        "I_o = " + str(I_o) + "\n" +
        "S_h2 = " + str(S_h2) + "\n" +
        "S_o2 = " + str(S_o2) + "\n" +
        "S_h2o = " + str(S_h2o) + "\n" +
        "S_oh = " + str(S_oh) + "\n" +
        "S_h = " + str(S_h) + "\n" +
        "S_o = " + str(S_o) + "\n" +
        "Cp_h2 = " + str(Cp_h2) + "\n" +
        "Cp_o2 = " + str(Cp_o2) + "\n" +
        "Cp_h2o = " + str(Cp_h2o) + "\n" +
        "Cp_oh = " + str(Cp_oh) + "\n" +
        "Cp_h = " + str(Cp_h) + "\n" +
        "Cp_o = " + str(Cp_o) + "\n"
        )

K_h2 = (2*I_h-I_h2)/R0/T
K_o2 = (2*I_o-I_o2)/R0/T
K_h2o = (2*I_h+I_o-I_h2o)/R0/T
K_oh = (I_o+I_h-I_oh)/R0/T

lg_KP_h2 = (((2*S_h-S_h2)/R0)-((2*I_h)-I_h2)/(R0*T))/2.31
lg_KP_o2 = (((2*S_o-S_o2)/R0)-((2*I_o)-I_o2)/(R0*T))/2.31
lg_KP_h2o = (((2*S_h+S_o-S_h2o)/R0)-((2*I_h+I_o)-I_h2o)/(R0*T))/2.31
lg_KP_oh = (((S_o+S_h-S_oh)/R0)-((I_h+I_o-I_oh))/(R0*T))/2.31

B_h = 2*p_h2+2*p_h2o+p_oh+p_h
B_o = 2*p_o2+p_h2o+p_oh+p_o
p_sum = p_h2 + p_o2 + p_h2o + p_oh + p_h + p_o

del_h2 = math.log10(p_h2)-2*math.log10(p_h)+lg_KP_h2
del_o2 = math.log10(p_o2)-2*math.log10(p_o)+lg_KP_o2
del_h2o = math.log10(p_h2o)-2*math.log10(p_h)-math.log10(p_o)+lg_KP_h2o
del_oh = math.log10(p_oh)-math.log10(p_o)-math.log10(p_h)+lg_KP_oh
del_h = math.log10(B_h)-math.log10(Mt)-math.log10(2)
del_o = math.log10(B_o)-math.log10(Mt)-math.log10(alpha)
del_p = math.log10(p_sum)-math.log10(P)

M_Coeff = Matrix([
    [1, 0, 0, 0, -2, 0, 0],
    [0, 1, 0, 0, 0, -2, 0],
    [0, 0, 1, 0, -2, -1, 0],
    [0, 0, 0, 1, -1, -1, 0],
    [2*p_h2, 0, 2*p_h2o, p_oh, p_h, 0, -B_h],
    [0, 2*p_o2, p_h2o, p_oh, 0, p_o, -B_o],
    [p_h2, p_o2, p_h2o, p_oh, p_h, p_o, 0]
])
sv_ch = Matrix([-del_h2, -del_o2, -del_h2o, -del_oh, -del_h*B_h, -del_o*B_o, -del_p*p_sum])
sol = Kramer(M_Coeff, sv_ch)

p_h2 = 10**(math.log10(p_h2)+sol[0])
p_o2 = 10**(math.log10(p_o2)+sol[1])
p_h2o = 10**(math.log10(p_h2o)+sol[2])
p_oh = 10**(math.log10(p_oh)+sol[3])
p_h = 10**(math.log10(p_h)+sol[4])
p_o = 10**(math.log10(p_o)+sol[5])
Mt = 10**(math.log10(Mt)+sol[6])

for i in range(8):
    B_h = 2*p_h2+2*p_h2o+p_oh+p_h
    B_o = 2*p_o2+p_h2o+p_oh+p_o
    p_sum = p_h2 + p_o2 + p_h2o + p_oh + p_h + p_o

    del_h2 = math.log10(p_h2)-2*math.log10(p_h)+math.log10(K_h2)
    del_o2 = math.log10(p_o2)-2*math.log10(p_o)+math.log10(K_o2)
    del_h2o = math.log10(p_h2o)-2*math.log10(p_h)-math.log10(p_o)+math.log10(K_h2o)
    del_oh = math.log10(p_oh)-math.log10(p_o)-math.log10(p_h)+math.log10(K_oh)
    del_h = math.log10(B_h)-math.log10(Mt)-math.log10(2)
    del_o = math.log10(B_o)-math.log10(Mt)-math.log10(alpha)
    del_p = math.log10(p_sum)-math.log10(P)
    if (i == 6 or i == 7):
        del_h2 = 0
        del_o2 = 0
        del_h2o = 0
        del_oh = 0
        del_h = 0
        del_o = 0
        del_p = 0
    i

    M_Coeff = Matrix([
        [1, 0, 0, 0, -2, 0, 0],
        [0, 1, 0, 0, 0, -2, 0],
        [0, 0, 1, 0, -2, -1, 0],
        [0, 0, 0, 1, -1, -1, 0],
        [2*p_h2, 0, 2*p_h2o, p_oh, p_h, 0, -B_h],
        [0, 2*p_o2, p_h2o, p_oh, 0, p_o, -B_o],
        [p_h2, p_o2, p_h2o, p_oh, p_h, p_o, 0]
    ])
    sv_ch = Matrix([-del_h2, -del_o2, -del_h2o, -del_oh, -del_h*B_h, -del_o*B_o, -del_p*p_sum])
    if (i == 6):
        sv_ch = Matrix([-K_h2, -K_o2, -K_h2o, -K_oh, 0, 0, 0])
    if (i == 7):
        sv_ch = Matrix([0, 0, 0, 0, 0, 0, 250])
    sol = Kramer(M_Coeff, sv_ch)

    if (i < 6):
        p_h2 = 10**(math.log10(p_h2)+sol[0])
        p_o2 = 10**(math.log10(p_o2)+sol[1])
        p_h2o = 10**(math.log10(p_h2o)+sol[2])
        p_oh = 10**(math.log10(p_oh)+sol[3])
        p_h = 10**(math.log10(p_h)+sol[4])
        p_o = 10**(math.log10(p_o)+sol[5])
        Mt = 10**(math.log10(Mt)+sol[6])
    if (i == 6):
        p_h2_8 = p_h2
        p_o2_8 = p_o2
        p_h2o_8 = p_h2o
        p_oh_8 = p_oh
        p_h_8 = p_h
        p_o_8 = p_o
        Mt_8 = Mt
        I_h2_8 = I_h2
        I_o2_8 = I_o2
        I_h2o_8 = I_h2o
        I_oh_8 = I_oh
        I_h_8 = I_h
        I_o_8 = I_o
        sol0_8 = sol[0]
        sol1_8 = sol[1]
        sol2_8 = sol[2]
        sol3_8 = sol[3]
        sol4_8 = sol[4]
        sol5_8 = sol[5]
        sol6_8 = sol[6]
    if (i == 7):
        sol6_9 = sol[6]
print(sol6_9)
        
chis = 1.008*2+16*0.6
R = R0*4184/(chis*Mt/P)
cp = (p_h2*Cp_h2 + p_o2*Cp_o2 + p_h2o*Cp_h2o + p_oh*Cp_oh + p_h*Cp_h + p_o*Cp_o)/Mt/chis*1000*4.184
I = (p_h2*I_h2 + p_o2*I_o2 + p_h2o*I_h2o + p_oh*I_oh + p_h*I_h + p_o*I_o)/Mt/chis*1000*4.184
c_pe = ((p_h2_8*I_h2_8*sol0_8+p_o2_8*I_o2_8*sol1_8+p_h2o_8*I_h2o_8*sol2_8+p_oh_8*I_oh_8*sol3_8+p_h_8*I_h_8*sol4_8+p_o_8*I_o_8*sol5_8)/Mt_8/chis*1000*4.184-I*sol6_8)/T +cp
k_r = c_pe/(c_pe - R*(1-sol6_8)**2/sol6_9)
f.write("Конечные данные: " + "\n" +
        "Количество итераций: 7" + "\n" +
        "p_h2 = " + str(p_h2) + "\n" +
        "p_o2 = " + str(p_o2) + "\n" +
        "p_h2o = " + str(p_h2o) + "\n" +
        "p_oh = " + str(p_oh) + "\n" +
        "p_h = " + str(p_h) + "\n" +
        "p_o = " + str(p_o) + "\n" +
        "I = " + str(I) + "\n" +
        "S = " + str((p_h2*S_h2 + p_o2*S_o2 + p_h2o*S_h2o + p_oh*S_oh + p_h*S_h + p_o*S_o)/Mt/chis*1000*4.184) + "\n" +
        "мю = " + str(chis*Mt/P) + "\n" +
        "R = " + str(R) + "\n" +
        "alpha = " + str(1/T*(1-sol6_8)) + "\n" +
        "betta_t = " + str(1/P*sol6_9) + "\n" +
        "cp = " + str(cp) + "\n" +
        "k_f = " + str(1/(1-R/cp)) + "\n" +
        "a_f = " + str(math.sqrt(1/(1-R/cp)*R*T)) + "\n" +
        "c_pe = " + str(c_pe) + "\n" +
        "k_r = " + str(k_r) + "\n" +
        "a = " + str(math.sqrt(k_r*R*T/sol6_9)))
f.close()
