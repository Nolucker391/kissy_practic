import math
import numpy as np

from sympy import Matrix


def strochka():
    print('#' + '-' * 80 + '#', file=file)


# def cramer_solver(coefficients, constants):
#     A = np.array(coefficients)
#     b = np.array(constants)
#     x = np.linalg.solve(A, b)
#     return x
def cramer_solver(coefficients, constants):
    """
    Функция крамера. Принимает в себя аргументы:
    коэффициенты системы уравнений и константы правой части уравнения.
    """
    det_All = coefficients.det()
    sol = []

    for i in range(coefficients.shape[0]):
        M_Coeff_i = coefficients.copy()
        M_Coeff_i[:, i] = constants
        sol.append(M_Coeff_i.det() / det_All)

    return sol


with open('result.txt', 'w') as file:
    """Исходные данные В4"""
    t = 3700
    p = 150  # 15000кПа
    r0 = 1.987
    alpha = 0.5

    o = 0
    b = [2, alpha]

    P_k = {
        'h2': 150,
        'o2': 150,
        'h2o': 150,
        'oh': 150,
        'h': 150,
        'o': 150,
        'Mt': 150
    }

    I_k = {
        'H2': 27525.9,
        'O2': 30243.8,
        'H2O': -17462,
        'OH': 37160.6,
        'H': 69001.8,
        'O': 76589.8
    }
    S_k = {
        'H2': 50.3573,
        'O2': 70.0068,
        'H2O': 71.5214,
        'OH': 63.263,
        'H': 39.9031,
        'O': 51.147
    }
    Cp_k = {
        'H2': 9.19141,
        'O2': 9.88457,
        'H2O': 14.2241,
        'OH': 9.09635,
        'H': 4.96821,
        'O': 5.05084
    }

    aij_ = [[2, 0, 2, 1, 1, 0],
            [0, 2, 1, 1, 0, 1]]

    strochka()
    print('Исходные данные по варианту N4 студента группы 23СЖ02 Губайдуллина Арслана', file=file)
    strochka()

    print('Температура T', '=', t, file=file)
    print('Давление', '=', p, file=file)
    print('Коэффициент избытка окислителя alpha', '=', alpha, file=file)

    for keys, value in I_k.items():
        print(f'I_{keys}', '=', value, 'ккал/кмоль', file=file)

    for keys, value in S_k.items():
        print(f'S_{keys}', '=', value, 'ккал/кмоль', file=file)

    for keys, value in Cp_k.items():
        print(f'Cp_{keys}', '=', value, 'ккал/кмоль', file=file)

    """Решение задачи:"""
    strochka()
    print('Результаты расчета:', file=file)
    strochka()

    K_h2 = (2 * I_k['H'] - I_k['H2']) / r0 / t
    K_o2 = (2 * I_k['O'] - I_k['O2']) / r0 / t
    K_h2o = (2 * I_k['H'] + I_k['O'] - I_k['H2O']) / r0 / t
    K_oh = (I_k['O'] + I_k['H'] - I_k['OH']) / r0 / t

    lg_KP_h2 = (((2 * S_k['H'] - S_k['H2']) / r0) - ((2 * I_k['H']) - I_k['H2']) / (r0 * t)) / 2.31
    lg_KP_o2 = (((2 * S_k['O'] - S_k['O2']) / r0) - ((2 * I_k['O']) - I_k['O2']) / (r0 * t)) / 2.31
    lg_KP_h2o = (((2 * S_k['H'] + S_k['O'] - S_k['H2O']) / r0) - ((2 * I_k["H"] + I_k['O']) - I_k['H2O']) / (
                r0 * t)) / 2.31
    lg_KP_oh = (((S_k['O'] + S_k['H'] - S_k['OH']) / r0) - ((I_k['H'] + I_k['O'] - I_k['OH'])) / (r0 * t)) / 2.31

    B_h = 2 * P_k['h2'] + 2 * P_k['h2o'] + P_k['oh'] + P_k['h']
    B_o = 2 * P_k['o2'] + P_k['h2o'] + P_k['oh'] + P_k['o']
    p_sum = P_k['h2'] + P_k['o2'] + P_k['h2o'] + P_k['oh'] + P_k['h'] + P_k['o']

    del_h2 = math.log10(P_k['h2']) - 2 * math.log10(P_k['h']) + lg_KP_h2
    del_o2 = math.log10(P_k['o2']) - 2 * math.log10(P_k['o']) + lg_KP_o2
    del_h2o = math.log10(P_k['h2o']) - 2 * math.log10(P_k['h']) - math.log10(P_k['o']) + lg_KP_h2o
    del_oh = math.log10(P_k['oh']) - math.log10(P_k['o']) - math.log10(P_k['h']) + lg_KP_oh
    del_h = math.log10(B_h) - math.log10(P_k['Mt']) - math.log10(2)
    del_o = math.log10(B_o) - math.log10(P_k['Mt']) - math.log10(alpha)
    del_p = math.log10(p_sum) - math.log10(p)

    coefficients = Matrix([
        [1, 0, 0, 0, -2, 0, 0],
        [0, 1, 0, 0, 0, -2, 0],
        [0, 0, 1, 0, -2, -1, 0],
        [0, 0, 0, 1, -1, -1, 0],
        [2 * P_k['h2'], 0, 2 * P_k['h2o'], P_k['oh'], P_k['h'], 0, -B_h],
        [0, 2 * P_k['o2'], P_k['h2o'], P_k['oh'], 0, P_k['o'], -B_o],
        [P_k['h2'], P_k['o2'], P_k['h2o'], P_k['oh'], P_k['h'], P_k['o'], 0]
    ])

    sv_ch = Matrix([-del_h2, -del_o2, -del_h2o, -del_oh, -del_h * B_h, -del_o * B_o, -del_p * p_sum])
    sol = cramer_solver(coefficients, sv_ch)

    p_h2 = 10 ** (math.log10(P_k['h2']) + sol[0])
    p_o2 = 10 ** (math.log10(P_k['o2']) + sol[1])
    p_h2o = 10 ** (math.log10(P_k['h2o']) + sol[2])
    p_oh = 10 ** (math.log10(P_k['oh']) + sol[3])
    p_h = 10 ** (math.log10(P_k['h']) + sol[4])
    p_o = 10 ** (math.log10(P_k['o']) + sol[5])
    Mt = 10 ** (math.log10(P_k['Mt']) + sol[6])

    for i in range(8):
        B_h = 2 * p_h2 + 2 * p_h2o + p_oh + p_h
        B_o = 2 * p_o2 + p_h2o + p_oh + p_o
        p_sum = p_h2 + p_o2 + p_h2o + p_oh + p_h + p_o

        del_h2 = math.log10(p_h2) - 2 * math.log10(p_h) + math.log10(K_h2)
        del_o2 = math.log10(p_o2) - 2 * math.log10(p_o) + math.log10(K_o2)
        del_h2o = math.log10(p_h2o) - 2 * math.log10(p_h) - math.log10(p_o) + math.log10(K_h2o)
        del_oh = math.log10(p_oh) - math.log10(p_o) - math.log10(p_h) + math.log10(K_oh)
        del_h = math.log10(B_h) - math.log10(Mt) - math.log10(2)
        del_o = math.log10(B_o) - math.log10(Mt) - math.log10(alpha)
        del_p = math.log10(p_sum) - math.log10(p)

        if (i == 6 or i == 7):
            del_h2 = 0
            del_o2 = 0
            del_h2o = 0
            del_oh = 0
            del_h = 0
            del_o = 0
            del_p = 0


        coefficients_second = Matrix([
            [1, 0, 0, 0, -2, 0, 0],
            [0, 1, 0, 0, 0, -2, 0],
            [0, 0, 1, 0, -2, -1, 0],
            [0, 0, 0, 1, -1, -1, 0],
            [2 * p_h2, 0, 2 * p_h2o, p_oh, p_h, 0, -B_h],
            [0, 2 * p_o2, p_h2o, p_oh, 0, p_o, -B_o],
            [p_h2, p_o2, p_h2o, p_oh, p_h, p_o, 0]
        ])
        sv_ch = Matrix([-del_h2, -del_o2, -del_h2o, -del_oh, -del_h * B_h, -del_o * B_o, -del_p * p_sum])

        if (i == 6):
            sv_ch = Matrix([-K_h2, -K_o2, -K_h2o, -K_oh, 0, 0, 0])
        if (i == 7):
            sv_ch = Matrix([0, 0, 0, 0, 0, 0, 250])
        sol = cramer_solver(coefficients_second, sv_ch)

        if (i < 6):
            p_h2 = 10 ** (math.log10(p_h2) + sol[0])
            p_o2 = 10 ** (math.log10(p_o2) + sol[1])
            p_h2o = 10 ** (math.log10(p_h2o) + sol[2])
            p_oh = 10 ** (math.log10(p_oh) + sol[3])
            p_h = 10 ** (math.log10(p_h) + sol[4])
            p_o = 10 ** (math.log10(p_o) + sol[5])
            Mt = 10 ** (math.log10(Mt) + sol[6])
        if (i == 6):
            p_h2_8 = p_h2
            p_o2_8 = p_o2
            p_h2o_8 = p_h2o
            p_oh_8 = p_oh
            p_h_8 = p_h
            p_o_8 = p_o
            Mt_8 = Mt
            I_h2_8 = I_k['H2']
            I_o2_8 = I_k['O2']
            I_h2o_8 = I_k['H2O']
            I_oh_8 = I_k['OH']
            I_h_8 = I_k['H']
            I_o_8 = I_k['O']
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

    chis = 1.008 * 2 + 16 * 0.6
    R = r0 * 4184 / (chis * Mt / p)
    cp = (
                 p_h2 * Cp_k['H2'] + p_o2 * Cp_k['O2'] + p_h2o * Cp_k['H2O'] + p_oh * Cp_k['OH'] + p_h * Cp_k[
             'H'] + p_o * Cp_k['O']) / Mt / chis * 1000 * 4.184
    I = (p_h2 * I_k['H2'] + p_o2 * I_k['O2'] + p_h2o * I_k['H2O'] + p_oh * I_k['OH'] + p_h * I_k['H'] + p_o * I_k[
        'O']) / Mt / chis * 1000 * 4.184
    c_pe = ((
                    p_h2_8 * I_h2_8 * sol0_8 + p_o2_8 * I_o2_8 * sol1_8 + p_h2o_8 * I_h2o_8 * sol2_8 + p_oh_8 * I_oh_8 * sol3_8 + p_h_8 * I_h_8 * sol4_8 + p_o_8 * I_o_8 * sol5_8) / Mt_8 / chis * 1000 * 4.184 - I * sol6_8) / t + cp
    k_r = c_pe / (c_pe - R * (1 - sol6_8) ** 2 / sol6_9)

    file.write(
        f"Конечные данные: \n" 
        f"p_h2 = {str(p_h2)}\n" 
        f"p_o2 = {str(p_o2)}\n" 
        f"p_h2o = {str(p_h2o)}\n" 
        f"p_oh = {str(p_oh)}\n" 
        f"p_h = {str(p_h)}\n" 
        f"p_o = {str(p_o)}\n"
        f"I = {str(I)}\n" 
        f"S = {str(
            (p_h2 * S_k['H2'] + p_o2 * S_k['O2'] + p_h2o * S_k['H2O'] + p_oh * 
             S_k['OH'] + p_h * S_k['H'] + p_o * S_k['O']) / Mt / chis * 1000 * 4.184)}\n"
        f"мю = {str(chis * Mt / p)}\n"
        f"R = {str(R)}\n" 
        f"alpha = {str(1 / t * (1 - sol6_8))}\n" 
        f"betta_t = {str(1 / p * sol6_9)}\n" 
        f"cp = {str(cp)}\n"
        f"k_f = {str(1 / (1 - R / cp))}\n"
        f"a_f = {str(math.sqrt(1 / (1 - R / cp) * R * t))}\n"
        f"c_pe = {str(c_pe)}\n"
        f"k_r = {str(k_r)}\n" 
        f"a = {str(math.sqrt(k_r * R * t / sol6_9))}"
    )
