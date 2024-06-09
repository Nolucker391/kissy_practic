import math
import sympy as sp
from sympy import Matrix


def strochka():
    print('#' + '-' * 80 + '#', file=file)


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


with open('result_prog.txt', 'w') as file:
    o = 0
    p = 250
    t = 2700
    ro = 1.987
    k = 0.8
    b = [2, k]
    p_ = [p for i in range(6)]
    psum = sum(p_)
    Mt = p
    I_k = [18579.6, 20605.7, -31302.2, 28228.5, 64033.6, 71581]
    S_k = [47.5422, 66.9743, 67.1675, 60.4515, 38.3377, 49.5693]
    Cp_k = [8.71414, 9.37875, 13.3858, 8.74102, 4.96821, 4.98706]
    aij_ = [[2, 0, 2, 1, 1, 0],
            [0, 2, 1, 1, 0, 1]]

    strochka()
    print('Исходные данные по варианту N4 студента группы 23СЖ02 Губайдуллина Арслана', file=file)
    strochka()

    print('Температура T', '=', t, file=file)
    print('Давление', '=', p, file=file)
    print('Коэффициент избытка окислителя alpha', '=', k, file=file)

    ko = ['H2', 'O2', 'H2O', 'OH', 'H', 'O']
    for i in range(6):
        print('I_', ko[i], '=', I_k[i], 'ккал/кмоль', file=file)
    for i in range(6):
        print('S_', ko[i], '=', S_k[i], 'ккал/кмоль/К', file=file)
    for i in range(6):
        print('Cp_', ko[i], '=', Cp_k[i], 'ккал/кмоль/К', file=file)
    strochka()

    print('Результаты расчета:', file=file)
    strochka()

    k_k = [0 for i in range(4)]
    for j in range(4):
        sum1 = 0
        sum2 = 0
        for i in range(2):
            sum1 += aij_[i][j] * S_k[i + 4]
            sum2 += aij_[i][j] * I_k[i + 4]
            k_k[j] = ((sum1 - S_k[j]) / ro - (sum2 - I_k[j]) / ro / t) / 2.31

    k_kk = [0 for i in range(4)]
    for j in range(4):
        sum3 = 0
        for i in range(2):
            sum3 += aij_[i][j] * I_k[i + 4]
            k_kk[j] = ((sum3 - I_k[j]) / ro / t)

    while True:
        p_sum = 0
        for i in range(6):
            p_sum += p_[i]

        B_ = [0 for i in range(2)]
        for i in range(2):
            for j in range(6):
                B_[i] += aij_[i][j] * p_[j]

        s_k = [0 for i in range(6)]
        for j in range(4):
            sum4 = 0
            for i in range(2):
                sum4 += aij_[i][j] * math.log10(p_[i + 4])
                s_k[j] = math.log10(p_[j]) - sum4 + k_k[j]
                s_k[i + 4] = math.log10(B_[i]) - math.log10(Mt) - math.log10(b[i])
        s_kp = math.log10(p_sum) - math.log10(p)
        coefficients = [[1, 0, 0, 0, -2, 0, 0],
                        [0, 1, 0, 0, 0, -2, 0],
                        [0, 0, 1, 0, -2, -1, 0],
                        [0, 0, 0, 1, -1, -1, 0],
                        [2 * p_[0], 0, 2 * p_[2], p_[3], p_[4], 0, -(B_[0])],
                        [0, 2 * p_[1], p_[2], p_[3], 0, p_[5], -(B_[1])],
                        [p_[0], p_[1], p_[2], p_[3], p_[4], p_[5], 0]]
        # coefficients = Matrix([
        #     [1, 0, 0, 0, -2, 0, 0],
        #     [0, 1, 0, 0, 0, -2, 0],
        #     [0, 0, 1, 0, -2, -1, 0],
        #     [0, 0, 0, 1, -1, -1, 0],
        #     [2 * p_[0], 0, 2 * p_[2], p_[3], p_[4], 0, -(B_[0])],
        #     [0, 2 * p_[1], p_[2], p_[3], 0, p_[5], -(B_[1])],
        #     [p_[0], p_[1], p_[2], p_[3], p_[4], p_[5], 0]
        # ])
        constants = Matrix([-s_k[0], -s_k[1], -s_k[2], -s_k[3], -s_k[4] * (B_[0]), -s_k[5] * (B_[1]), -s_kp * p_sum])

        solution = cramer_solver(Matrix(coefficients), constants)

        for i in range(6):
            p_[i] = 10 ** (math.log10(p_[i]) + solution[i])
        Mt = 10 ** (math.log10(Mt) + solution[6])

        svch = sum(constants) - (-s_kp * p_sum)
        # object_matrix_list = coefficients.tolist()
        # print(object_matrix_list)
        o += 1
        print('приближение', o, file=file)

        for i in range(len(coefficients)):
            for j in range(len(coefficients[i])):
                print("{:22}".format(coefficients[i][j]), end=" ", file=file)
            if i < len(solution):
                print("{:.10f}".format(constants[i]), file=file)
        if abs(svch) < 0.000001:
            break
        # if coefficients:
        #     num_columns = coefficients.shape[1]
        #     for i in range(coefficients.rows):
        #         for j in range(num_columns):
        #             element_str = "{:22}".format(str(coefficients[i, j]))
        #             print(element_str, end=" ", file=file)
        #     if i < len(solution):
        #         print("{:.10f}".format(constants[i]), file=file)
        # if abs(svch) < 0.000001:
        #     break

    # coefficients1 = Matrix([
    #     [1, 0, 0, 0, -2, 0, 0],
    #     [0, 1, 0, 0, 0, -2, 0],
    #     [0, 0, 1, 0, -2, -1, 0],
    #     [0, 0, 0, 1, -1, -1, 0],
    #     [2 * p_[0], 0, 2 * p_[2], p_[3], p_[4], 0, -(B_[0])],
    #     [0, 2 * p_[1], p_[2], p_[3], 0, p_[5], -(B_[1])],
    #     [p_[0], p_[1], p_[2], p_[3], p_[4], p_[5], 0]
    # ])
    coefficients1 = [[1, 0, 0, 0, -2, 0, 0],
                     [0, 1, 0, 0, 0, -2, 0],
                     [0, 0, 1, 0, -2, -1, 0],
                     [0, 0, 0, 1, -1, -1, 0],
                     [2 * p_[0], 0, 2 * p_[2], p_[3], p_[4], 0, -(B_[0])],
                     [0, 2 * p_[1], p_[2], p_[3], 0, p_[5], -(B_[1])],
                     [p_[0], p_[1], p_[2], p_[3], p_[4], p_[5], 0]]
    constants1 = Matrix([-k_kk[0], -k_kk[1], -k_kk[2], -k_kk[3], 0, 0, 0])
    solution1 = cramer_solver(Matrix(coefficients1), constants1)

    # coefficients2 = Matrix([
    #     [1, 0, 0, 0, -2, 0, 0],
    #     [0, 1, 0, 0, 0, -2, 0],
    #     [0, 0, 1, 0, -2, -1, 0],
    #     [0, 0, 0, 1, -1, -1, 0],
    #     [2 * p_[0], 0, 2 * p_[2], p_[3], p_[4], 0, -(B_[0])],
    #     [0, 2 * p_[1], p_[2], p_[3], 0, p_[5], -(B_[1])],
    #     [p_[0], p_[1], p_[2], p_[3], p_[4], p_[5], 0]
    # ])
    coefficients2 = [[1, 0, 0, 0, -2, 0, 0],
                     [0, 1, 0, 0, 0, -2, 0],
                     [0, 0, 1, 0, -2, -1, 0],
                     [0, 0, 0, 1, -1, -1, 0],
                     [2 * p_[0], 0, 2 * p_[2], p_[3], p_[4], 0, -(B_[0])],
                     [0, 2 * p_[1], p_[2], p_[3], 0, p_[5], -(B_[1])],
                     [p_[0], p_[1], p_[2], p_[3], p_[4], p_[5], 0]]
    constants2 = Matrix([0, 0, 0, 0, 0, 0, p_sum])
    solution2 = cramer_solver(Matrix(coefficients2), constants2)

    mut = 1.008 * 2 + k * 16
    sum5 = 0
    for j in range(6):
        sum5 += I_k[j] * p_[j]
    I = 4184 * sum5 / Mt / mut
    sum6 = 0
    for j in range(6):
        sum6 += S_k[j] * p_[j]
    S = 4184 * sum6 / Mt / mut
    mu = Mt * mut / p
    r = 8314 / mu
    alpha = (1 / t) * (1 - solution1[6])
    beta = (1 / p) * (solution2[6])
    sum7 = 0
    for j in range(6):
        sum7 += Cp_k[j] * p_[j]
    cpf = 4184 * sum7 / Mt / mut
    kf = 1 / (1 - (r / cpf))
    af = (t * r * kf) ** 0.5
    sum8 = 0
    for j in range(6):
        sum8 += I_k[j] * p_[j] * solution1[j]
    cpe = cpf + (1 / t) * (4184 * sum8 / Mt / mut - I * solution1[6])
    kr = cpe / (cpe - ((r * (1 - solution1[6]) ** 2) / solution2[6]))
    a = (t * r * kr / solution2[6]) ** 0.5

    for i in range(6):
        print('p', ko[i], '=', p_[i], 'кг/см2', file=file)
    print('%-6s %-5d' % ('Число итераций -', o), file=file)
    bo = '%9s %s %20f %s'
    print(bo % ('I', '=', I, 'Дж/кг/К'), file=file)
    print(bo % ('S', '=', S, 'Дж/кг/К'), file=file)
    print(bo % ('Mu', '=', mu, 'кг/кмоль'), file=file)
    print(bo % ('R', '=', r, 'Дж/кг/К'), file=file)
    print(bo % ('alpha_p', '=', alpha, '1/К'), file=file)
    print(bo % ('beta_t', '=', beta, 'м2/H'), file=file)
    print(bo % ('cp_f', '=', cpf, 'Дж/кг/К'), file=file)
    print('%9s %s %20f' % ('k_f', '=', kf), file=file)
    print(bo % ('a_f', '=', af, 'м/с'), file=file)
    print(bo % ('c_p_e', '=', cpe, 'Дж/кг/К'), file=file)
    print('%9s %s %20f' % ('k_R', '=', kr), file=file)
    print(bo % ('a', '=', a, 'м/с'), file=file)
    strochka()
